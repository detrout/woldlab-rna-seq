#
# Copyright 2021 California Institute of Technology
#
# Licensed under a BSD like license that includes a
# non-endorsement clause.
#

import codecs
import csv
import gzip
import hashlib
from io import BytesIO, StringIO
import netrc
import os
from pathlib import Path
import requests
import shutil
import stat
from subprocess import run
import tarfile
import time
from urllib.parse import urlparse

## Example configuration file.
## This configuration example is based on a 10x multiome experiment
## cDNA
## List all read1 and read2 fastq files in the same order.
#read1:
#  - ENCFF150FBF
#  - ENCFF385IAW
## barcode cell+UMI
#read2:
#  - ENCFF351VBS
#  - ENCFF503CCI
## include_intron tells star to include intronic genomic regions in
## gene counts
#include_intron: True
## the multiome RNA-seq experiments appear to be forward stranded
#stranded: "Forward"
## parameter to set umi length Chromium V2 uses 10, other versions use 12.
#umi_length: 12
## experiment id is added to archive file metadata
#experiment_accession: "ENCSR724KET"
## experiment description added to archive file metadata
#description: "snRNA on human adrenal gland. Please note that this experiment is part of 10x multiome and has a corresponding scATAC experiment (ENCSR693GAD)"
## library id is added to the archive file metadata and is intended to
## separate replicates
#library_accession: "ENCLB002DZK"
## id to represent the particular pipeline software version and index versions
#analysis_step_version: "diane-trout:2"
#
## we need to define a target directory where the genome index is going to be
## this allows using --singularity-args "-B <index dir>:$(pwd)/genome:ro" to
## share the index
#genome_dir: "genome"
#genome_index_url: "https://woldlab.caltech.edu/~diane/genome/GRCh38-V29-male-2.7.8a.tar.gz"
#allow_list_url: "https://woldlab.caltech.edu/~diane/genome/737K-arc-v1.txt.gz"
#star_container: "https://woldlab.caltech.edu/~diane/containers/star-dev_EoI_2.7.9a-2021-09-10.sif"
#
## These are guesses based on my current experiments and might need to change
#mem_mb: 65536
#disk_mb: 51200


MISSING_ERROR = "Please define before upload"
DEFAULT_10X_CB_LENGTH=16
DEFAULT_MEM_MB = 1000
SOLO_ROOT = Path("Solo.out")
MULTIREAD_NAME = {
    "Unique": "matrix.mtx",
    "Rescue": "UniqueAndMult-Rescue.mtx",
    "EM": "UniqueAndMult-EM.mtx",
}

configfile: "config.yaml"

for key in ['genome_dir']:
    if key in config and config[key].startswith('~'):
        config[key] = str(Path(config[key]).expanduser())


def compute_allow_list_name(url):
    parts = urlparse(url)
    filename = Path(parts.path).name
    return filename.replace(".gz", "")


config['allow_file'] = compute_allow_list_name(config['allow_list_url'])


try:
    authdb = netrc.netrc()
    username, _, password = authdb.hosts['www.encodeproject.org']
    auth = (username, password) if username is not None else None
except FileNotFoundError:
    auth = None


def generate_read_argument(config, read):
    argument = []
    read_num = int(read[-1])
    for accession in config[read]:
        argument.append(
            "{accession}_R{read_num}.fastq.gz".format(
                accession=accession,
                read_num=read_num,
            ))
    return argument

#######
# Functions for making filtered .mtx file

def read_barcode_lineno_map(filename):
    """Build a map of barcodes to line number from filename
    """
    barcodes = {}
    with open(filename, "rt") as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        for i, line in enumerate(reader):
            barcodes[line[0]] = i + 1

    return barcodes


def compute_raw_to_filtered_map(filtered_barcodes, raw_barcodes):
    """Generate a raw index to filtered index mapping"""
    raw_to_filtered_mapping = {}
    for filtered_barcode in filtered_barcodes:
        filtered_index = filtered_barcodes[filtered_barcode]
        raw_index = raw_barcodes[filtered_barcode]
        raw_to_filtered_mapping[raw_index] = filtered_index
    return raw_to_filtered_mapping


def filter_mtx(raw_barcode_filename, raw_matrix_filename, filtered_barcode_filename):
    """generating a filtered mtx file from raw and filtered values

    We read the the raw barcodes, raw matrix, and filtered barcodes
    and use that to generate a mapping between the raw barcode index
    values and the filtered barcode index values found in the mtx file.
    """
    raw_barcodes = read_barcode_lineno_map(raw_barcode_filename)
    filtered_barcodes = read_barcode_lineno_map(filtered_barcode_filename)
    raw_to_filtered_mapping = compute_raw_to_filtered_map(filtered_barcodes, raw_barcodes)

    header = True
    results = []
    with open(raw_matrix_filename, "rt") as instream:
        # copy comments
        for line in instream:
            if line.startswith("%"):
                yield line
            elif header:
                # After the comment comes the one header line
                total_rows, total_columns, total_counts = [int(x) for x in line.rstrip().split()]
                assert total_columns == len(raw_barcodes)
                header = False
            else:
                # row, column, count
                row, column, count = line.rstrip().split()
                row = int(row)
                column = int(column)
                if column in raw_to_filtered_mapping:
                    new_column = raw_to_filtered_mapping[column]
                    results.append((row, new_column, count))

        rs = sorted(results, key=lambda row: (row[1], row[0]))
        total_columns = len(filtered_barcodes)
        total_counts = len(rs)
        yield "{} {} {}\n".format(total_rows, total_columns, total_counts)
        for row, column, count in rs:
            yield "{} {} {}\n".format(row, column, count)


#######
# Functions for making manifest
#


def compute_md5sums(filenames):
    BLOCK = 2 ** 20
    results = []
    for f in filenames:
        with open(f, "rb") as instream:
            md5 = hashlib.md5()
            readable = instream.readable()
            while readable:
                read_block = instream.read(BLOCK)
                if len(read_block) == 0:
                    readable = False
                else:
                    md5.update(read_block)
        results.append((f, md5.hexdigest()))
    return results


def create_metadata(config, output_type, md5s):
    metadata = {
        "type": "MatrixMarketGeneArchive_v1",
        "output_type": output_type,
        "experiment_accession": config.get("experiment_accession", MISSING_ERROR),
        "description": config.get("description", MISSING_ERROR),
        "library_accession": config.get("library_accession", MISSING_ERROR),
        "analysis_step_version": config.get("analysis_step_version", MISSING_ERROR),
    }

    for filename, md5 in md5s:
        metadata[filename] = "md5sum:{}".format(md5)

    return metadata


def write_metadata(outstream, config):
    writer = csv.writer(outstream, delimiter="\t")
    writer.writerow(["name", "value"])
    for key in config:
        writer.writerow([key, config[key]])
    return outstream


####
# functions for making archive file
def make_list_of_archive_files(solo_root, quantification="GeneFull", multiread="Unique", matrix="raw"):
    archive_files = []

    archive_files.append(solo_root / quantification / matrix / "barcodes.tsv")
    archive_files.append(solo_root / quantification / matrix / "features.tsv")

    if quantification == "SJ":
        if multiread != "Unique":
            raise ValueError("Splice junctions do not support multread assignment")
        if matrix != "raw":
            raise ValueError("Splice junctions are only available as raw")

    archive_files.append(solo_root / quantification / matrix / MULTIREAD_NAME[multiread])
    return archive_files


def update_tarinfo(info, filename):
    stat_info = os.stat(filename)
    info.size = stat_info[stat.ST_SIZE]
    info.mode = stat_info[stat.ST_MODE]
    info.mtime = time.time()
    info.uid = stat_info[stat.ST_UID]
    info.gid = stat_info[stat.ST_GID]
    info.type = tarfile.REGTYPE


def make_output_type_term(quantification="GeneFull", multiread="Unique", matrix="raw"):
    assert quantification in ["Gene", "GeneFull", "GeneFull_Ex50pAS", "SJ"]
    assert multiread in ["Unique", "EM"]
    assert matrix in ["filtered", "raw"]

    gene_term = {
        "Gene": "gene",
        "GeneFull": "gene",
        "GeneFull_Ex50pAS": "gene",
        "SJ": "splice junction",
    }[quantification]

    multiread_term = {
        "Unique": "unique",
        "EM": "EM",
    }[multiread]

    matrix_term = matrix

    output_type = "sparse {multiread} {quantification} {count_matrix}".format(
        multiread=multiread_term,
        quantification=gene_term,
        count_matrix=matrix_term,
    )
    return output_type


def archive_star(solo_root, quantification="GeneFull", multiread="Unique", matrix="raw"):
    assert quantification in ["Gene", "GeneFull", "GeneFull_Ex50pAS", "SJ"]
    assert multiread in ["Unique", "Rescue", "EM"]
    assert matrix in ["filtered", "raw"]

    archive_files = make_list_of_archive_files(solo_root, quantification, multiread, matrix)
    output_type = make_output_type_term(quantification, multiread, matrix)

    md5s = compute_md5sums(archive_files)
    manifest = create_metadata(config, output_type, md5s)
    manifest_buffer = BytesIO(write_metadata(StringIO(), manifest).getvalue().encode("utf-8"))

    tar_name = "{}_{}_{}.tar.gz".format(quantification, multiread, matrix)
    with tarfile.open(tar_name, "w:gz") as archive:
        info = tarfile.TarInfo("manifest.tsv")
        update_tarinfo(info, archive_files[0])
        info.size = len(manifest_buffer.getvalue())
        archive.addfile(info, manifest_buffer)
        for filename in archive_files:
            info = tarfile.TarInfo(str(filename.relative_to(solo_root)))
            update_tarinfo(info, filename)
            with open(filename, "rb") as instream:
                archive.addfile(info, instream)


###
# Main rules
def get_gene_model(config):
    return "GeneFull_Ex50pAS" if config['include_intron'] else "Gene"


rule ALL:
    input:
        "Log.final.out",
        "{}_Unique_filtered.tar.gz".format(get_gene_model(config)),
        "{}_EM_filtered.tar.gz".format(get_gene_model(config)),
        "{}_Unique_raw.tar.gz".format(get_gene_model(config)),
        "{}_EM_raw.tar.gz".format(get_gene_model(config)),
        "SJ_Unique_raw.tar.gz",


rule get_encode_fastq:
    output:
        "{accession}_R{read}.fastq.gz"
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB,
    run:
        url = "https://www.encodeproject.org/files/{accession}/@@download/{accession}.fastq.gz".format(
            accession=wildcards.accession)
        with requests.get(url, auth=auth, stream=True) as instream:
            instream.raise_for_status()
            with open(output[0], "wb") as outstream:
                shutil.copyfileobj(instream.raw, outstream)


rule genome:
    output:
        genome_dir = temp(directory(config['genome_dir']))
    resources:
        mem_mb = DEFAULT_MEM_MB,
    threads: 1
    run:
        genome_dir = Path(output.genome_dir)
        print("genome_dir")
        if not Path(Path(genome_dir).parts[0]).exists():
            print("Creating ", genome_dir.parts[0])
            os.mkdir(genome_dir.parts[0])
        if "genome_index_url" in config:
            print("Downloading? {}".format(config["genome_index_url"]))
            with requests.get(config["genome_index_url"], stream=True) as instream:
                tar = tarfile.open(fileobj=instream.raw, mode="r:*")
                # this is a vulnerability, only use trusted tar files
                # See https://docs.python.org/3/library/tarfile.html#tarfile.TarFile.extractall
                tar.extractall(path=genome_dir.parts[0])


rule allow_list:
    output:
        allow_file = config['allow_file']
    params:
        allow_list_url = config['allow_list_url']
    resources:
        mem_mb = DEFAULT_MEM_MB,
    run:
        with requests.get(params.allow_list_url, stream=True) as instream:
            instream.raise_for_status()
            with open(output.allow_file, "wb") as outstream:
                if (params.allow_list_url.endswith('.gz') or
                   instream.headers.get("Content-Encoding") == "gzip"):
                    instream = gzip.GzipFile(fileobj=instream.raw)
                else:
                    instream = instream.raw
                shutil.copyfileobj(instream, outstream)


rule star_solo_10x:
    input:
        sequence_reads = generate_read_argument(config, "read2"),
        barcode_reads = generate_read_argument(config, "read1"),
        genome_index = config['genome_dir'],
        allow_list = config['allow_file'],
    params:
        stranded = config['stranded'],
        gene_model = get_gene_model(config),
        sequence_reads = ",".join(generate_read_argument(config, "read2")),
        barcode_reads = ",".join(generate_read_argument(config, "read1")),
        umi_length = int(config["umi_length"]),
        cb_length = int(config.get("cb_length", DEFAULT_10X_CB_LENGTH)),
        star_tmp = temp(directory("_STARtmp")),
    resources:
        mem_mb = config['mem_mb'],
        mem_bytes = config['mem_mb'] * (2 ** 20),
        disk_mb = config['disk_mb'],
    threads: 16
    log: "star_solo_10x.out"
    output:
        aligned_bam = "Aligned.sortedByCoord.out.bam",
        log_final = "Log.final.out",
        log_progress = "Log.progress.out",
        splice_junctions = "SJ.out.tab",
        barcode_stats = SOLO_ROOT / "Barcodes.stats",
        features_stats = SOLO_ROOT / get_gene_model(config) / "Features.stats",
        umis = SOLO_ROOT / get_gene_model(config) / "UMIperCellSorted.txt",
        filtered_barcodes = SOLO_ROOT / get_gene_model(config) / "filtered" / "barcodes.tsv",
        filtered_features = SOLO_ROOT / get_gene_model(config) / "filtered" / "features.tsv",
        filtered_uniq_matrix = SOLO_ROOT / get_gene_model(config) / "filtered" / "matrix.mtx",
        raw_barcodes = SOLO_ROOT / get_gene_model(config) / "raw" / "barcodes.tsv",
        raw_features = SOLO_ROOT / get_gene_model(config) / "raw" / "features.tsv",
        raw_unique_matrix = SOLO_ROOT / get_gene_model(config) / "raw" / "matrix.mtx",
        raw_em_matrix = SOLO_ROOT / get_gene_model(config) / "raw" / "UniqueAndMult-EM.mtx",
        sj_barcodes = SOLO_ROOT / "SJ" / "raw" / "barcodes.tsv",
        sj_features = SOLO_ROOT / "SJ" / "raw" / "features.tsv",
        sj_matrix = SOLO_ROOT / "SJ" / "raw" / "matrix.mtx",
    singularity:
        config['star_container']
    shell:
        "STAR --genomeDir {input.genome_index} \
           --readFilesIn {params.sequence_reads} {params.barcode_reads} \
           --readFilesCommand zcat \
           --runThreadN {threads} \
           --genomeLoad NoSharedMemory \
           --outFilterMultimapNmax 20 \
           --alignSJoverhangMin 8 \
           --alignSJDBoverhangMin 1 \
           --outFilterMismatchNmax 999 \
           --outFilterMismatchNoverReadLmax 0.04 \
           --alignIntronMin 20 \
           --alignIntronMax 1000000 \
           --alignMatesGapMax 1000000 \
           --outSAMheaderCommentFile COfile.txt \
           --outSAMheaderHD @HD VN:1.4 SO:coordinate \
           --outSAMunmapped Within \
           --outFilterType BySJout \
           --outSAMattributes NH HI AS NM MD CB CR CY UB UR UY gx gn \
           --outSAMstrandField intronMotif \
           --outSAMtype BAM SortedByCoordinate \
           --sjdbScore 1 \
           --clipAdapterType CellRanger4 \
           --outFilterScoreMin 30 \
           --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
           --soloUMIdedup 1MM_CR \
           --soloUMIfiltering MultiGeneUMI_CR \
           --soloType CB_UMI_Simple \
           --soloCellFilter EmptyDrops_CR \
           --soloUMIlen {params.umi_length} \
           --soloCBlen {params.cb_length} \
           --soloBarcodeReadLength 0 \
           --soloCBwhitelist {input.allow_list} \
           --soloStrand {params.stranded} \
           --soloFeatures {params.gene_model} SJ \
           --soloMultiMappers Unique EM \
           --limitBAMsortRAM {resources.mem_bytes} \
           --outTmpDir {params.star_tmp} \
           --outFileNamePrefix ./ 2>&1 >> {log} \
"

rule filter_em_matrix:
    input:
        filtered_barcode_tsv = SOLO_ROOT / get_gene_model(config) / "filtered" / "barcodes.tsv",
        raw_barcode_tsv = SOLO_ROOT / get_gene_model(config) / "raw" / "barcodes.tsv",
        raw_em_matrix_mtx = SOLO_ROOT / get_gene_model(config) / "raw" / "UniqueAndMult-EM.mtx",
    output:
        filtered_em_mtx = SOLO_ROOT / get_gene_model(config) / "filtered" / "UniqueAndMult-EM.mtx",
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    run:
        with open(output.filtered_em_mtx, "wt") as outstream:
            for line in filter_mtx(input.raw_barcode_tsv, input.raw_em_matrix_mtx, input.filtered_barcode_tsv):
                outstream.write(line)


rule to_archive:
    input:
        # At least snakemake 5.4.0 wants input files to be a list of strings
        # make_list_of_archive_files returns a list of Paths so we need to
        # convert them
        lambda wildcards: [str(x) for x in make_list_of_archive_files(
            SOLO_ROOT, wildcards.gene_model, wildcards.multiread, wildcards.matrix)]
    output:
        "{gene_model}_{multiread}_{matrix}.tar.gz"
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB,
    run:
        archive_star(SOLO_ROOT, wildcards.gene_model, wildcards.multiread, wildcards.matrix)
