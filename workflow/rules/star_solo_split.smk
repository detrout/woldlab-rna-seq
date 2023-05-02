import gzip
import os
import requests
import shutil
import tarfile
from urllib.parse import urlparse
from encoded_client.encoded import ENCODED

from woldrnaseq.snakeutils import (
    compute_inclusion_list_name,
)

DEFAULT_MEM_MB = 1000


rule download_inclusion_list:
    output:
        allow_file = directory(temp(compute_inclusion_list_name(config['inclusion_list_url'])))
    params:
        inclusion_list_url = config['inclusion_list_url']
    resources:
        mem_mb = DEFAULT_MEM_MB,
    run:
        server = ENCODED("www.encodeproject.org")
        with server.get_response(params.inclusion_list_url, stream=True) as response:
            encoding = response.headers.get("Content-Encoding")
            if (params.inclusion_list_url.endswith(".tar.gz")):
                with tarfile.open(mode="r:gz", fileobj=response.raw) as tar:
                    tar.extractall(".")
            elif (params.inclusion_list_url.endswith('.gz') or encoding == "gzip"):
                with open(output.allow_file, "wb") as outstream:
                    decompress = gzip.GzipFile(fileobj=response.raw)
                    shutil.copyfileobj(decompress, outstream)
            else:
                with open(output.allow_file, "wb") as outstream:
                    shutil.copyfileobj(response.raw, outstream)


def generate_star_read_argument(config, read):
    """Convert accessions into list of {accession}_R{read_num}.fastq.gz

    this makes formatting for the star paired readFilesIn argument a bit
    simpler.
    """
    argument = []
    read_num = int(read[-1])
    for accession in config[read]:
        argument.append(
            "{accession}_R{read_num}.fastq.gz".format(
                accession=accession,
                read_num=read_num,
            ))
    return argument


rule star_solo_splitseq:
    input:
        sequence_reads = generate_star_read_argument(config, "read2"),
        barcode_reads = generate_star_read_argument(config, "read1"),
        genome_index = config['genome_dir'],
        inclusion_list = compute_inclusion_list_name(config['inclusion_list_url']),
        cofile = "COfile.txt",
    params:
        stranded = config['stranded'],
        gene_model = get_gene_model(),
        sequence_reads = ",".join(generate_star_read_argument(config, "read2")),
        barcode_reads = ",".join(generate_star_read_argument(config, "read1")),
        star_tmp = temp(directory("_STARtmp")),
    resources:
        mem_mb = config['mem_mb'],
        mem_bytes = config['mem_mb'] * (2 ** 20),
        disk_mb = config['disk_mb'],
    threads: 16
    log: "star_solo_splitseq.out"
    output:
        aligned_bam = "Aligned.sortedByCoord.out.bam",
        log_final = "Log.final.out",
        log_progress = "Log.progress.out",
        log_out = "Log.out",
        splice_junctions = "SJ.out.tab",
        barcode_stats = SOLO_ROOT / "Barcodes.stats",
        features_stats = SOLO_ROOT / get_gene_model() / "Features.stats",
        umis = SOLO_ROOT / get_gene_model() / "UMIperCellSorted.txt",
        gene_summary = SOLO_ROOT / get_gene_model() / "Summary.csv",
        raw_barcodes = temp(SOLO_ROOT / get_gene_model() / "raw_2bc" / "barcodes.tsv"),
        raw_features = temp(SOLO_ROOT / get_gene_model() / "raw_2bc" / "features.tsv"),
        raw_unique_matrix = temp(SOLO_ROOT / get_gene_model() / "raw_2bc" / "matrix.mtx"),
        raw_em_matrix = temp(SOLO_ROOT / get_gene_model() / "raw_2bc" / "UniqueAndMult-EM.mtx"),
        sj_feature_stats = SOLO_ROOT / "SJ" / "Features.stats",
        sj_summary = SOLO_ROOT / "SJ" / "Summary.csv",
        sj_barcodes = temp(SOLO_ROOT / "SJ" / "raw_2bc" / "barcodes.tsv"),
        sj_features = temp(SOLO_ROOT / "SJ" / "raw_2bc" / "features.tsv"),
        sj_matrix = temp(SOLO_ROOT / "SJ" / "raw_2bc" / "matrix.mtx"),
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
           --clip5pAdapterSeq AAGCAGTGGTATCAACGCAGAGTGAATGGG \
           --outFilterScoreMin 30 \
           --soloUMIdedup 1MM_CR \
           --soloUMIfiltering MultiGeneUMI_CR \
           --soloType CB_UMI_Complex \
           --soloCellFilter EmptyDrops_CR \
           --soloCBmatchWLtype EditDist_2 \
           --soloCBposition 0_10_0_17 0_48_0_55 0_78_0_85 \
           --soloUMIposition 0_0_0_9 \
           --soloCBwhitelist {input.inclusion_list}/CB23.txt {input.inclusion_list}/CB23.txt {input.inclusion_list}/CB1.txt \
           --soloStrand {params.stranded} \
           --soloFeatures {params.gene_model} SJ \
           --soloMultiMappers Unique EM \
           --limitBAMsortRAM {resources.mem_bytes} \
           --outTmpDir {params.star_tmp} \
           --outFileNamePrefix ./ 2>&1 >> {log} ; \
        rm -rf {params.star_tmp} ; \
        rm -rfv Solo.out/{params.gene_model}/filtered/ ; \
        mv -v Solo.out/{params.gene_model}/raw/* Solo.out/{params.gene_model}/raw_2bc/ ; \
        mv -v Solo.out/SJ/raw/* Solo.out/SJ/raw_2bc/ ; \
"
# I really have no idea where the raw_2bc directory above is being created.


rule merge_raw_cells:
    input:
        raw_barcodes = SOLO_ROOT / get_gene_model() / "raw_2bc" / "barcodes.tsv",
        raw_features = SOLO_ROOT / get_gene_model() / "raw_2bc" / "features.tsv",
        raw_unique_matrix = SOLO_ROOT / get_gene_model() / "raw_2bc" / "matrix.mtx",
        raw_em_matrix = SOLO_ROOT / get_gene_model() / "raw_2bc" / "UniqueAndMult-EM.mtx",
    output:
        raw_barcodes = temp(SOLO_ROOT / get_gene_model() / "raw" / "barcodes.tsv"),
        raw_features = temp(SOLO_ROOT / get_gene_model() / "raw" / "features.tsv"),
        raw_unique_matrix = temp(SOLO_ROOT / get_gene_model() / "raw" / "matrix.mtx"),
        raw_em_matrix = temp(SOLO_ROOT / get_gene_model() / "raw" / "UniqueAndMult-EM.mtx"),
    resources:
        mem_mb = 8192
    threads: 1
    run:
        from woldrnaseq.splitseq_merger import write_merged_splitseq_matrix
        destination_dir = Path(output.raw_unique_matrix).parent
        write_merged_splitseq_matrix(input.raw_unique_matrix, destination_dir)
        write_merged_splitseq_matrix(input.raw_em_matrix, destination_dir)


rule merge_sj_cells:
    input:
        sj_barcodes = SOLO_ROOT / "SJ" / "raw_2bc" / "barcodes.tsv",
        sj_features = SOLO_ROOT / "SJ" / "raw_2bc" / "features.tsv",
        sj_matrix = SOLO_ROOT / "SJ" / "raw_2bc" / "matrix.mtx",
    output:
        sj_barcodes = temp(SOLO_ROOT / "SJ" / "raw" / "barcodes.tsv"),
        sj_features = temp(SOLO_ROOT / "SJ" / "raw" / "features.tsv"),
        sj_matrix = temp(SOLO_ROOT / "SJ" / "raw" / "matrix.mtx"),
    resources:
        mem_mb = 8192
    threads: 1
    run:
        from woldrnaseq.splitseq_merger import write_merged_splitseq_matrix
        sj_dir = Path(output.sj_matrix).parent
        write_merged_splitseq_matrix(input.sj_matrix, sj_dir)


rule filter_merged_unique_cells:
    input:
        raw_barcodes = SOLO_ROOT / get_gene_model() / "raw" / "barcodes.tsv",
        raw_features = SOLO_ROOT / get_gene_model() / "raw" / "features.tsv",
        raw_unique_matrix = SOLO_ROOT / get_gene_model() / "raw" / "matrix.mtx",
    output:
        filtered_barcodes = temp(SOLO_ROOT / get_gene_model() / "filtered" / "barcodes.tsv"),
        filtered_features = temp(SOLO_ROOT / get_gene_model() / "filtered" / "features.tsv"),
        filtered_unique_matrix = temp(SOLO_ROOT / get_gene_model() / "filtered" / "matrix.mtx"),
    params:
        input_directory = lambda wildcards, input: Path(input.raw_unique_matrix).parent,
        output_directory = lambda wildcards, output: Path(output.filtered_unique_matrix).parent,
        star_tmp = temp(directory("_STARtmp")),
    resources:
        mem_mb = 8192
    threads: 1
    log: "filter_merged_unique_cells.out"
    singularity:
        config['star_container']
    shell:
        "STAR --runMode soloCellFiltering \
            {params.input_directory}  \
            {params.output_directory}/ \
            --soloCellFilter EmptyDrops_CR ; \
        rm -rf {params.star_tmp}"
