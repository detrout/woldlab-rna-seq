#
# Copyright 2021 California Institute of Technology
#
# Licensed under a BSD like license that includes a
# non-endorsement clause.
#
# To run use something like:
# snakemake -j 1 --cores 8 --use-conda \
#    --snakefile ${woldlab-rna-seq}/workflow/process-encode-split-subpool.snakefile
# adjusting the number of jobs and cores to reflect your environment.


import os
import datetime
from pathlib import Path
import pandas

from woldrnaseq.snakeutils import (
    compute_inclusion_list_name,
)
from mex_gene_archive.starsolo import (
    archive_star_solo,
    MULTIREAD_NAME,
)

from snakemake.utils import min_version
min_version("6.0")

## Example configuration file.
## This configuration example is for processing multiple subpools
## of a split-seq experiment and then merging them into one result.
##
##
# library:
##   Repeat subpool/sublibrary ids for each subpool that needs
##   to be processed and merged.
#    {library_id}:
##       List all read1 and read2 fastq files in the same order.
##       Read 1 should contain the real sequence
#        read1:
#          - ENCFF150FBF
#          - ENCFF385IAW
##       Read 2 should contain the barcode cell+UMI
#        read2:
#          - ENCFF351VBS
#          - ENCFF503CCI
## include_intron tells star to include intronic genomic regions in
## gene counts, for single nucleus experiments this should be true.
#include_intron: True
## The split-seq experiments are forward stranded.
#stranded: "Forward"
## biosample accession id
# biosample: "biosampe_accession"
## experiment id is added to archive file metadata
#experiment_accession: "ENCSR724KET"
## experiment description added to archive file metadata
#description: "snRNA on human adrenal gland. Please note that this experiment is part of 10x multiome and has a corresponding scATAC experiment (ENCSR693GAD)"
## id to represent the particular pipeline software version and index versions
#analysis_step_version: "diane-trout:2"
#
## we need to define a target directory where the genome index is going to be
## this allows using --singularity-args "-B <index dir>:$(pwd)/genome:ro" to
## share the index
## This is now the default value. If you want you can override, but you don't need to.
#genome_dir: "genome"
## If you provide an accession for the genome index the code will attempt to download
## the genome index from that location (if needed)
#genome_accession: "TSTFF984371"
## If you only provide a genome_index_url, the code will attempt to use the filename
## portion to determine the accession
#genome_index_url: "https://woldlab.caltech.edu/~diane/genome/GRCh38-V29-male-2.7.8a.tar.gz"
## If the genome_accession is set the code will pull the assembly and genome_annotation
## from the ENCODE portal
#assembly: "GRCh38"
#genome_annotation: "V29"
## Given the accession it will determine the inclusion_list_url
##inclusion_accession: "encode:parse-biosciences-v1"
## if there is no inclusion_accession, it will guess the accession from the inclusion_list_url
#inclusion_list_url: "https://woldlab.caltech.edu/~diane/genome/737K-arc-v1.txt.gz"
## The pipeline now defaults to the container used for my analysis (The version of STAR after 2.7.9a)
# You'll probably get better performance if you put a copy somewhere local
#star_container: "https://woldlab.caltech.edu/~diane/containers/star-dev_EoI_2.7.9a-2021-09-10.sif"
## This is a container with python and scanpy in it. It also defaults to my container.
# You'll probably get better performance if you put a copy somewhere local
#scanpy_container: "https://woldlab.caltech.edu/~diane/containers/bullseye-scanpy-1.8.2.sif"

#
## These are guesses based on my current experiments and might need to change
#mem_mb: 65536
#disk_mb: 51200

configfile: "config.yaml"

include:
    "rules/star_solo_config.smk"

# Set defaults for this version of the ENCODE scRNA-seq
# pipeline
config.setdefault(
    "alignment_step_run",
    "barbara-wold:starsolo-split-alignment-step-run"
)
config.setdefault(
    "quantification_step_run",
    "barbara-wold:starsolo-split-quantification-step-run",
)
config.setdefault(
    "merge_quantification_matrix_step_run",
    "barbara-wold:starsolo_merging_step_v0_run"
)


def get_gene_model(config):
    return "GeneFull_Ex50pAS" if config['include_intron'] else "Gene"


def get_submit_host(config):
    default_host = "test.encodedcc.org" #"www.encodeproject.org"
    return config.get("encode_portal_host", default_host)


def generate_star_read_argument(config, library_id, read):
    """Convert accessions into list of {accession}_R{read_num}.fastq.gz

    this makes formatting for the star paired readFilesIn argument a bit
    simpler.
    """
    argument = []
    read_num = int(read[-1])

    library = config["library"].get(library_id)
    if library is None:
        raise ValueError("Unable to find {} in {}".format(library_id, config["library"].keys()))
    
    if read not in library:
        raise ValueError("Read {} not found under library {}".format(read, libary_id))

    for accession in library[read]:
        argument.append(
            "{library_id}/{accession}_R{read_num}.fastq.gz".format(
                library_id=library_id,
                accession=accession,
                read_num=read_num,
            ))
    return argument


def list_split_subpool_targets(config):
    # Add in merged targets
    default = [
        "{}_Unique_raw_merged.tar.gz".format(get_gene_model(config)),
        "{}_EM_raw_merged.tar.gz".format(get_gene_model(config)),        
        "metadata.{}.csv".format(get_submit_host(config)),
    ]
    if config["automatic_submission"]:
        default.extend([
            "posted.{}.csv".format(get_submit_host(config)),
        ])

    # Add in per library targets?
    for library_id in config["library"]:
        library = config["library"][library_id]
        library_dir = Path(library_id)

        default.extend([
            library_dir / "Log.final.out",
            library_dir / "Aligned.sortedByCoord.out.bam",
            library_dir / "{}_Unique_raw.tar.gz".format(get_gene_model(config)),
            library_dir / "{}_EM_raw.tar.gz".format(get_gene_model(config)),
            library_dir / "SJ_Unique_raw.tar.gz",
            library_dir / "umi_per_cell.png",
            library_dir / "qc_metric_violin.{}_EM_raw.png".format(get_gene_model(config)),
            library_dir / "pct_count_mt.{}_EM_raw.png".format(get_gene_model(config)),
            library_dir / "n_genes_by_counts.{}_EM_raw.png".format(get_gene_model(config)),
            library_dir / "metadata.{}.csv".format(get_submit_host(config)),
        ])
        if config["automatic_submission"]:
            default.extend([
                library_dir / "posted.{}.csv".format(get_submit_host(config)),
                library_dir / "Log.final.out.{}.qc-upload".format(get_submit_host(config)),
                library_dir / SOLO_ROOT / get_gene_model(config) / "Summary.csv.{}.qc-upload".format(
                    get_submit_host(config)
                ),
                library_dir / "pct_count_mt.{}_EM_raw.png.{}.qc-upload".format(
                    get_gene_model(config), get_submit_host(config)
                ),
            ])
    return default

rule ALL:
    input:
        list_split_subpool_targets(config)


include:
    "rules/provide_star_index.smk"

include:
    "rules/star_co_file.smk"

rule get_fastq:
    output:
        "{library_id}/{accession}_R{read}.fastq.gz"        
    resources:
        mem_mb = 1000
    threads: 1
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/getfastqs"


rule download_inclusion_list:
    output:
        allow_file = directory(temp(compute_inclusion_list_name(config['inclusion_list_url'])))
    params:
        inclusion_list_url = config['inclusion_list_url']
    resources:
        mem_mb = DEFAULT_MEM_MB,
    threads: 1
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

        
rule star_solo_splitpool:
    input:
        sequence_reads = lambda wildcards: generate_star_read_argument(config, wildcards.library_id, "read2"),
        barcode_reads = lambda wildcards: generate_star_read_argument(config, wildcards.library_id, "read1"),
        genome_index = config['genome_dir'],
        inclusion_list = compute_inclusion_list_name(config['inclusion_list_url']),
    params:
        stranded = config['stranded'],
        gene_model = get_gene_model(config),
        library_id = lambda wildcards: wildcards.library_id,
        sequence_reads = lambda wildcards: ",".join(generate_star_read_argument(config, wildcards.library_id, "read2")),
        barcode_reads = lambda wildcards: ",".join(generate_star_read_argument(config, wildcards.library_id, "read1")),
    resources:
        mem_mb = config['mem_mb'],
        mem_bytes = config['mem_mb'] * (2 ** 20),
        disk_mb = config['disk_mb'],
    threads: 16
    #log: "star_solo_splitseq.out"
    output:
        aligned_bam = "{library_id}/Aligned.sortedByCoord.out.bam",
        log_final = "{library_id}/Log.final.out",
        log_progress = "{library_id}/Log.progress.out",
        log_out = "{library_id}/Log.out",
        splice_junctions = "{library_id}/SJ.out.tab",
        cofile = temp("{library_id}/COfile.txt"),
        barcode_stats = "{library_id}/Solo.out/Barcodes.stats",
        features_stats = expand(
            "{library_id}/Solo.out/{gene_model}/Features.stats",
            gene_model=get_gene_model(config),
            allow_missing=True),
        umis = expand(
            "{library_id}/Solo.out/{gene_model}/UMIperCellSorted.txt",
            gene_model=get_gene_model(config),
            allow_missing=True),
        gene_summary = expand(
            "{library_id}/Solo.out/{gene_model}/Summary.csv",
            gene_model=get_gene_model(config),
            allow_missing=True),
        raw_barcodes = temp(expand(
            "{library_id}/Solo.out/{gene_model}/raw/barcodes.tsv",
            gene_model=get_gene_model(config),
            allow_missing=True)),
        raw_features = temp(expand(
            "{library_id}/Solo.out/{gene_model}/raw/features.tsv",
            gene_model=get_gene_model(config),
            allow_missing=True)),
        raw_unique_matrix = temp(expand(
            "{library_id}/Solo.out/{gene_model}/raw/matrix.mtx",
            gene_model=get_gene_model(config),
            allow_missing=True)),
        raw_em_matrix = temp(expand(
            "{{library_id}}/Solo.out/{gene_model}/raw/UniqueAndMult-EM.mtx",
            gene_model=get_gene_model(config),
            allow_missing=True)),
        sj_feature_stats = "{library_id}/Solo.out/SJ/Features.stats",
        sj_summary = "{library_id}/Solo.out/SJ/Summary.csv",
        sj_barcodes = temp("{library_id}/Solo.out/SJ/raw/barcodes.tsv"),
        sj_features = temp("{library_id}/Solo.out/SJ/raw/features.tsv"),
        sj_matrix = temp("{library_id}/Solo.out/SJ/raw/matrix.mtx"),
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/star_solo_splitseq"


rule to_archive:
    #
    # This rule generates the various mex archive files
    input:
        barcodes = "{library_id}/Solo.out/{gene_model}/{matrix}/barcodes.tsv",
        features = "{library_id}/Solo.out/{gene_model}/{matrix}/features.tsv",
        matrix = lambda wildcards: "{library_id}/Solo.out/{gene_model}/{matrix}/" + MULTIREAD_NAME[wildcards.multiread],
    output:
        "{library_id}/{gene_model}_{multiread}_{matrix}.tar.gz"
    wildcard_constraints:
        gene_model = "(GeneFull_Ex50pAS)|(GeneFull)|(Gene)|(SJ)",
        multiread = "(EM)|(Unique)",
        matrix = "(raw)|(filtered)",
    threads: 1
    resources:
        mem_mb = 256,
    run:
        temp_config = config.copy()
        temp_config["library_accession"] = wildcards.library_id
        archive_star_solo(
            Path(wildcards.library_id) / "Solo.out",
            temp_config,
            wildcards.gene_model,
            wildcards.multiread,
            wildcards.matrix
        )


rule generate_umi:
    input:
        expand("{library_id}/Solo.out/{gene_model}/UMIperCellSorted.txt",
               gene_model=get_gene_model(config),
               allow_missing=True),
    output:
        "{library_id}/umi_per_cell.png"
    resources:
        mem_mb = DEFAULT_MEM_MB
    run:
        from woldrnaseq.plots.star_solo_barcodes import main
        main(["--umi-per-cell", input[0], "-o", output[0]])


rule generate_pool_single_cell_mex_qc_plots:
    input:
        count_matrix = "{library_id}/{gene_model}_{multiread}_{matrix}.tar.gz",
        genome_dir = config['genome_dir'],
    output:
        qc_violin = "{library_id}/qc_metric_violin.{gene_model}_{multiread}_{matrix}.png",
        pct_mt = "{library_id}/pct_count_mt.{gene_model}_{multiread}_{matrix}.png",
        genes_by_count = "{library_id}/n_genes_by_counts.{gene_model}_{multiread}_{matrix}.png",
    threads: 1
    resources:
        mem_mb = 1000
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/single_cell_mex_qc_plots"

include:
    "rules/create_md5s.smk"


rule prepare_encode_subpool_metadata:
    input:
        bam = "{library_id}/Aligned.sortedByCoord.out.bam",
        bam_md5 = "{library_id}/Aligned.sortedByCoord.out.bam.md5",
        gene_unique_raw = expand("{library_id}/{gene_model}_Unique_raw.tar.gz",
                                 gene_model=get_gene_model(config),
                                 allow_missing=True),
        gene_unique_raw_md5 = expand("{library_id}/{gene_model}_Unique_raw.tar.gz.md5",
                                     gene_model=get_gene_model(config),
                                     allow_missing=True),
        gene_multi_raw = expand("{library_id}/{gene_model}_EM_raw.tar.gz",
                                gene_model=get_gene_model(config),
                                allow_missing=True),
        gene_multi_raw_md5 = expand("{library_id}/{gene_model}_EM_raw.tar.gz.md5",
                                    gene_model=get_gene_model(config),
                                    allow_missing=True),
        sj_unique_raw = "{library_id}/SJ_Unique_raw.tar.gz",
        sj_unique_raw_md5 = "{library_id}/SJ_Unique_raw.tar.gz.md5",
    output:
        expand("{library_id}/metadata.{host}.csv",
               host=get_submit_host(config),
               allow_missing=True)
    threads: 1
    resources:
        mem_mb = 1000
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/prepare_encode_subpool_metadata"

rule submit_star_solo_subpool_data:
    input:
        bam = "{library_id}/Aligned.sortedByCoord.out.bam",
        gene_unique_raw = expand("{library_id}/{gene_model}_Unique_raw.tar.gz",
                                 gene_model=get_gene_model(config),
                                 allow_missing=True),
        gene_multi_raw = expand("{library_id}/{gene_model}_EM_raw.tar.gz",
                                gene_model=get_gene_model(config),
                                allow_missing=True),
        sj_unique_raw = "{library_id}/SJ_Unique_raw.tar.gz",
        metadata = expand("{library_id}/metadata.{host}.csv",
                          host=get_submit_host(config),
                          allow_missing=True),
    output:
        # these names need to match what's computed by encoded_client.submission.make_upload_filename
        bam = "{{library_id}}_Aligned.sortedByCoord.out.bam.{host}.upload".format(
            host=get_submit_host(config)),
        gene_unique_raw = "{{library_id}}_{gene_model}_Unique_raw.tar.gz.{host}.upload".format(
            gene_model=get_gene_model(config),
            host=get_submit_host(config)),
        gene_multi_raw = "{{library_id}}_{gene_model}_EM_raw.tar.gz.{host}.upload".format(
            gene_model=get_gene_model(config),
            host=get_submit_host(config)),
        sj_unique_raw = "{{library_id}}_SJ_Unique_raw.tar.gz.{host}.upload".format(
            host=get_submit_host(config)),
        posted = expand("{library_id}/posted.{host}.csv",
                        host=get_submit_host(config),
                        allow_missing=True),
    params:
        submit_host = get_submit_host(config),
    log: "{library_id}/submit_star_solo_subpool_data.log"
    threads: 1
    resources:
        mem_mb = 1000
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/post_encode_star_solo_files"
        

rule post_star_qc:
    input:
        log_final = "{library_id}/Log.final.out",
        posted = expand("{library_id}/posted.{host}.csv",
                        host=get_submit_host(config),
                                 allow_missing=True),
    output:
        expand("{library_id}/Log.final.out.{host}.qc-upload",
               host=get_submit_host(config),
               allow_missing=True)
    params:
        submit_host = get_submit_host(config),
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    log: "{library_id}/post_star_qc.log"
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/post_encode_star_qc_metrics"


rule post_star_solo_qc_metrics:
    input:
        gene_summary = expand("{library_id}/Solo.out/{gene_model}/Summary.csv",
                              gene_model=get_gene_model(config),
                              allow_missing=True),
        umi_plot = "{library_id}/umi_per_cell.png",
        posted = expand("{library_id}/posted.{host}.csv",
                        host=get_submit_host(config),
                        allow_missing=True),
    output:
        expand("{library_id}/Solo.out/{gene_model}/Summary.csv.{host}.qc-upload",
               gene_model=get_gene_model(config),
               host=get_submit_host(config),
               allow_missing=True)
    params:
        submit_host = get_submit_host(config),
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    log: "{library_id}/post_star_solo_qc.log"
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/post_encode_star_solo_qc_metrics"


rule post_single_cell_mex_qc_plots:
    input:
        archive = "{library_id}/{gene_model}_{multiread}_{matrix}.tar.gz",
        pct_mt_plot = "{library_id}/pct_count_mt.{gene_model}_{multiread}_{matrix}.png",
        genes_by_count_plot = "{library_id}/n_genes_by_counts.{gene_model}_{multiread}_{matrix}.png",
        counts_violin_plot = "{library_id}/qc_metric_violin.{gene_model}_{multiread}_{matrix}.png",
        posted = expand("{library_id}/posted.{host}.csv",
                        host=get_submit_host(config),
                        allow_missing=True),
    output:
        expand("{library_id}/pct_count_mt.{gene_model}_{multiread}_{matrix}.png.{host}.qc-upload",
               host=get_submit_host(config),
               allow_missing=True)
    params:
        submit_host = get_submit_host(config),
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    log: "{library_id}/post_count_matrix_qc_{gene_model}_{multiread}_{matrix}.log"
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/post_encode_star_solo_gene_count_plots"
        

rule generate_merged_mex_archives:
    input:
        expand("{library_id}/{gene_model}_{{multiread}}_raw.tar.gz",
               library_id=config["library"].keys(),
               gene_model=get_gene_model(config)),
    output:
        matrix="{gene_model}_{{multiread}}_raw_merged.tar.gz".format(gene_model=get_gene_model(config)),
        correlations="{gene_model}_{{multiread}}_spearman.tsv".format(gene_model=get_gene_model(config)),
        correlation_plot="{gene_model}_{{multiread}}_spearman.png".format(gene_model=get_gene_model(config)),
    params:
        gene_model=get_gene_model(config)
    threads: 1
    resources:
        mem_mb = 8196
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/create_merged_mex_archive"

rule prepare_encode_merged_splitseq_metadata:
    input:
        gene_unique_raw = expand("{gene_model}_Unique_raw_merged.tar.gz",
                                 gene_model=get_gene_model(config)),
        gene_unique_raw_md5 = expand("{gene_model}_Unique_raw_merged.tar.gz.md5",
                                     gene_model=get_gene_model(config)),
        gene_multi_raw = expand("{gene_model}_EM_raw_merged.tar.gz",
                                gene_model=get_gene_model(config)),
        gene_multi_raw_md5 = expand("{gene_model}_EM_raw_merged.tar.gz.md5",
                                    gene_model=get_gene_model(config)),
    output:
        expand("metadata.{host}.csv",
               host=get_submit_host(config))
    threads: 1
    resources:
        mem_mb = 1000
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/prepare_encode_merged_splitseq_metadata"


        
rule submit_star_solo_merged_data:
    input:
        gene_unique_raw = rules.prepare_encode_merged_splitseq_metadata.input.gene_unique_raw,
        gene_multi_raw = rules.prepare_encode_merged_splitseq_metadata.input.gene_multi_raw,
        metadata = rules.prepare_encode_merged_splitseq_metadata.output[0],
    output:
        gene_unique_raw = "{gene_model}_Unique_raw_merged.tar.gz.{host}.upload".format(
            gene_model=get_gene_model(config),
            host=get_submit_host(config)),
        gene_multi_raw = "{gene_model}_EM_raw_merged.tar.gz.{host}.upload".format(
            gene_model=get_gene_model(config),
            host=get_submit_host(config)),
        posted = expand("posted.{host}.csv", host=get_submit_host(config), allow_missing=True)
    params:
        submit_host = get_submit_host(config),
    log: "submit_star_solo_merged_data.log"
    threads: 1
    resources:
        mem_mb = 1000
    wrapper:
        "https://raw.githubusercontent.com/detrout/woldrnaseq-wrappers/0.0.2/post_encode_star_solo_files"
