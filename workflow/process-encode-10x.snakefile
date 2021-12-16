#
# Copyright 2021 California Institute of Technology
#
# Licensed under a BSD like license that includes a
# non-endorsement clause.
#

from pathlib import Path

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
## This is now the default value. If you want you can override, but you don't need to.
#genome_dir: "genome"
## The code will try to figure out the genome_accession if you only provide a genome_index_url
#genome_accession: "TSTFF984371"
## The code will attempt to figure out genome_index_url from the ENCODE DCC genome_accession
#genome_index_url: "https://woldlab.caltech.edu/~diane/genome/GRCh38-V29-male-2.7.8a.tar.gz"
## If the genome_accession is set the code will pull these values from the ENCODE portal
#assembly: "GRCh38"
#genome_annotation: "V29"
## Given the accession it will determine the inclusion_list_url
##inclusion_accession: "737K-arc-v1(GEX)"
## if there is no inclusion_accession, it will guess the accession from the inclusion_list_url
#inclusion_list_url: "https://woldlab.caltech.edu/~diane/genome/737K-arc-v1.txt.gz"
## The pipeline now defaults to the container used for my analysis (The version of STAR after 2.7.9a)
#star_container: "https://woldlab.caltech.edu/~diane/containers/star-dev_EoI_2.7.9a-2021-09-10.sif"
# This is a container with python and scanpy in it. It also defaults to my container.
#scanpy_container: "https://woldlab.caltech.edu/~diane/containers/bullseye-scanpy-1.8.2.sif"

#
## These are guesses based on my current experiments and might need to change
#mem_mb: 65536
#disk_mb: 51200

include:
    "rules/star_solo_10x_config.smk"



def list_default_targets():
    default = [
        "Log.final.out",
        "Aligned.sortedByCoord.out.bam",
        "{}_Unique_filtered.tar.gz".format(get_gene_model()),
        "{}_EM_filtered.tar.gz".format(get_gene_model()),
        "{}_Unique_raw.tar.gz".format(get_gene_model()),
        "{}_EM_raw.tar.gz".format(get_gene_model()),
        "SJ_Unique_raw.tar.gz",
        "umi_per_cell.png",
        "pct_count_mt.{}_EM_filtered.png".format(get_gene_model()),
        "n_genes_by_counts.{}_EM_filtered.png".format(get_gene_model()),
        "metadata.{}.csv".format(get_submit_host()),
    ]
    if automatic_submission:
        default.extend([
            "posted.{}.csv".format(get_submit_host()),
            "Log.final.out.{}.qc-upload".format(get_submit_host()),
            SOLO_ROOT / get_gene_model() / "Summary.csv.{}.qc-upload".format(
                get_submit_host()
            ),
            "pct_count_mt.{}_EM_filtered.png.{}.qc-upload".format(
                get_gene_model(), get_submit_host()
            ),
        ])
    return default


rule ALL:
    input:
        list_default_targets()

include:
    "rules/provide_star_index.smk"

include:
    "rules/fastqs.smk"

include:
    "rules/star_solo_10x.smk"

include:
    "rules/archive_solo_counts.smk"

include:
    "rules/prepare-qc.smk"

include:
    "rules/encode-submit.smk"
