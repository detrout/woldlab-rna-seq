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
#genome_dir: "genome"
#genome_index_url: "https://woldlab.caltech.edu/~diane/genome/GRCh38-V29-male-2.7.8a.tar.gz"
#allow_list_url: "https://woldlab.caltech.edu/~diane/genome/737K-arc-v1.txt.gz"
#star_container: "https://woldlab.caltech.edu/~diane/containers/star-dev_EoI_2.7.9a-2021-09-10.sif"
#
## These are guesses based on my current experiments and might need to change
#mem_mb: 65536
#disk_mb: 51200


DEFAULT_10X_CB_LENGTH = 16
DEFAULT_MEM_MB = 1000
SOLO_ROOT = Path("Solo.out")

configfile: "config.yaml"

for key in ['genome_dir']:
    if key in config and config[key].startswith('~'):
        config[key] = str(Path(config[key]).expanduser())


###
# Main rules
def get_gene_model():
    return "GeneFull_Ex50pAS" if config['include_intron'] else "Gene"


rule ALL:
    input:
        "Log.final.out",
        "Aligned.sortedByCoord.out.bam",
        "{}_Unique_filtered.tar.gz".format(get_gene_model()),
        "{}_EM_filtered.tar.gz".format(get_gene_model()),
        "{}_Unique_raw.tar.gz".format(get_gene_model()),
        "{}_EM_raw.tar.gz".format(get_gene_model()),
        "SJ_Unique_raw.tar.gz",
        "posted.csv",  #.format(config["experiment_accession"]),

include:
    "rules/provide_star_index.smk"

include:
    "rules/star_solo_10x.smk"

include:
    "rules/archive_solo_counts.smk"

include:
    "rules/fastqs.smk"

include:
    "rules/encode-submit.smk"
