
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

from encoded_client.encoded import ENCODED, HTTPError
from encoded_client.metadata import compute_dcc_file_accession_from_url
from pathlib import Path

DEFAULT_MEM_MB = 1000
SOLO_ROOT = Path("Solo.out")
UMI_PER_CELL_PLOT_NAME = "umi_per_cell.png"


def get_gene_model(config):
    return "GeneFull_Ex50pAS" if config['include_intron'] else "Gene"


def get_submit_host(config):
    default_host = "test.encodedcc.org"   #"www.encodeproject.org"
    return config.get("encode_portal_host", default_host)


def attributes_included(config, attributes):
    not_found = []
    for name in attributes:
        if name not in config:
            not_found.append(name)

    return not_found


def get_dcc_accession(server, accession):
    metadata = None
    try:
        metadata = server.get_json(accession)
    except HTTPError as e:
        if e.response.status_code == 404:
            print("Accession {} not found, discovery will fail".format(accession))
        else:
            print("Other error {} for {}".format(e.response.status_code, accession))
            raise e

    return metadata


def update_genome_annotation_info(config):
    all_genome_attributes = [
        "genome_accession", "genome_index_url", "assembly", "genome_annotation"
    ]
    if len(attributes_included(config, all_genome_attributes)) == 0:
        return

    if "genome_accession" in config:
        genome_accession = config["genome_accession"]
    elif "genome_index_url" in config:
        index_url = config["genome_index_url"]
        genome_accession = compute_dcc_file_accession_from_url(index_url)
        config["genome_accession"] = genome_accession
    else:
        raise ValueError("genome_accession or genome_index_url are required parameters")

    server = ENCODED(get_submit_host(config))
    metadata = get_dcc_accession(server, genome_accession)
    if metadata is not None:
        config.setdefault("assembly", metadata["assembly"])
        config.setdefault("genome_annotation", metadata["genome_annotation"])
        config.setdefault("genome_index_url", server.prepare_url(metadata["href"]))
    else:
        print("Unable to retrieve accession please make sure assembly and genome_index_url are set")


def update_exclusion_info(config):
    all_inclusion_attributes = ["inclusion_accession", "inclusion_list_url"]
    if len(attributes_included(config, all_inclusion_attributes)) == 0:
        return

    if "inclusion_list_url" in config:
        config.setdefault(
            "inclusion_accession",
            compute_dcc_file_accession_from_url(config["inclusion_list_url"]))
    elif "inclusion_accession" in config:
        accession = config["inclusion_accession"]
        server = ENCODED(get_submit_host(config))
        metadata = get_dcc_accession(server, accession)
        if metadata is not None:
            config.setdefault(
                "inclusion_list_url",
                server.prepare_url(metadata["href"]))
        else:
            print("Please set inclusion_list_url")

# username expand any filenames still in the config file
for key in ['genome_dir']:
    if key in config and config[key].startswith('~'):
        config[key] = str(Path(config[key]).expanduser())

config.setdefault(
    "star_container",
    "https://woldlab.caltech.edu/~diane/containers/star-bash-2.7.10a.sif"
)
config.setdefault(
    "scanpy_container",
    "https://woldlab.caltech.edu/~diane/containers/bullseye-scanpy-1.8.2.sif"
)
config.setdefault(
    "alias_prefix",
    Path(config["lab"]).name
)
config.setdefault(
    "genome_dir",
    "genome"
)
update_genome_annotation_info(config)
update_exclusion_info(config)


try:
    automatic_submission = bool(config.get("automatic_submission", False))
except ValueError:
    print("Unable to parse '{}'".format(config.get("automatic_submission")))
    automatic_submission = False
