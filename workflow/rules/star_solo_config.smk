from encoded_client.encoded import ENCODED, HTTPError
from encoded_client.metadata import compute_dcc_file_accession_from_url
from pathlib import Path

DEFAULT_10X_CB_LENGTH = 16
DEFAULT_MEM_MB = 1000
SOLO_ROOT = Path("Solo.out")
UMI_PER_CELL_PLOT_NAME = "umi_per_cell.png"


def get_gene_model():
    return "GeneFull_Ex50pAS" if config['include_intron'] else "Gene"


def get_submit_host():
    default_host = "test.encodedcc.org" #"www.encodeproject.org"
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

    server = ENCODED(get_submit_host())
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
        server = ENCODED(get_submit_host())
        metadata = get_dcc_accession(server, accession)
        if metadata is not None:
            config.setdefault(
                "inclusion_list_url",
                server.prepare_url(metadata["href"]))
        else:
            print("Please set inclusion_list_url")


configfile: "config.yaml"


# username expand any filenames still in the config file
for key in ['genome_dir']:
    if key in config and config[key].startswith('~'):
        config[key] = str(Path(config[key]).expanduser())

config.setdefault(
    "star_container",
    # unreleased
    "https://woldlab.caltech.edu/~diane/containers/star-bash-dev_EoI-head.sif"
    # to old
    #"https://woldlab.caltech.edu/~diane/containers/star-dev_EoI_2.7.9a-2021-09-10.sif"
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
