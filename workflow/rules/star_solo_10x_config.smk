from encoded_client.encoded import ENCODED
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


def update_genome_annotation_info():
    if "genome_accession" in config:
        genome_accession = config["genome_acession"]
    elif "genome_index_url" in config:
        index_url = config["genome_index_url"]
        genome_accession = compute_dcc_file_accession_from_url(index_url)
        config["genome_accession"] = genome_accession
    else:
        raise ValueError("genome_accession or genome_index_url are required parameters")


    if not('assembly' in config or 'genome_annotation' in config):
        try:
            server = ENCODED(get_submit_host())
            index_metadata = server.get_json(genome_accession)
        except HTTPError as e:
            if e.response.status_code == 404:
                print("{} not found, please set assembly and genone_annotation parameters".format(accession))
            else:
                print("Other error {} for {}".format(e.response.status_code, accession))

        if "assembly" not in config:
            config["assembly"] = index_metadata["assembly"]
        if "genome_annotation" not in config:
            config["genome_annotation"] = index_metadata["genome_annotation"]

configfile: "config.yaml"


# username expand any filenames still in the config file
for key in ['genome_dir']:
    if key in config and config[key].startswith('~'):
        config[key] = str(Path(config[key]).expanduser())

# Set defaults for this version of the ENCODE scRNA-seq
# pipeline
config.setdefault(
    "alignment_step_run",
    "barbara-wold:starsolo-alignment-step-run"
)
config.setdefault(
    "quantification_step_run",
    "barbara-wold:starsolo-quantification-step-run",
)
config.setdefault(
    "star_container",
    "https://woldlab.caltech.edu/~diane/containers/star-dev_EoI_2.7.9a-2021-09-10.sif"
)
config.setdefault(
    "scanpy_container",
    "https://woldlab.caltech.edu/~diane/containers/bullseye-scanpy-1.8.2.sif"
)
config.setdefault(
    "alias_prefix",
    Path(config["lab"]).name
)
update_genome_annotation_info()

try:
    automatic_submission = bool(config.get("automatic_submission", False))
except ValueError:
    print("Unable to parse '{}'".format(config.get("automatic_submission")))
    automatic_submission = False
