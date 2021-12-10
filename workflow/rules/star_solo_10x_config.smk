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

try:
    automatic_submission = bool(config.get("automatic_submission", False))
except ValueError:
    print("Unable to parse '{}'".format(config.get("automatic_submission")))
    automatic_submission = False
