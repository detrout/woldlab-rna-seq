"""Functions that should be shared between multiple snakemake rules
"""
from pathlib import Path

from .models import (
    genome_name_from_library,
    sanitize_name
)


def get_library_by_analysis_dir(libraries, analysis_dir):
    """Return library row from libraries table given a target analysis directory

    Parameters:
      - libraries (pandas.DataFrame): table describing library metadata
      - analysis_dir (str): target directory name

    Returns:
      - (pandas.Series): row of current library metadata
    """
    library = libraries.set_index("analysis_dir").loc[analysis_dir]
    return library


def get_genome_dir(config, libraries, analysis_dir):
    """Compute the genome_dir from library files

    Parameters:
      - config (dict): dictionary of default settings
      - libraries (pandas.DataFrame): table describing library metadata
      - analysis_dir (str): target directory name

    Returns:
      - (str) directory name containing index files
    """
    library = get_library_by_analysis_dir(libraries, analysis_dir)
    genome_name = genome_name_from_library(library)
    genome_dir = Path(config["genome_dir"]) / genome_name
    return str(genome_dir)


def get_library_stranded(config, libraries, analysis_dir):
    """Compute the genome_dir from library files

    Parameters:
      - config (dict): dictionary of default settings
      - libraries (pandas.DataFrame): table describing library metadata
      - analysis_dir (str): target directory name

    Returns:
      - (str) stranded setting for library
    """
    library = get_library_by_analysis_dir(libraries, analysis_dir)
    return library.get("stranded", "unstranded")


def get_genome_cache(config, libraries, analysis_dir):
    """Compute the genome_cach3e from library files

    Parameters:
      - config (dict): dictionary of default settings
      - libraries (pandas.DataFrame): table describing library metadata
      - analysis_dir (str): target directory name

    Returns:
      - (str) directory name containing index files
    """
    library = get_library_by_analysis_dir(libraries, analysis_dir)
    genome_name = genome_name_from_library(library)
    genome_dir = Path(config["genome_dir"]) / genome_name
    cache_basename = str(genome_name) + ".h5"
    cache_filename = genome_dir / cache_basename
    if not cache_filename.exists():
        raise FileNotFoundError("Unable to find cache file, please create with gff2table")
    return str(cache_filename)


def get_rsem_gene_results(config, libraries, experiments, experiment_name=None):
    """Return the list of rsem result files
    """
    quoted = [sanitize_name(x) for x in experiments.index]
    if experiment_name is not None:
        location = quoted.index(experiment_name)
        experiment = experiments.iloc[location]
        current_libraries = libraries.reindex(experiment.replicates)
    else:
        current_libraries = libraries

    dependencies = []
    for i, library in current_libraries.reset_index().iterrows():
        # currently implicitly defined in new_workflow/Snakefile
        dependencies.extend(library_quantification_target(library))

    return dependencies


def library_star_log_target(library):
    return [
        "{analysis_dir}/Log.final.out".format(**library)
    ]


def library_bam_target(library):
    genome_bam = "{analysis_dir}/{library_id}-{genome_name}_genome.bam".format(
        **library)
    transcript_bam = "{analysis_dir}/{library_id}-{genome_name}_anno.bam".format(
        **library)
    return [
        genome_bam,
        genome_bam + ".bai",
        transcript_bam,
    ]


def library_quantification_target(library):
    return [
        "{analysis_dir}/{library_id}-{genome_name}_anno_rsem.genes.results".format(**library),
        "{analysis_dir}/{library_id}-{genome_name}_anno_rsem.isoforms.results".format(
            **library),
    ]


def library_samstats_target(library):
    return [
        "{analysis_dir}/{library_id}-{genome_name}_genome.samstats".format(
            **library)
    ]


def library_coverage_target(library):
    return [
        "{analysis_dir}/{library_id}-{genome_name}.coverage".format(
            **library)
    ]


def library_distribution_target(library):
    return [
        "{analysis_dir}/{library_id}-{genome_name}.sam_reads_genes".format(
            **library)
    ]


def library_bigwigs_target(library):
    if library.stranded.lower() == "unstranded":
        # unstranded bigwigs
        return [
            "{analysis_dir}/{library_id}-{genome_name}_uniq.bw".format(**library),
            "{analysis_dir}/{library_id}-{genome_name}_uniq.bw".format(**library),
        ]
    else:
        # stranded bigwigs
        return [
            "{analysis_dir}/{library_id}-{genome_name}_minusUniq.bw".format(**library),
            "{analysis_dir}/{library_id}-{genome_name}_plusUniq.bw".format(**library),
            "{analysis_dir}/{library_id}-{genome_name}_minusAll.bw".format(**library),
            "{analysis_dir}/{library_id}-{genome_name}_plusAll.bw".format(**library),
        ]


def all_experiment_quantifications(experiments):
    targets = []
    for experiment_name, experiment in experiments.iterrows():
        attributes = {
            "experiment_name": sanitize_name(experiment_name)
        }
        targets.append("{experiment_name}_gene_FPKM.csv".format(**attributes))
        targets.append("{experiment_name}_gene_TPM.csv".format(**attributes))

    return targets


def all_experiment_correlations(experiments):
    targets = []

    for experiment_name, experiment in experiments.iterrows():
        attributes = {
            "experiment_name": sanitize_name(experiment_name)
        }
        targets.append("{experiment_name}_correlation.h5".format(**attributes))

    return targets
