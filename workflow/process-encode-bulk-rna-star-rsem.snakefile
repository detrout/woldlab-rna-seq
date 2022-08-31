import pandas
from pathlib import Path
from woldrnaseq.common import (
    find_fastqs,
    sanitize_name,
)
from woldrnaseq.models import (
    load_library_tables,
    load_experiments,
    genome_name_from_library,
)
from woldrnaseq.snakeutils import (
    library_star_log_target,
    library_bam_target,
    library_quantification_target,
    library_samstats_target,
    library_coverage_target,
    library_distribution_target,
    library_bigwigs_target,
    all_experiment_quantifications,
    all_experiment_correlations,
)

from configparser import SafeConfigParser

htsw = SafeConfigParser()
htsw.read([
    Path("~/.htsworkflow.ini").expanduser(),
    Path("/etc/htsworkflow.ini")])

if "analysis" in htsw:
    for term in ["genome_dir", "star_dir", "rsem_dir", "georgi_dir", "ucsc_tools_dir"]:
        if term in htsw["analysis"]:
            config.setdefault(term, htsw.get("analysis", term))


config["star_container"] = str(Path("~/public_html/containers/star-2.5.1b.sif").expanduser())
config["rsem_container"] = str(Path("~/public_html/containers/rsem-1.2.31-v0.0.3.sif").expanduser())
config["ucsc_tools_container"] = str(Path("~/public_html/containers/ucsc-tools-433.sif").expanduser())
config.setdefault("experiments", ["experiments.tsv"])
config.setdefault("libraries", ["libraries.tsv"])

for term in ["experiments", "libraries"]:
    if isinstance(config[term], str):
        config[term] = config[term].split(',')

experiments = load_experiments(config["experiments"])
libraries = load_library_tables(config["libraries"])

#libraries["read_1"] = find_fastqs(libaries
libraries["read_1"] = pandas.Series(dict(find_fastqs(libraries, "read_1")))
if "read_2" in libraries.columns:
    libraries["read_2"] = pandas.Series(dict(find_fastqs(libraries, "read_2")))



def default_targets():
    targets = []

    # per library targets
    for i, library in libraries.reset_index().iterrows():
        targets.extend(library_star_log_target(library))
        targets.extend(library_bam_target(library))
        targets.extend(library_quantification_target(library))
        targets.extend(library_samstats_target(library))
        targets.extend(library_coverage_target(library))
        targets.extend(library_distribution_target(library))
        targets.extend(library_bigwigs_target(library))

    # experiment targets
    targets.extend(all_experiment_quantifications(experiments))
    targets.extend(all_experiment_correlations(experiments))
    targets.append("report.html")

    return targets

def genome_name_re(libraries):
    def escape_re(name):
        name = name.replace("\\", "\\\\")
        name = name.replace("+", "\\+")
        name = name.replace("*", "\\*")
        name = name.replace("[", "\\[")
        name = name.replace("]", "\\]")
        return name
    print("escape", escape_re("a+-[foo]"))
    print("Unique", libraries["genome_name"].unique())
    names = [escape_re(name) for name in libraries["genome_name"].unique()]
    names = "|".join(names)
    return names
    
wildcard_constraints:
    genome_name = genome_name_re(libraries)

include: "../../woldlab-rna-seq/workflow/rules/bulk-rna-seq-align.smk"
include: "../../woldlab-rna-seq/workflow/rules/index-bam.smk"
include: "../../woldlab-rna-seq/workflow/rules/bulk-rna-seq-quantification.smk"
include: "../../woldlab-rna-seq/workflow/rules/gene-coverage.smk"
include: "../../woldlab-rna-seq/workflow/rules/qc-samstats.smk"
include: "../../woldlab-rna-seq/workflow/rules/read-distribution.smk"
include: "../../woldlab-rna-seq/workflow/rules/bam2bigwig.smk"
include: "../../woldlab-rna-seq/workflow/rules/makersemcsv.smk"
include: "../../woldlab-rna-seq/workflow/rules/compute_correlations.smk"
include: "../../woldlab-rna-seq/workflow/rules/bulk-report.smk"


rule ALL:
    input:
        default_targets()
