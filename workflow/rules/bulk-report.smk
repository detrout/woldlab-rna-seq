import sys

from woldrnaseq.snakeutils import (
    library_star_log_target,
    library_quantification_target,
    library_coverage_target,
    library_distribution_target,
    all_experiment_correlations,
)
    

def get_report_dependencies(config, libraries, experiments):
    dependencies = []
    for i, library in libraries.reset_index().iterrows():
        # currently implicitly defined in new_workflow/Snakefile
        dependencies.extend(library_star_log_target(library))
        dependencies.extend(library_quantification_target(library))
        dependencies.extend(library_samstats_target(library))
        dependencies.extend(library_coverage_target(library))
        dependencies.extend(library_distribution_target(library))

    dependencies.extend(all_experiment_correlations(experiments))
    return dependencies


rule generate_report:
    input:
        get_report_dependencies(config, libraries, experiments)
    output:
        "report.html"
    params:
        formatted_libraries = ["-l {}".format(x) for x in config["libraries"]],
        formatted_experiments = ["-e {}".format(x) for x in config["experiments"]],
    shell:
        """{sys.executable} -m woldrnaseq.report \
             -q TPM \
             {params.formatted_libraries} \
             {params.formatted_experiments}"""
