import sys

from woldrnaseq.snakeutils import (
    get_genome_dir,
    get_rsem_gene_results,
    all_experiment_correlations,
)


rule compute_correlations:
    input:
        lambda wildcards: get_rsem_gene_results(
            config, libraries, experiments),
        config["libraries"],
        config["experiments"],
    output:
        all_experiment_correlations(experiments)
    params:
        formatted_libraries = ["-l {}".format(x) for x in config["libraries"]],
        formatted_experiments = ["-e {}".format(x) for x in config["experiments"]],
    resources:
        # should scale with size of max (genome * number of experiments)
        mem_mb = 4000
    threads: 1
    shell:
        "{sys.executable} -m woldrnaseq.madqc -q TPM {params.formatted_libraries} {params.formatted_experiments}"
