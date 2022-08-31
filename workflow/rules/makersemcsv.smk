import sys
from pathlib import Path

from woldrnaseq.snakeutils import (
    get_genome_dir,
    get_rsem_gene_results,
    all_experiment_quantifications,
)


rule makersemcsv_genes:
    input:
        # This is wrong. This depends on just some of the rsem files
        # but the program generates all of them.
        lambda wildcards: get_rsem_gene_results(
            config, libraries, experiments)
    output:
        all_experiment_quantifications(experiments)
    params:
        formatted_libraries = ["-l {}".format(x) for x in config["libraries"]],
        formatted_experiments = ["-e {}".format(x) for x in config["experiments"]],
        genome_dir = config["genome_dir"]
    shell:
        """
        {sys.executable} -m woldrnaseq.makersemcsv -q TPM -q FPKM \
           {params.formatted_libraries} {params.formatted_experiments} \
           --genome-dir {params.genome_dir} --add-names
        """
