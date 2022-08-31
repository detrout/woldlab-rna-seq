import sys
from woldrnaseq.snakeutils import (
    get_genome_cache,
    get_library_stranded
)


rule compute_read_distribution:
    input:
        genome = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam",
        genome_index = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam.bai",
        cache = lambda wildcards: get_genome_cache(
            config, libraries, wildcards.analysis_dir),
    output:
        "{analysis_dir}/{analysis_name}-{genome_name}.sam_reads_genes"
    params:
        strand = lambda wildcards: get_library_stranded(
            config, libraries, wildcards.analysis_dir)
    threads: 8
    resources:
        mem_mb = 16000
    shell:
        """{sys.executable} -m woldrnaseq.compute_read_distribution \
           --strand {params.strand} \
           --gtf-cache {input.cache} \
           -o {output} {input.genome}"""
