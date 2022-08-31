import sys
from woldrnaseq.snakeutils import get_genome_cache


rule gene_coverage:
    input:
        genome_bam = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam",
        genome_bai = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam.bai",
        cache = lambda wildcards: get_genome_cache(
            config, libraries, wildcards.analysis_dir),
    output:
        "{analysis_dir}/{analysis_name}-{genome_name}.coverage",
    threads: 1
    resources:
        mem_mb = 8000
    log:
        "{analysis_dir}/{analysis_name}-{genome_name}.gene_coverage.out"
    shell:
        """{sys.executable} -m woldrnaseq.gene_coverage_wig_gtf \
            --gtf {input.cache} \
            {input.genome_bam} \
            --output {output} \
            --gene-normalization max \
            --expression-threshold 1.0 2>&1 > {log}
        """
