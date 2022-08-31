import sys

rule compute_read_coverage:
    input:
        genome_bam = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam",
        cache = lambda wildcards: get_genome_cache(
            config, libraries, wildcards.analysis_dir),
    output:
        coverage = "{analysis_dir}/{analysis_name}-{genome_name}.coverage",
    resources:
        mem_mb = 4000
    threads: 1
    shell:
        """{sys.executable} -m woldrnaseq.gene_coverage_wig_gtf \
              --gtf {input.cache} \
              --output {output.coverage} \
              --gene-normalization max \
              --expression-threshold 1.0 \
              {input.genome_bam}
        """
