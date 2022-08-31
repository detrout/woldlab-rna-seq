
rule index_genome_bam:
    input:
        "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam"
    output:
        "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam.bai"
    resources:
        mem_mb = 16000,
    threads: 10
    singularity:
        config["rsem_container"]
    shell:
        "samtools index -b -@ {threads}  {input} {output}"

