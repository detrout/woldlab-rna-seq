
from woldrnaseq.snakeutils import get_genome_dir

def rsem_paired_argument(library_id):
    library = libraries.loc[library_id]
    if "read_2" in library and len(library["read_2"]) > 0:
        return "--paired-end"
    else:
        return ""


def rsem_strand_probability(library_id):
    library = libraries.loc[library_id]
    if "stranded" not in libraries.columns:
        return '--forward-prob 0.5'
    else:
        stranded = library["stranded"]
        if stranded == 'forward':
            return '--forward-prob 1'
        elif stranded == 'reverse':
            return '--forward-prob 0'
        else:
            return '--forward-prob 0.5'


rule rsem_quantification:
    input:
        isoforms = "{analysis_dir}/{analysis_name}-{genome_name}_anno.bam",
        genome_index = lambda wildcards: get_genome_dir(
            config, libraries, wildcards.analysis_dir),
    output:
        "{analysis_dir}/{analysis_name}-{genome_name}_anno_rsem.genes.results",
        "{analysis_dir}/{analysis_name}-{genome_name}_anno_rsem.isoforms.results",
    wildcard_constraints:
        analysis_name = r"[^/]+",

    resources:
        mem_mb = 30000
    threads: 16
    params:
        prefix="{analysis_dir}/{analysis_name}-{genome_name}_anno_rsem",
        paired_end=lambda wildcards: rsem_paired_argument(wildcards.analysis_name),
        strand_probability=lambda wildcards: rsem_strand_probability(wildcards.analysis_name),
        curdir = lambda wildcards, input: Path(input.isoforms).parent,
        seed = 12345,
        tempdir = lambda wildcards, input: temp(directory(
            Path(input.isoforms).parent / "{}_{}.tmp".format(
                wildcards.analysis_name, wildcards.genome_name))),
    singularity:
        config["rsem_container"],
    log:
        "{analysis_dir}/{analysis_name}-{genome_name}_quantification.out"
    shell:
        """rsem-calculate-expression --bam --estimate-rspd --calc-ci -p {threads} \
           --seed {params.seed} \
           --no-bam-output --ci-memory {resources.mem_mb} \
           {params.paired_end} {params.strand_probability} \
           --temporary-folder  {params.tempdir} \
           {input.isoforms} \
           {input.genome_index}/rsem \
           {params.prefix} \
        """
