import sys
from pathlib import Path

samstats = Path("~/proj/GeorgiScripts/SAMstats.py").expanduser()


def get_samstats_paired_args(libraries, analysis_dir):
    library = libraries.set_index("analysis_dir").loc[analysis_dir]
    if "read_2" in libraries.columns:
        return "-paired"
    else:
        return ""


rule samstats:
    input:
        genome_bam = "{analysis_dir}/{analysis_name}-{genome_name}_genome.bam",
        chrom_info = "{analysis_dir}/{analysis_name}-{genome_name}_genome.chrominfo",
        samtools = "/usr/bin/samtools"
    params:
        paired_args = lambda wildcards: get_samstats_paired_args(libraries, wildcards.analysis_dir)
    output:
        "{analysis_dir}/{analysis_name}-{genome_name}_genome.samstats",
    wildcard_constraints:
        analysis_name = r"[^/]+"
    resources:
        mem_mb = 4096
    threads: 1
    shell:
        "{sys.executable} {samstats} {input.genome_bam} {output} -bam {input.chrom_info} {input.samtools} {params.paired_args}"
