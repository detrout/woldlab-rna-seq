from woldrnaseq import models
from pathlib import Path
import yaml

config = models.config_from_condor_initialdir()


def generate_cram_targets(config):
    targets = []
    cram_suffixes = ["_anno.cram", "_genome.cram"]

    genome_dir = Path(config["genome_dir"]) / config["genome_triple"]
    genome_config = yaml.open(genome_dir / "config.yaml", "r")
    config["genome_fasta"] = genome_dir / genome_config['fasta']
    config["transcript_fasta"] = genome_dir / "rsem.transcripts.fa"

    targets.append(str(config["genome_fasta"]) + '.fai')
    targets.append(str(config["transcript_fasta"]) + '.fai')

    for suffix in suffixes:
        targets.append("{analysis_name}-{genome_triple}{suffix}".format(
            analysis_name=config["analysis_name"],
            genome_triple=config["genome_triple"],
            suffix=suffix))
        if os.path.exists(targets[-1] + '.bai'):
            pass
    return targets


rule ALL:
    input:
        generate_cram_targets(config)

rule index_fasta:
    input:
        "{fasta_basename}.fa"
    output:
        "{fasta_basename}.fa.fai"
    threads: 1
    log: "{fasta_basename}.faidx.log"
    shell:
        "samtools faidx {input} > {log}"

rule cram_genome:
    input:
        bam = "{bam}_genome.bam"
        fasta = config['genome_fasta']
    output:
        "{bam}_genome.cram"
    log: "{bam}_genome_cram.log"
    threads: 8
    shell:
        "samtools view -T {input.fasta} -@{threads} {input.bam} -o {output}"

rule cram_transcript:
    input:
        bam = "{bam}_anno.bam"
        fasta = config['transcript_fasta']
    output:
        "{bam}_anno.cram"
    log: "{bam}_genome_cram.log"
    threads: 8
    shell:
        "samtools view -T {input.fasta} -@{threads} {input.bam} -o {output}"
