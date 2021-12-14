import yaml
from pathlib import Path
import subprocess
from woldrnaseq import gff2table
from woldrnaseq.common import (
    get_star_version,
    get_rsem_version,
)

with open("config.yaml", "rt") as instream:
    config = yaml.load(instream)

# expand ~
for p in ["output_dir", "star_dir", "rsem_dir"]:
    config[p] = Path(config[p]).expanduser()


def save_bamcomment(filename, config, version):
    with open(filename, "wt") as outstream:
        for key in ["refid", "annid", "spikeid"]:
            outstream.write("@CO\t{key}:{value}\n".format(
                key=key.upper(), value=config[key]))
        outstream.write("@CO\tINDEXERVER:{}\n".format(version))


def get_gffcache_name(config):
    return config['name'] + ".h5"


rule ALL:
    input:
        get_gffcache_name(config),
        "SA",
        "rsem.idx.fa",
        "star_bamCommentLines.txt",
        "rsem_bamCommentLines.txt",


rule star_comment:
    input:
        "config.yaml"
    output:
        "star_bamCommentLines.txt"
    params:
        version = get_star_version(config['star_dir'])
    run:
        save_bamcomment(output[0], config, params.version)

rule rsem_comment:
    input:
        "config.yaml"
    output:
        "rsem_bamCommentLines.txt"
    params:
        version = get_rsem_version(config['rsem_dir'])
    run:
        save_bamcomment(output[0], config, params.version)

rule star_index:
    input:
        fasta = config['fasta'],
        gtf = config['gtf'],
        star_dir = config['star_dir'],
    output:
        "chrLength.txt",
        "chrName.txt",
        "chrNameLength.txt",
        "chrStart.txt",
        "exonGeTrInfo.tab",
        "exonInfo.tab",
        "Genome",
        "genomeParameters.txt",
        "geneInfo.tab",
        "Log.out",
        "SA",
        "SAindex",
        "sjdbInfo.txt",
        "sjdbList.fromGTF.out.tab",
        "sjdbList.out.tab",
        "transcriptInfo.tab",
    threads: 12
    resources:
        mem_mb = 40000
    params:
        output_dir = config["output_dir"]
    shell:
        "{input.star_dir}/STAR --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {params.output_dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100"

rule rsem_index:
    input:
        fasta = config['fasta'],
        gtf = config['gtf'],
        output_dir = config['output_dir'],
        rsem_dir = config['rsem_dir'],
    output:
        "rsem.chrlist",
        "rsem.grp",
        "rsem.idx.fa",
        "rsem.n2g.idx.fa",
        "rsem.seq",
        "rsem.ti",
        "rsem.transcripts.fa",
    threads: 1
    resources:
        mem_mb = 40000
    shell:
        "{input.rsem_dir}/rsem-prepare-reference --gtf {input.gtf} {input.fasta} {params.output_dir}/rsem"


rule gffcache:
    input:
        gtf = config['gtf'],
    output:
        get_gffcache_name(config),
    threads: 1
    resources:
        mem_mb = 2000
    run:
        gff2table.main(["-o", output[0], input.gtf, "--verbose"])
