import yaml
import gzip
from pathlib import Path
import requests
import shutil
import subprocess
from urllib.parse import urlparse
from woldrnaseq import gff2table
from woldrnaseq.common import (
    get_star_version,
    get_rsem_version,
)
from xopen import xopen
DEFAULT_MEM_MB = 1000

configfile: "config.yaml"

# expand ~
for p in ["output_dir", "star_dir", "rsem_dir"]:
    if p in config:
        config[p] = Path(config[p]).expanduser()

config.setdefault("output_dir", Path.cwd())
config.setdefault("name", Path(config["output_dir"]).name)


def get_star_command(config):
    star = "STAR"
    star_dir = config.get("star_dir")
    if star_dir is not None:
        return star_dir / star
    else:
        return star

def get_rsem_prepare_reference_command(config):
    rsem = "rsem-prepare-reference"
    rsem_dir = config.get("rsem_dir")
    if rsem_dir is not None:
        return rsem_dir / rsem
    else:
        return rsem


def save_bamcomment(filename, config, version):
    with open(filename, "wt") as outstream:
        for key in ["refid", "annid", "spikeid"]:
            outstream.write("@CO\t{key}:{value}\n".format(
                key=key.upper(), value=config[key]))
        outstream.write("@CO\tINDEXERVER:{}\n".format(version))


def get_gff_archive_name(config):
    return config["name"] + ".gtf.gz"


def get_gffcache_name(config):
    return config['name'] + ".h5"


def get_star_archive_name(config):
    return "{}-star.tar.gz".format(config["name"])


def get_rsem_archive_name(config):
    return "{}-rsem.tar.gz".format(config["name"])


def download(source, destination):
    parsed = urlparse(source)
    destination = Path(destination)
    if parsed.scheme in ['https', 'http']:
        print("Downloading {} to {}".format(source, destination))
        with open(destination, "wb") as outstream:
            with requests.get(source, stream=True) as instream:
                instream.raise_for_status()
                shutil.copyfileobj(instream.raw, outstream)
    elif parsed.scheme in ['file', '']:
        filename = Path(parsed.path)
        print("Linking {} to {}".format(destination, filename))
        destination.symlink_to(filename)
    else:
        print("Unrecognized scheme {}".format(parsed.scheme))

rule ALL:
    input:
        "SA",
        "rsem.idx.fa",
        "star_bamCommentLines.txt",
        "rsem_bamCommentLines.txt",
        get_gff_archive_name(config),
        get_gffcache_name(config),
        get_star_archive_name(config),
        get_rsem_archive_name(config),


rule star_comment:
    input:
        "config.yaml"
    output:
        "star_bamCommentLines.txt"
    resources:
        mem_mb = 100
    threads: 1
    run:
        save_bamcomment(output[0], config, config["star_version"])

rule rsem_comment:
    input:
        "config.yaml"
    output:
        "rsem_bamCommentLines.txt"
    resources:
        mem_mb = 100
    threads: 1
    run:
        save_bamcomment(output[0], config, config["rsem_version"])

rule download:
    output:
        fasta = temp(Path(config['fasta']).name),
        spikeins = [temp(Path(x).name) for x in config["spikeins"]],
        trna = temp(Path(config["trna"]).name),
        annotation = temp(Path(config["gtf"]).name),
    resources:
        mem_mb = 100
    threads: 1
    run:
        for name in ["fasta", "trna", "gtf"]:
            download(config[name], Path(config[name]).name)
        for spike in config["spikeins"]:
            download(spike, Path(spike).name)

rule prepare_reference:
    input:
        fasta = Path(config['fasta']).name,
        spikeins = [Path(x).name for x in config["spikeins"]],
    output:
        temp(config["output_dir"] / "{}.fasta".format(config["name"]))
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    run:
        files = [input.fasta]
        if isinstance(input.spikeins, list):
            files.extend(input.spikeins)
        else:
            files.append(input.spikeins)
        with open(output[0], "wt") as outstream:
            for filename in files:
                with xopen(filename, "rt") as instream:
                    shutil.copyfileobj(instream, outstream)

rule prepare_annotation:
    input:
        trna = Path(config["trna"]).name,
        gtf = Path(config["gtf"]).name,
        spikeins = [Path(x).name for x in config["spikeins"]],
    output:
        temp(config["output_dir"] / "{}.gtf".format(config["name"]))
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    run:
        from woldrnaseq.merge_encode_annotations import main as annotation_main
        args = ["--trna", input.trna, input.gtf]
        if isinstance(input.spikeins, list):
            for spikein in input.spikeins:
                args.append("--spikein")
                args.append(spikein)
        else:
            args.extend(["--spikein", config["spikeins"]])
        args.extend(["--output", output[0]])
        print("merge", args)
        annotation_main(args)


rule star_index:
    input:
        fasta = rules.prepare_reference.output[0],
        gtf = rules.prepare_annotation.output[0],
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
        output_dir = config["output_dir"],
        star_cmd = get_star_command(config),
    singularity:
        config["star_container"]
    shell:
        "{params.star_cmd} --runThreadN {threads} \
            --runMode genomeGenerate \
            --genomeDir {params.output_dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang 100"

rule star_archive:
    input:
        rules.star_comment.output,
        rules.star_index.output,
    output:
        get_star_archive_name(config)
    resources:
        mem_mb = 512
    threads: 1
    shell:
        "tar czvf {output} {input}"

rule rsem_index:
    input:
        fasta = rules.prepare_reference.output[0],
        gtf = rules.prepare_annotation.output[0],
    params:
        output_dir = config['output_dir'],
        rsem_cmd = get_rsem_prepare_reference_command(config)
    singularity:
        config["rsem_container"],
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
        "{params.rsem_cmd} --gtf {input.gtf} {input.fasta} {params.output_dir}/rsem"


rule rsem_archive:
    input:
        rules.rsem_comment.output,
        rules.rsem_index.output,
    output:
        get_rsem_archive_name(config)
    resources:
        mem_mb = 512
    threads: 1
    shell:
        "tar czvf {output} {input}"


rule gff_archive:
    input:
        gtf = rules.prepare_annotation.output[0],
    output:
        get_gff_archive_name(config)
    threads: 1
    resources:
        mem_mb = 256
    run:
        with open(input.gtf, "rb") as instream:
            with gzip.GzipFile(output[0], "wb") as outstream:
                shutil.copyfileobj(instream, outstream)

rule gffcache:
    input:
        gtf = rules.prepare_annotation.output[0],
    output:
        get_gffcache_name(config),
    threads: 1
    resources:
        mem_mb = 2000
    run:
        gff2table.main(["-o", output[0], input.gtf, "--verbose"])
