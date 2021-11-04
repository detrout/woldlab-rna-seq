from pathlib import Path
import requests
import tarfile

DEFAULT_MEM_MB = 1000


#
# Rules run under singularity have their singularity home directory
# set to the snakemake working directory so can only access files
# in subdirectories beneath them. So my old habit of putting an
# already extracted index somewhere else doesn't work under singularity.
#
# This rule is intended to either make sure there is a genome directory
# for singularity-args to have an already unpacked genome index mounted
# into or to download a genome index tarfile into.
#


rule genome:
    output:
        genome_dir = temp(directory(config.get("genome_dir", "genome")))
    resources:
        mem_mb = DEFAULT_MEM_MB,
    threads: 1
    message: "Preparing genome index in {output.genome_dir}"
    run:
        print(singularity_args)
        genome_dir = Path(output.genome_dir)
        if not Path(Path(genome_dir).parts[0]).exists():
            os.mkdir(genome_dir.parts[0])
        if "genome_index_url" in config:
            with requests.get(config["genome_index_url"], stream=True) as instream:
                tar = tarfile.open(fileobj=instream.raw, mode="r:*")
                # this is a vulnerability, only use trusted tar files
                # See https://docs.python.org/3/library/tarfile.html#tarfile.TarFile.extractall
                tar.extractall(path=genome_dir.parts[0])
