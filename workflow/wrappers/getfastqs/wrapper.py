from snakemake import shell
from pathlib import Path
import os
import re
import shutil
import requests

DEFAULT_MEM_MB = 1000


def get_encode_fastq(config, target):
    from encoded_client.encoded import ENCODED
    filename_re = re.compile("(?P<accession>(ENC|TST)FF[0-9A-Z]+)_(?P<read>[RI][1-3])\\.fastq\\.(?P<archive>gz|bz2|xz)")
    host = config.get("encode_portal_host", "www.encodeproject.org")
    server = ENCODED(host)
    match = filename_re.match(target.name)
    if match is None:
        raise ValueError(
            "Unrecognized filename expecting {{accession}}_{{read}}.fastq.(gz|bz2|xz). Got {}".format(target))
    accession = match.group("accession")
    archive = match.group("archive")
    path = "/files/{accession}/@@download/{accession}.fastq.{archive}".format(
        accession=accession, archive=archive)
    url = server.prepare_url(path)
    with requests.get(url, auth=server.auth, stream=True) as instream:
        instream.raise_for_status()
        with open(target, "wb") as outstream:
            shutil.copyfileobj(instream.raw, outstream)


def get_sra_fastq(config, target):
    targets = {
        "1": Path(target.replace("_R2.fastq.gz", "_R1.fastq.gz")),
        "2": Path(target.replace("_R1.fastq.gz", "_R2.fastq.gz")),
    }
    target = Path(target)

    if not target.exists():
        shell("fastq-dump --split-files -v --gzip {wildcards.accession}")

        for read in ["1", "2"]:
            source = "{}_{}.fastq.gz".format(snakemake.wildcards.accession, read)
            os.rename(source, targets[read])


for target in snakemake.output:
    target = Path(target)
    if target.name.startswith("ENCFF"):
        get_encode_fastq(snakemake.config, target)
    elif target.name.startswith("TSTFF"):
        get_encode_fastq(snakemake.config, target)
    elif target.name.startswith("SRX"):
        get_sra_fastq(snakemake.config, target)

