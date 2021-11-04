
from pathlib import Path

DEFAULT_MEM_MB = 1000


def get_netrc_auth(host):
    authdb = netrc.netrc()
    username, _, password = authdb.hosts[host]
    auth = (username, password)
    return auth


def get_auth(config):
    host = config.get("encode_portal_url", "www.encodeproject.org")
    home_netrc = Path("~/.netrc").expanduser()

    if home_netrc.exists():
        auth = get_netrc_auth(host)
    else:
        username = os.environ.get("DCC_API_KEY")
        password = os.environ.get("DCC_SECRET_KEY")
        auth = (username, password)

    if auth == (None, None):
        return None
    else:
        return auth


rule get_encode_fastq:
    output:
        "{accession}_R{read}.fastq.gz"
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB,
    wildcard_constraints:
        accession = "ENCFF.*"
    run:
        url = "https://www.encodeproject.org/files/{accession}/@@download/{accession}.fastq.gz".format(
            accession=wildcards.accession)
        with requests.get(url, auth=auth, stream=True) as instream:
            instream.raise_for_status()
            with open(output[0], "wb") as outstream:
                shutil.copyfileobj(instream.raw, outstream)


rule get_sra_fastq:
    output:
        "{accession}_R1.fastq.gz",
        "{accession}_R2.fastq.gz",
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    wildcard_constraints:
        accession = "SRX.*"
    shell:
        """
        fastq-dump --split-files -v --gzip {wildcards.accession}
        mv {wildcards.accession}_1.fastq.gz {wildcards.accession}_R1.fastq.gz
        mv {wildcards.accession}_2.fastq.gz {wildcards.accession}_R2.fastq.gz
        """


