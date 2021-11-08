from urllib.parse import urljoin
from pathlib import Path

DEFAULT_MEM_MB = 1000


rule get_encode_fastq:
    output:
        "{accession}_R{read}.fastq.gz"
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB,
    wildcard_constraints:
        accession = "ENCFF.*"
    run:
        from encoded_client.encoded import ENCODED
        host = config.get("encode_portal_host", "www.encodeproject.org")
        server = ENCODED(host)
        path = "/files/{accession}/@@download/{accession}.fastq.gz".format(
            accession=wildcards.accession)
        url = server.prepare_url(path)
        with requests.get(url, auth=server.auth, stream=True) as instream:
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


