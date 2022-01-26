import sys


rule create_md5s:
    input:
        "{filename}"
    output:
        "{filename}.md5"
    params:
        python = sys.executable,
    threads: 1
    resources:
        mem_mb = 100
    shell:
        "{params.python} -m encoded_client.hashfile {input}"


