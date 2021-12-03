import sys
import pandas
from encoded_client.hashfile import make_md5sum
from encoded_client.encoded import ENCODED
from encoded_client.submission import process_files
from encoded_client.metadata import generate_star_solo_processed_metadata

# What are we uploading?
# bam file
# the 5 mex archives
# QC Summary.csv, Features.stats, UMIperCellSorted.txt & plot

rule prepare_md5:
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


rule prepare_star_solo_10x_submission_metadata:
    input:
        bam = "Aligned.sortedByCoord.out.bam",
        bam_md5 = "Aligned.sortedByCoord.out.bam.md5",
        gene_unique_filtered = "{}_Unique_filtered.tar.gz".format(get_gene_model()),
        gene_unique_filtered_md5 = "{}_Unique_filtered.tar.gz.md5".format(get_gene_model()),
        gene_multi_filtered = "{}_EM_filtered.tar.gz".format(get_gene_model()),
        gene_multi_filtered_md5 = "{}_EM_filtered.tar.gz.md5".format(get_gene_model()),
        gene_unique_raw = "{}_Unique_raw.tar.gz".format(get_gene_model()),
        gene_unique_raw_md5 = "{}_Unique_raw.tar.gz.md5".format(get_gene_model()),
        gene_multi_raw = "{}_EM_raw.tar.gz".format(get_gene_model()),
        gene_multi_raw_md5 = "{}_EM_raw.tar.gz.md5".format(get_gene_model()),
        sj_unique_raw = "SJ_Unique_raw.tar.gz",
        sj_unique_raw_md5 = "SJ_Unique_raw.tar.gz",
    output:
        "metadata.csv"
    threads: 1
    resources:
        mem_mb = 100
    run:
        metadata = generate_star_solo_metadata(config, {
            "alignments": input.bam,
            "sparse gene count matrix of unique reads": input.gene_unique_filtered,
            "sparse gene count matrix of all reads": input.gene_multi_filtered,
            "unfiltered sparse gene count matrix of unique reads": input.gene_unique_raw,
            "unfiltered sparse gene count matrix of all reads": input.gene_multi_raw,
            "unfiltered sparse splice junction count matrix of unique reads": input.sj_unique_raw,
        })
        metadata = pandas.DataFrame(metadata)
        metadata.to_csv(output[0], index=False)


rule submit_processed_data:
    input:
        bam = "Aligned.sortedByCoord.out.bam",
        gene_unique_filtered = "{}_Unique_filtered.tar.gz".format(get_gene_model()),
        gene_multi_filtered = "{}_EM_filtered.tar.gz".format(get_gene_model()),
        gene_unique_raw = "{}_Unique_raw.tar.gz".format(get_gene_model()),
        gene_multi_raw = "{}_EM_raw.tar.gz".format(get_gene_model()),
        sj_unique_raw = "SJ_Unique_raw.tar.gz",
        metadata = "metadata.csv"
    output:
        bam = "Aligned.sortedByCoord.out.bam.upload",
        gene_unique_filtered = "{}_Unique_filtered.tar.gz.upload".format(get_gene_model()),
        gene_multi_filtered = "{}_EM_filtered.tar.gz.upload".format(get_gene_model()),
        gene_unique_raw = "{}_Unique_raw.tar.gz.upload".format(get_gene_model()),
        gene_multi_raw = "{}_EM_raw.tar.gz.upload".format(get_gene_model()),
        sj_unique_raw = "SJ_Unique_raw.tar.gz.upload",
        metadata_posted = "posted.csv",
    threads: 1
    resources:
        mem_mb = 100
    run:
        host = config.get("encode_portal_host", "www.encodeproject.org")
        server = ENCODED(host)
        metadata = pandas.read_csv(input.metadata, index_col=None)
        uploaded = process_files(server, metadata, dry_run=False)
        print("Processed {} files".format(len(uploaded)))
        print(metadata)
        metadata.to_csv(output.metadata_posted, index=False)
