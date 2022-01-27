import sys
import logging
import pandas
from encoded_client.encoded import ENCODED, DCCValidator, make_attachment, HTTPError
from encoded_client.submission import process_files
from encoded_client.metadata import generate_star_solo_processed_sheet

# What are we uploading?
# bam file
# the 5 mex archives
# QC Summary.csv, Features.stats, UMIperCellSorted.txt & plot


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
        sj_unique_raw_md5 = "SJ_Unique_raw.tar.gz.md5",
    output:
        "metadata.{}.csv".format(get_submit_host())
    threads: 1
    resources:
        mem_mb = 100
    run:
        metadata = generate_star_solo_processed_sheet(config, {
            "alignments": input.bam,
            "sparse gene count matrix of unique reads": input.gene_unique_filtered,
            "sparse gene count matrix of all reads": input.gene_multi_filtered,
            "unfiltered sparse gene count matrix of unique reads": input.gene_unique_raw,
            "unfiltered sparse gene count matrix of all reads": input.gene_multi_raw,
            "unfiltered sparse splice junction count matrix of unique reads": input.sj_unique_raw,
        })
        # need to add quality metrics too.
        # Star metric https://www.encodeproject.org/profiles/star_quality_metric (from log file
        # generic metric https://www.encodeproject.org/profiles/generic_quality_metric
        #   Summary.csv as text/plain? text/tab-separated-value?
        #   UMIperCellSorted as image/png
        #   hopefully simple umap.
        # the first 3 metrics should all be attached to the bam file
        # the last UMAP would be attached to which ever count matrix was processed
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
        metadata = rules.prepare_star_solo_10x_submission_metadata.output[0],
    output:
        bam = "Aligned.sortedByCoord.out.bam.{}.upload".format(get_submit_host()),
        gene_unique_filtered = "{}_Unique_filtered.tar.gz.{}.upload".format(get_gene_model(), get_submit_host()),
        gene_multi_filtered = "{}_EM_filtered.tar.gz.{}.upload".format(get_gene_model(), get_submit_host()),
        gene_unique_raw = "{}_Unique_raw.tar.gz.{}.upload".format(get_gene_model(), get_submit_host()),
        gene_multi_raw = "{}_EM_raw.tar.gz.{}.upload".format(get_gene_model(), get_submit_host()),
        sj_unique_raw = "SJ_Unique_raw.tar.gz.{}.upload".format(get_submit_host()),
        metadata_posted = "posted.{}.csv".format(get_submit_host()),
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    log: "submit_processed_data.log"
    run:
        logger = logging.getLogger("submit_processed_data")
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.FileHandler(log[0]))
        logger.addHandler(logging.StreamHandler(sys.stderr))
        host = get_submit_host()
        server = ENCODED(host)
        metadata = pandas.read_csv(input.metadata, index_col=None)
        # catch previous submission
        for i, row in metadata.iterrows():
            if pandas.isnull(row["accession"]):
                try:
                    file = server.get_json("md5:{}".format(row["md5sum"]))
                    metadata[i, "accession"] = row["accession"]
                    metadata[i, "uuid"] = row["uuid"]
                    upload_file = Path("{}.{}.upload".format(row["submitted_file_name"], host))
                    if not upload_file.exists():
                        upload_file.touch()
                except HTTPError as e:
                    if e.response.status_code != 404:
                        logger.warning("Unexpected status code {}".format(e))
        uploaded = process_files(server, metadata, dry_run=False)
        logger.info("Processed {} files".format(len(uploaded)))
        logger.info(metadata)
        metadata.to_csv(output.metadata_posted, index=False)


# Should the metadata be in a single file?
# That makes it harder to see what's been processed and pull submitted
# accessions between rules.
