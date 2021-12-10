import sys
import pandas
from encoded_client.hashfile import make_md5sum
from encoded_client.encoded import ENCODED, DCCValidator, make_attachment
from encoded_client.submission import process_files
from encoded_client.metadata import generate_star_solo_processed_metadata

from woldrnaseq.models import (
    load_star_final_log,
    load_star_solo_quality_metric,
)
# What are we uploading?
# bam file
# the 5 mex archives
# QC Summary.csv, Features.stats, UMIperCellSorted.txt & plot


filename_to_output_type = {
    "Aligned.sortedByCoord.out.bam": "alignments",
    "GeneFull_Ex50pAS_Unique_filtered.tar.gz": "sparse gene count matrix of unique reads",
    "GeneFull_Ex50pAS_EM_filtered.tar.gz": "sparse gene count matrix of all reads",
    "GeneFull_Ex50pAS_Unique_raw.tar.gz": "unfiltered sparse gene count matrix of unique reads",
    "GeneFull_Ex50pAS_EM_raw.tar.gz": "unfiltered sparse gene count matrix of all reads",
    "SJ_Unique_raw.tar.gz": "unfiltered sparse splice junction count matrix of unique reads",
}

def prepare_star_qc_metric(config, metric_of, filename):
    star_log = load_star_final_log(filename)
    star_quality_metric = {
        "Mapping speed, Million of reads per hour": star_log.loc[("", "Mapping speed, Million of reads per hour")],
        "Number of input reads": star_log.loc[("", "Number of input reads")],
        "Average input read length": star_log.loc[("", "Average input read length")],

        "Average mapped length": star_log.loc[("UNIQUE READS", "Average mapped length")],
        "Deletion average length": star_log.loc[("UNIQUE READS", "Deletion average length")],
        "Deletion rate per base": "{:0.2f}".format(star_log.loc[("UNIQUE READS", "Deletion rate per base")]),
        "Insertion average length": star_log.loc[("UNIQUE READS", "Insertion average length")],
        "Insertion rate per base": "{:0.2f}".format(star_log.loc[("UNIQUE READS", "Insertion rate per base")]),
        "Mismatch rate per base, %": "{:0.2f}%".format(star_log.loc[("UNIQUE READS", "Mismatch rate per base, %")]),
        "Number of splices: AT/AC": star_log.loc[("UNIQUE READS", "Number of splices: AT/AC")],
        "Number of splices: Annotated (sjdb)": star_log.loc[("UNIQUE READS", "Number of splices: Annotated (sjdb)")],
        "Number of splices: GC/AG": star_log.loc[("UNIQUE READS", "Number of splices: GC/AG")],
        "Number of splices: GT/AG": star_log.loc[("UNIQUE READS", "Number of splices: GT/AG")],
        "Number of splices: Non-canonical": star_log.loc[("UNIQUE READS", "Number of splices: Non-canonical")],
        "Number of splices: Total": star_log.loc[("UNIQUE READS", "Number of splices: Total")],
        "Uniquely mapped reads %": "{:0.2f}".format(star_log.loc[("UNIQUE READS", "Uniquely mapped reads %")]),
        "Uniquely mapped reads number": star_log.loc[("UNIQUE READS", "Uniquely mapped reads number")],

        "Number of reads mapped to multiple loci": star_log.loc[("MULTI-MAPPING READS", "Number of reads mapped to multiple loci")],
        "% of reads mapped to multiple loci": "{:0.2f}%".format(star_log.loc[("MULTI-MAPPING READS", "% of reads mapped to multiple loci")]),
        "Number of reads mapped to too many loci": star_log.loc[("MULTI-MAPPING READS", "Number of reads mapped to too many loci")],
        "% of reads mapped to too many loci": "{:0.2f}%".format(star_log.loc[("MULTI-MAPPING READS", "% of reads mapped to too many loci")]),

        "% of reads unmapped: too many mismatches": "{:0.2f}%".format(star_log.loc[("UNMAPPED READS", "% of reads unmapped: too many mismatches")]),
        "% of reads unmapped: other": "{:0.2f}%".format(star_log.loc[("UNMAPPED READS", "% of reads unmapped: other")]),
        "% of reads unmapped: too short": "{:0.2f}%".format(star_log.loc[("UNMAPPED READS", "% of reads unmapped: too short")]),

        "% of chimeric reads": "{:0.2f}%".format(star_log.loc[("CHIMERIC READS", "% of chimeric reads")]),
        "Number of chimeric reads": star_log.loc[("CHIMERIC READS", "Number of chimeric reads")],

        # run parameters not from the log file
        "quality_metric_of": metric_of,
        "step_run": config["alignment_step_run"],
        "lab": config["lab"],
        "award": config["award"],
    }
    return star_quality_metric


def prepare_star_solo_qc(config, metric_of, filename, umi_plot):
    summary = load_star_solo_quality_metric(filename)

    star_solo_metrics = {
        "assay_term_name": "single-cell RNA sequencing assay",
        "barcode_rank_plot": make_attachment(umi_plot),
        "estimated_number_of_cells": summary["Estimated Number of Cells"],
        "fraction_of_unique_reads_in_cells": summary["Fraction of Unique Reads in Cells"],
        "mean_UMI_per_cell": summary["Mean UMI per Cell"],
        "mean_reads_per_cell": summary["Mean Reads per Cell"],
        "median_UMI_per_cell": summary["Median UMI per Cell"],
        "median_reads_per_cell": summary["Median Reads per Cell"],
        "number_of_reads": summary["Number of Reads"],
        "q30_bases_in_CB_UMI": summary["Q30 Bases in CB+UMI"],
        "q30_bases_in_rna_read": summary["Q30 Bases in RNA read"],
        "reads_mapped_to_genome_unique": summary["Reads Mapped to Genome: Unique"],
        "reads_mapped_to_genome_unique_and_multiple": summary["Reads Mapped to Genome: Unique+Multiple"],
        "reads_with_valid_barcodes": summary["Reads With Valid Barcodes"],
        "sequencing_saturation": summary["Sequencing Saturation"],
        "umis_in_cells": summary["UMIs in Cells"],

        # run parameters not from the log file
        "quality_metric_of": metric_of,
        "step_run": config["quantification_step_run"],
        "lab": config["lab"],
        "award": config["award"],
    }

    if "Mean Gene per Cell" in summary:
        star_solo_metrics.update({
            "mode": "Gene",
            "reads_mapped_to_gene_unique_and_multiple_gene": summary["Reads Mapped to Gene: Unique+Multiple"],
            "reads_mapped_to_gene_unique_gene": summary["Reads Mapped to Gene: Unique"],
            "unique_reads_in_cells_mapped_to_gene": summary["Unique Reads in Cells Mapped to Gene"],
            "mean_gene_per_cell": summary["Mean Gene per Cell"],
            "median_gene_per_cell": summary["Median Gene per Cell"],
            "total_gene_detected": summary["Total Gene Detected"],
        })
    elif "Mean GeneFull_Ex50pAS per Cell" in summary:
        star_solo_metrics.update({
            "mode": "GeneFull_Ex50pAS",
            "reads_mapped_to_genefull_ex50pas_unique_and_multiple_gene_ex50pas": summary["Reads Mapped to GeneFull_Ex50pAS: Unique+Multiple GeneFull_Ex50pAS"],
            "reads_mapped_to_genefull_ex50pas_unique_genefull_ex50pas": summary["Reads Mapped to GeneFull_Ex50pAS: Unique GeneFull_Ex50pAS"],
            "unique_reads_in_cells_mapped_to_genefull_ex50pas": summary["Unique Reads in Cells Mapped to GeneFull_Ex50pAS"],
            "mean_genefull_ex50pas_per_cell": summary["Mean GeneFull_Ex50pAS per Cell"],
            "median_genefull_ex50pas_per_cell": summary["Median GeneFull_Ex50pAS per Cell"],
            "total_genefull_ex50pas_detected": summary["Total GeneFull_Ex50pAS Detected"],
        })
    elif "Mean GeneFull_Ex50pAS per Cell" in summary:
        star_solo_metrics.update({
            "mode": "GeneFull",
            "reads_mapped_to_genefull_unique_and_multiple_genefull": summary["Reads Mapped to GeneFull: Unique+Multiple GeneFull"],
            "reads_mapped_to_genefull_unique_genefull": summary["Reads Mapped to GeneFull: Unique GeneFull"],
            "unique_reads_in_cells_mapped_to_genefull": summary["Unique Reads in Cells Mapped to GeneFull"],
            "mean_genefull_per_cell": summary["Mean GeneFull_Ex50pAS per Cell"],
            "median_genefull_per_cell": summary["Median GeneFull_Ex50pAS per Cell"],
            "total_genefull_detected": summary["Total GeneFull_Ex50pAS Detected"],
        })
    else:
        raise ValueError("Unknown mode? {}".format(sorted(summary.keys())))

    return star_solo_metrics


def prepare_sc_count_matrix_qc_metric(config, metric_of, pct_mt_plot, gene_by_count_plot):
    sc_count_metric = {
        "assay_term_name": "single-cell RNA sequencing assay",
        "total_counts_vs_pct_mitochondria": make_attachment(pct_mt_plot),
        "total_counts_vs_genes_by_count": make_attachment(gene_by_count_plot),

        # run parameters not from the log file
        "quality_metric_of": metric_of,
        "step_run": config["quantification_step_run"],
        "lab": config["lab"],
        "award": config["award"],
    }
    validator.validate(sc_count_metric, "scrna_seq_counts_summary_quality_metric")


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
        "metadata.{}.csv".format(get_submit_host())
    threads: 1
    resources:
        mem_mb = 100
    run:
        metadata = generate_star_solo_processed_metadata(config, {
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
        metadata = "metadata.csv"
    output:
        bam = "Aligned.sortedByCoord.out.bam.upload",
        gene_unique_filtered = "{}_Unique_filtered.tar.gz.{}.upload".format(get_gene_model(), get_submit_host()),
        gene_multi_filtered = "{}_EM_filtered.tar.gz.{}.upload".format(get_gene_model(), get_submit_host()),
        gene_unique_raw = "{}_Unique_raw.tar.gz.{}.upload".format(get_gene_model(), get_submit_host()),
        gene_multi_raw = "{}_EM_raw.tar.gz.{}.upload".format(get_gene_model(), get_submit_host()),
        sj_unique_raw = "SJ_Unique_raw.tar.gz.{}.upload".format(get_submit_host()),
        metadata_posted = "posted.{}.csv".format(get_submit_host()),
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


rule post_star_qc:
    input:
        log_final = "Log.final.out",
        metadata_posted = "posted.{}.csv".format(get_submit_host()),
    output:
        "Log.final.out.{}.qc-upload".format(get_submit_host())
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    run:
        host = get_submit_host()
        server = ENCODED(host)
        uploaded = pandas.read_csv(input.metadata_posted)
        accession = uploaded[uploaded["file_format"] == "bam"]["accession"].to_list()
        qc = prepare_star_qc_metric(config, accession, input.log_final)
        validator = DCCValidator(server)
        validator.validate(qc, "/star_quality_metric/")
        results = server.post_json("/star_quality_metric/", qc)
        with open(output[0], "wt") as outstream:
            outstream.write(str(results))


def get_summary_csv_path():
    return str(SOLO_ROOT / get_gene_model() / "Summary.csv")


rule post_star_solo_qc:
    input:
        gene_summary = get_summary_csv_path(),
        umi_plot = UMI_PER_CELL_PLOT_NAME,
        metadata_posted = "posted.{}.csv".format(get_submit_host()),
    output:
        "{}.{}.qc-upload".format(get_summary_csv_path(), get_submit_host())
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    run:
        host = get_submit_host()
        server = ENCODED(host)
        uploaded = pandas.read_csv(input.metadata_posted)
        accession = uploaded[uploaded["file_format"] == "bam"]["accession"].to_list()
        qc = prepare_star_solo_qc(config, accession, input.gene_summary, input.umi_plot)
        validator = DCCValidator(server)
        validator.validate(qc, "/star_quality_metric/")
        results = server.post_json("/star_quality_metric/", qc)
        with open(output[0], "wt") as outstream:
            outstream.write(str(results))


rule post_count_matrix_qc:
    input:
        pct_mt_plot = "pct_count_mt.{gene_model}_{multiread}_{matrix}.png",
        genes_by_count_plot = "n_genes_by_counts.{gene_model}_{multiread}_{matrix}.png",
        metadata_posted = "posted.{}.csv".format(get_submit_host()),
    output:
        "pct_count_mt.{{gene_model}}_{{multiread}}_{{matrix}}.png.{}.qc-upload".format(get_submit_host())
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    run:
        host = get_submit_host()
        server = ENCODED(host)
        uploaded = pandas.read_csv(input.metadata_posted)
        output_type = filename_to_output_type["{}_EM_filtered.tar.gz".format(get_gene_model())]
        accession = uploaded[uploaded["output_type"] == output_type]["accession"].to_list()
        qc = prepare_sc_count_matrix_qc_metric(config, accession, input.pct_mt_plot, input.gene_by_count_plot)
        validator = DCCValidator(server)
        validator.validate(qc, "/scrna_seq_counts_summary_quality_metric/")
        results = server.post_json("/scrna_seq_counts_summary_quality_metric/", qc)
        with open(output[0], "wt") as outstream:
            outstream.write(str(results))


# Should the metadata be in a single file?
# That makes it harder to see what's been processed and pull submitted
# accessions between rules.
