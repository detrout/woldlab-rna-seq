import sys
import logging
import pandas
from encoded_client.encoded import ENCODED, DCCValidator, make_attachment


filename_to_output_type = {
    "Aligned.sortedByCoord.out.bam": "alignments",
    "Gene_Unique_filtered.tar.gz": "sparse gene count matrix of unique reads",
    "Gene_EM_filtered.tar.gz": "sparse gene count matrix of all reads",
    "Gene_Unique_raw.tar.gz": "unfiltered sparse gene count matrix of unique reads",
    "Gene_EM_raw.tar.gz": "unfiltered sparse gene count matrix of all reads",
    "GeneFull_Ex50pAS_Unique_filtered.tar.gz": "sparse gene count matrix of unique reads",
    "GeneFull_Ex50pAS_EM_filtered.tar.gz": "sparse gene count matrix of all reads",
    "GeneFull_Ex50pAS_Unique_raw.tar.gz": "unfiltered sparse gene count matrix of unique reads",
    "GeneFull_Ex50pAS_EM_raw.tar.gz": "unfiltered sparse gene count matrix of all reads",
    "SJ_Unique_raw.tar.gz": "unfiltered sparse splice junction count matrix of unique reads",
}


def prepare_sc_count_matrix_qc_metric(config, metric_of, pct_mt_plot, gene_by_count_plot, genes_by_count_plot):
    sc_count_metric = {
        "assay_term_name": "single-cell RNA sequencing assay",
        "total_counts_vs_pct_mitochondria": make_attachment(pct_mt_plot),
        "total_counts_vs_genes_by_count": make_attachment(gene_by_count_plot),
        "counts_violin_plot": make_attachment(genes_by_count_plot),

        # run parameters not from the log file
        "quality_metric_of": metric_of,
        "step_run": config["quantification_step_run"],
        "lab": config["lab"],
        "award": config["award"],
    }
    return sc_count_metric


rule generate_umi:
    input:
        #rules.star_solo_10x.output.umis
        SOLO_ROOT / get_gene_model() / "UMIperCellSorted.txt",
    output:
        UMI_PER_CELL_PLOT_NAME
    resources:
        mem_mb = DEFAULT_MEM_MB
    singularity:
        config['scanpy_container']
    shell:
        """python3 -m woldrnaseq.plots.star_solo_barcodes \
                    --umi-per-cell {input} -o {output}
        """


rule generate_count_qc_metrics:
    input:
        count_matrix = "{gene_model}_{multiread}_{matrix}.tar.gz",
        genome_dir = config['genome_dir'],
    output:
        qc_violin = "qc_metric_violin.{gene_model}_{multiread}_{matrix}.png",
        pct_mt = "pct_count_mt.{gene_model}_{multiread}_{matrix}.png",
        genes_by_count = "n_genes_by_counts.{gene_model}_{multiread}_{matrix}.png",
    singularity:
        config['scanpy_container']
    params:
        library_accession = lambda wildcards: getattr(wildcards["library_id"], config["library_accession"]),
    resources:
        mem_mb = config["mem_mb"]
    shell:
        """python3 -m woldrnaseq.plots.scrna_matrix_qc \
                     --gene-info {input.genome_dir}/geneInfo.tab \
                     --title "Library {params.library_accession}" \
                     {input.count_matrix}
        """

rule post_count_matrix_qc:
    input:
        pct_mt_plot = "pct_count_mt.{gene_model}_{multiread}_{matrix}.png",
        genes_by_count_plot = "n_genes_by_counts.{gene_model}_{multiread}_{matrix}.png",
        counts_violin_plot = "qc_metric_violin.{gene_model}_{multiread}_{matrix}.png",
        metadata_posted = "posted.{}.csv".format(get_submit_host(config)),
    output:
        "pct_count_mt.{{gene_model}}_{{multiread}}_{{matrix}}.png.{}.qc-upload".format(get_submit_host(config))
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB
    log: "post_count_matrix_qc_{gene_model}_{multiread}_{matrix}.log"
    run:
        logger = logging.getLogger("post_count_matrix_qc")
        logger.setLevel(logging.INFO)
        logger.addHandler(logging.FileHandler(log[0]))
        logger.addHandler(logging.StreamHandler(sys.stderr))
        host = get_submit_host(config)
        server = ENCODED(host)
        uploaded = pandas.read_csv(input.metadata_posted)
        file_qced = "{}_EM_{}.tar.gz".format(get_gene_model(), wildcards.matrix)
        logger.info("Using {} as QC file".format(file_qced))
        output_type = filename_to_output_type[file_qced]
        accession = uploaded[uploaded["output_type"] == output_type]["accession"].to_list()
        assert len(accession) > 0 and not pandas.isnull(accession)
        qc = prepare_sc_count_matrix_qc_metric(
            config,
            accession,
            input.pct_mt_plot,
            input.genes_by_count_plot,
            input.counts_violin_plot
        )
        try:
            validator = DCCValidator(server)
            validator.validate(qc, "scrna_seq_counts_summary_quality_metric")
            results = server.post_json("/scrna_seq_counts_summary_quality_metric/", qc)
            with open(output[0], "wt") as outstream:
                outstream.write(str(results))
        except Exception as e:
            logger.error(e)
            raise e
