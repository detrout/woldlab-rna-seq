
include:
    "star_solo_10x_config.smk"

rule test_all:
    input:
        "umi_per_cell.png",
        "pct_count_mt.png",
        "n_genes_by_counts.png",

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
        library_accession = config["library_accession"],
    resources:
        mem_mb = config["mem_mb"]
    shell:
        """python3 -m woldrnaseq.plots.scrna_matrix_qc \
                     --gene-info {input.genome_dir}/geneInfo.tab \
                     --title "Library {params.library_accession}" \
                     {input.count_matrix}
        """
