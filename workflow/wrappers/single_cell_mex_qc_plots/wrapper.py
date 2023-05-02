from pathlib import Path
from woldrnaseq.plots.scrna_matrix_qc import (
    read_mex_archive_as_anndata,
    calculate_qc,
    generate_violin_plot,
    read_gene_info_names,
    make_pct_mt_scatter,
    make_gene_by_count_scatter
)


count_matrix = Path(snakemake.input.count_matrix)
genome_info = Path(snakemake.input.genome_dir) / "geneInfo.tab"

if "library_id" in snakemake.wildcards:
    library_accession = snakemake.wildcards["library_id"]
elif "library_accession" in snakemake.params:
    library_accession = snakemake.wildcards["library_accession"]
elif len(count_matrix.parts) > 1:
    library_accession = count_matrix.parts[-2]
else:
    raise ValueError("Please specify a library id or accession")

title = "Library {}".format(library_accession)

target_dir = Path(snakemake.output.qc_violin).parent

gene_id_to_name = read_gene_info_names(genome_info)

adata = read_mex_archive_as_anndata(count_matrix)

calculate_qc(adata, gene_id_to_name)

name = count_matrix.name[:-len("tar.gz")]

f = generate_violin_plot(adata)
f.savefig(target_dir / "qc_metric_violin.{}png".format(name))

f = make_pct_mt_scatter(adata, title)
f.savefig(target_dir / "pct_count_mt.{}png".format(name))

f = make_gene_by_count_scatter(adata, title)
f.savefig(target_dir / "n_genes_by_counts.{}png".format(name))
