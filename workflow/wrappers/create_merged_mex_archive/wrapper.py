import anndata
from mex_gene_archive.anndata import write_anndata_as_mex_archive
from mex_gene_archive.reader import read_mex_archive_as_anndata
import numpy
import pandas
from pathlib import Path
from matplotlib import pyplot
import seaborn
from woldrnaseq.madqc import compute_all_vs_all_scores


def merge_mex_metadata(collection, accession=""):
    if isinstance(collection, dict):
        collection = collection.values()

    terms_to_skip = {"software_version", "arguments"}
    found_metadata = {}
    library_ids = []
    for adata in collection:
        for key in adata.uns:
            if key == "library_accession":
                library_ids.append(adata.uns[key])
            elif key not in terms_to_skip:
                found_metadata.setdefault(key, set()).add(adata.uns[key])

    merged_metadata = {}
    for key in found_metadata:
        if len(found_metadata[key]) > 1:
            print(
                "warn {} has multiple values: {}".format(
                    key, ",".join(found_metadata[key])
                )
            )
        else:
            merged_metadata[key] = found_metadata[key].pop()

    merged_metadata["software_version"] = "mex_archive"
    merged_metadata["arguments"] = "merge {}".format(",".join(library_ids))
    merged_metadata["library_accession"] = accession

    return merged_metadata


def load_subpools(config, name, gene_id_to_name):
    subpools = {}
    for library in config["library"]:
        mex = Path(library) / name
        adata = read_mex_archive_as_anndata(mex)
        calculate_qc(adata, gene_id_to_name)
        subpools[library] = adata
    return subpools


subpools = [read_mex_archive_as_anndata(f) for f in snakemake.input]

pseudo_bulk = {}
for filename, adata in zip(snakemake.input, subpools):
    filename = Path(filename)
    library_id = filename.parent.name
    pseudo_bulk[library_id] = numpy.asarray(adata.X.sum(axis=0))[0]
pseudo_bulk = pandas.DataFrame(pseudo_bulk)

scores = compute_all_vs_all_scores(pseudo_bulk)

scores["rafa_spearman"].to_csv(snakemake.output.correlations, sep="\t")

fig = pyplot.figure()
ax = fig.add_subplot(1, 1, 1)
seaborn.heatmap(scores["rafa_spearman"], annot=True, ax=ax)
fig.savefig(snakemake.output.correlation_plot)

merged_metadata = merge_mex_metadata(subpools, snakemake.config["biosample"])

merged = anndata.concat(subpools)
merged.uns = merged_metadata
merged.var["gene_symbols"] = subpools[0].var["gene_symbols"]

write_anndata_as_mex_archive(
    merged, snakemake.params.gene_model, snakemake.wildcards.multiread, "raw"
)

filename = Path(
    "{}_{}_{}.tar.gz".format(
        snakemake.params.gene_model, snakemake.wildcards.multiread, "raw"
    )
)

if not filename.exists():
    raise FileNotFoundError("Missing {}".format(filename))

filename.rename(snakemake.output.matrix)
