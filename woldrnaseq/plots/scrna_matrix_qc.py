#!/usr/bin/python3
from argparse import ArgumentParser
import matplotlib
from matplotlib import pyplot
from mex_gene_archive.reader import read_mex_archive_as_anndata
from pathlib import Path
import re
import scanpy as sc
import sys

matplotlib.use("Agg")

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)
    target_dir = Path(args.output_dir)

    gene_id_to_name = read_gene_info_names(args.gene_info)
    adata = read_mex_archive_as_anndata(args.filename[0])

    calculate_qc(adata, gene_id_to_name)

    f = generate_violin_plot(adata)
    f.savefig(target_dir / "qc_metric_violin.png")

    f = make_pct_mt_scatter(adata, args.title)
    f.savefig(target_dir / "pct_count_mt.png")

    f = make_gene_by_count_scatter(adata, args.title)
    f.savefig(target_dir / "n_genes_by_counts.png")

    return 0


def read_gene_info_names(gene_info):
    gene_id_name = {}
    with open(gene_info, "rt") as instream:
        _ = next(instream).rstrip()
        for row in instream:
            fields = row.rstrip().split("\t")
            if len(fields) < 3:
                raise RuntimeError(
                    "The geneInfo file does not contain gene names please use a newer version of STAR"
                )
            gene_id_name[fields[0]] = fields[1]

    return gene_id_name


def calculate_qc(adata, gene_id_to_name):
    mt_prefix = re.compile("^mt[-_]", re.I)
    adata.var["gene_name"] = [
        gene_id_to_name.get(x, x) for x in adata.var_names
    ]
    adata.var["mt"] = adata.var["gene_name"].apply(lambda x: mt_prefix.match(x) is not None)
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )


def generate_violin_plot(adata):
    f = sc.pl.violin(
        adata,
        ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    return f


def make_pct_mt_scatter(adata, title):
    f = pyplot.figure()
    ax = f.add_subplot(1, 1, 1)
    sc.pl.scatter(
        adata, x="total_counts", y="pct_counts_mt", title=title, ax=ax, show=False
    )
    return f


def make_gene_by_count_scatter(adata, title):
    f = pyplot.figure()
    ax = f.add_subplot(1, 1, 1)
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        title=title,
        ax=ax,
        show=False,
    )
    return f


def make_parser():
    parser = ArgumentParser()
    parser.add_argument("filename", nargs=1)
    parser.add_argument(
        "--gene-info", required=True, help="Specify path to STAR geneInfo"
    )
    parser.add_argument(
        "-o", "--output-dir", default=".", help="output directory to write plots too"
    )
    parser.add_argument("--title", default=None, help="set title for plots")
    return parser


if __name__ == "__main__":
    sys.exit(main())
