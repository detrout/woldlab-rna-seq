#
# This contains rules to archive the STAR Solo output files
# into mex_gene_archive files.
#

from mex_gene_archive.filter import filter_mtx
from mex_gene_archive.starsolo import (
    archive_star_solo,
    make_list_of_archive_files,
)


rule filter_em_matrix:
    #
    # The version of STAR Solo I was testing was not writing out the smaller filtered
    # EM count matrices since was supposed to be relatively easy to compute from
    # the raw matrix and the filtered barcodes
    #
    # Doing it efficiently and in a way that is binary identical (for the
    # unique matrix being generated in both cases) took a bit.
    #
    input:
        filtered_barcode_tsv = SOLO_ROOT / get_gene_model() / "filtered" / "barcodes.tsv",
        raw_barcode_tsv = SOLO_ROOT / get_gene_model() / "raw" / "barcodes.tsv",
        raw_em_matrix_mtx = SOLO_ROOT / get_gene_model() / "raw" / "UniqueAndMult-EM.mtx",
    output:
        filtered_em_mtx = SOLO_ROOT / get_gene_model() / "filtered" / "UniqueAndMult-EM.mtx",
    threads: 1
    resources:
        mem_mb = config.get("mem_mb", DEFAULT_MEM_MB)
    run:
        with open(output.filtered_em_mtx, "wt") as outstream:
            for line in filter_mtx(input.raw_barcode_tsv,
                                   input.raw_em_matrix_mtx,
                                   input.filtered_barcode_tsv):
                outstream.write(line)


rule to_archive:
    #
    # This rule generates the various
    input:
        # At least snakemake 5.4.0 wants input files to be a list of strings
        # make_list_of_archive_files returns a list of Paths so we need to
        # convert them
        lambda wildcards: [str(x) for x in make_list_of_archive_files(
            SOLO_ROOT, wildcards.gene_model, wildcards.multiread, wildcards.matrix)]
    output:
        "{gene_model}_{multiread}_{matrix}.tar.gz"
    wildcard_constraints:
        gene_model = "(GeneFull_Ex50pAS)|(GeneFull)|(Gene)|(SJ)",
        multiread = "(EM)|(Unique)",
        matrix = "(raw)|(filtered)",
    threads: 1
    resources:
        mem_mb = DEFAULT_MEM_MB,
    run:
        archive_star_solo(
            SOLO_ROOT,
            config,
            wildcards.gene_model,
            wildcards.multiread,
            wildcards.matrix
        )
