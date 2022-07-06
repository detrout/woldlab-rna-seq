#!/usr/bin/python3

"""Attempt to determine how well reads align to the RNA models.

For RNA-seq you expect most of the reads to line up with gene models,
for fully spliced messages most of the reads should align with exons,
for less processed transcripts there should be more reads cross into
introns
"""

from argparse import ArgumentParser
from collections import namedtuple
import logging
import pandas
import pysam
import sys
import time

from .common import (
    add_debug_arguments,
    configure_logging,
)
from .gff2table import tokenize as gff_tokenize

logger = logging.getLogger(__name__)

GENE = "G"
EXON = "E"

DistributionTotals = namedtuple(
    "DistributionTotals",
    ["Exonic", "Intronic", "Intergenic", "Spikeins"])

ReferenceDistribution = namedtuple(
    "ReferenceDistribution",
    ["Exonic", "Intronic", "Intergenic"])


def count_exonic_genomic_reads(bam, gtf_name, stranded="unstranded"):
    assert stranded in ["forward", "reverse", "unstranded"]
    if stranded == "unstranded":
        check_strand = False
    else:
        check_strand = True

    with pysam.AlignmentFile(bam, "rb") as bam:
        exon_read_count = 0
        gene_read_count = 0
        intergenic_read_count = 0
        other_read_count = 0
        spikein_read_count = 0
        t0 = time.monotonic()

        for reference_name in bam.references:
            gtf_cache = read_tabix_as_pandas(gtf_name, reference_name)
            gene_plus, gene_minus = build_gene_locations(
                gtf_cache, reference_name, check_strand=check_strand
            )
            counts = count_exonic_genomic_reads_for_reference(
                bam, reference_name, gene_plus, gene_minus, stranded
            )
            if not is_spike(reference_name):
                exon_read_count += counts["Exonic"]
                gene_read_count += counts["Intronic"]
                intergenic_read_count += counts["Intergenic"]
                other_read_count += counts["Other"]
            else:
                spikein_read_count += sum(counts.values())

        logger.info(
            "exons {}, genes {}, intergenic count {}, other{}, spikes {} in {} seconds".format(
                exon_read_count,
                gene_read_count,
                intergenic_read_count,
                other_read_count,
                spikein_read_count,
                time.monotonic() - t0,
            )
        )

    return DistributionTotals(
        Exonic=exon_read_count,
        Intronic=gene_read_count,
        Intergenic=intergenic_read_count,
        Spikeins=spikein_read_count,
    )


def count_exonic_genomic_reads_for_reference(
    bam, reference_name, gene_plus, gene_minus, stranded
):
    assert stranded in ["forward", "reverse", "unstranded"]
    exon_read_count = 0
    gene_read_count = 0
    intergenic_read_count = 0
    other_read_count = 0
    t0 = time.monotonic()

    # by using 3 values in the tuple we can take the -1 / 1 flag we
    # create from reverse add 1 to it and use it to lookup the right
    # dictionary in this table
    #          strand     -1        0       1
    gene_map = {
        "forward": (gene_minus, None, gene_plus),
        "unstranded": (gene_plus, None, gene_plus),
        "reverse": (gene_plus, None, gene_minus),
    }
    print("+", gene_plus)
    print("-", gene_minus)
    for read in bam.fetch(reference_name):
        strand = -1 if read.is_reverse else 1
        print(
            "read: {} {} {} {} {} ".format(
                read.query_name, read.reference_name, read.pos, strand, read.cigarstring
            )
        )
        exon_base_count = 0
        gene_base_count = 0
        base_count = 0

        gene_model = gene_map[stranded][strand + 1]
        print("reference pos:", read.get_reference_positions())
        for pos in read.get_reference_positions():
            base_count += 1
            code = gene_model.get(pos, None)
            if code == "E":
                exon_base_count += 1
            elif code == "G":
                gene_base_count += 1

        print("ega", exon_base_count, gene_base_count, base_count)
        weight = 1 / read.get_tag("NH")
        if exon_base_count / base_count >= 0.5:
            exon_read_count += weight
        elif gene_base_count / base_count >= 0.5:
            gene_read_count += weight
        else:
            intergenic_read_count += weight

    tnow = time.monotonic()
    print(
        "{}: {} {} {} in {:.4}s".format(
            reference_name,
            exon_read_count,
            gene_read_count,
            intergenic_read_count,
            other_read_count,
            tnow - t0,
        )
    )
    return ReferenceDistribution(
        Exonic=exon_read_count,
        Intronic=gene_read_count,
        Intergenic=intergenic_read_count,
    )


def build_gene_locations(gtf, reference_name, check_strand=True):
    """construct dictionaries of what type various genomic locations are.

    G means is within the gene model.
    E means is within an exon in the gene model.

    returns a plus and minus strand dictionary in 0-based coordinates

    If check_strand is False all genes will be in the plus dictionary

    For chromsome 1 of the gencode GRCh38 comprehensive set this takes a bit
    over 5 gigabytes of memory.
    """
    gene_plus = {}
    gene_minus = {}

    for i, row in gtf.iterrows():
        if check_strand and row.strand == -1:
            gene_locations = gene_minus
        else:
            gene_locations = gene_plus

        if row["type"] == "gene":
            for pos in range(row["start"] - 1, row["stop"]):
                gene_locations.setdefault(pos, GENE)

        if row["type"] == "exon":
            for pos in range(row["start"] - 1, row["stop"]):
                gene_locations[pos] = EXON

    logger.debug(
        "{} plus={}, minus={}".format(reference_name, len(gene_plus), len(gene_minus))
    )
    return (gene_plus, gene_minus)


def read_filtered_gtf_cache(gtf_name, reference_name):
    """Read just the relevant parts of the gtf cache"""
    logger.info("Reading cache {}".format(gtf_name))
    store = pandas.HDFStore(gtf_name)
    gtf = store.select(
        "gtf",
        where='chromosome="{}"'.format(reference_name),
        columns=["chromosome", "type", "start", "stop", "strand", "gene_id"],
    )
    store.close()
    logger.debug("gtf shape for {}: {}".format(reference_name, gtf.shape))
    return gtf
def read_tabix(gtf_name, reference_name, sep=' '):
    logger.info("Reading gtf file {}".format(gtf_name))
    strand_map = {
        "+": 1,
        ".": 0,
        "-": -1,
    }
    with pysam.TabixFile(gtf_name) as tab:
        if reference_name not in tab.contigs:
            logger.info("{} is not in the GTF/GFF file".format(reference_name))
            return []
        else:
            for read in tab.fetch(reference_name):
                fields = read.split("\t")
                chromosome = fields[0]
                record_type = fields[2]
                start = int(fields[3])
                stop = int(fields[4])
                strand = strand_map[fields[6]]
                tokens = gff_tokenize(fields[8], sep)
                gene_id = None
                gene_name = None
                for term in tokens:
                    name = term
                    attribute_sep = next(tokens)
                    assert attribute_sep == sep
                    value = next(tokens)
                    if name == "gene_id":
                        gene_id = value
                    elif name == "gene_name":
                        gene_name = value

                    if gene_id is not None and gene_name is not None:
                        break

                    try:
                        field_sep = next(tokens)
                        assert field_sep == ";"
                    except StopIteration:
                        break

                yield (chromosome, start, stop, strand, record_type, gene_id, gene_name)


def read_tabix_as_pandas(gtf_name, reference_name, sep=' '):
    return pandas.DataFrame(
        read_tabix(gtf_name, reference_name, sep=sep),
        columns=["chromosome", "start", "stop", "strand", "type", "gene_id", "gene_name"])


def is_spike(reference_name):
    if reference_name.startswith("ERCC"):
        return True
    elif reference_name == "phiX174":
        return True
    else:
        return False


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    counts = count_exonic_genomic_reads(args.filename[0], args.gtf_cache, stranded=args.strand)

    counts = pandas.Series(counts)
    counts.index.name = "#Class"
    counts.name = "Counts"
    fractions = counts / counts.sum()
    fractions.name = "Fraction"

    try:
        if args.output is not None:
            outstream = open(args.output, "wt")
        else:
            outstream = sys.stdout
        fractions.to_csv(outstream, sep="\t")
    finally:
        if args.output is not None:
            outstream.close()

    return 0


def make_parser():
    parser = ArgumentParser()

    parser.add_argument("filename", nargs=1, help="name of bam file to read")
    parser.add_argument(
        "--gtf-cache", required=True, help="Location of .h5 cache of parse gtf file"
    )
    parser.add_argument(
        "--strand", choices=["forward", "unstranded", "reverse"], default="unstranded"
    )
    parser.add_argument("-o", "--output", help="Location to write result file")
    add_debug_arguments(parser)
    return parser


if __name__ == "__main__":
    sys.exit(main())
