#!/usr/bin/python3
"""Make refFlat file from gff/gtf

Picard CollectRnaSeqMetrics requires a ref_flat to run and there
wasn't an obvious way to generate ref_flat files from gff/gtf files.

See:
https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.0.0/picard_analysis_CollectRnaSeqMetrics.php

"""
from argparse import ArgumentParser
import os
import sys

from woldrnaseq.gff2table import parse_attributes


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.output:
        outstream = open(args.output, 'wt')
    else:
        outstream = sys.stdout

    with open(args.gtf, 'rt') as instream:
        convert_to_refflat(instream, outstream, sep=args.sep)

    if args.output:
        outstream.close()

    return 0


def make_parser():
    parser = ArgumentParser()
    parser.add_argument('gtf')
    parser.add_argument('-o', '--output', help='name to write output to')
    parser.add_argument('--sep', default=' ',
                        help='Separator between fields in 8th column')
    return parser


def convert_to_refflat(instream, outstream, sep):
    exon_starts = []
    exon_stops = []
    start = None
    stop = None
    transcript_id = None
    gene_name = None
    chromosome = None
    strand = None

    for line in instream:
        records = line.split('\t')
        record_type = records[2]
        attributes = dict(parse_attributes(records[8], sep))

        if record_type == 'transcript' or transcript_id != attributes.get('transcript_id', None):
            if len(exon_starts) > 0:
                write_refflat_record(
                    outstream,
                    transcript_id,
                    gene_name,
                    chromosome,
                    strand,
                    start,
                    stop,
                    exon_starts,
                    exon_stops)
                exon_starts = []
                exon_stops = []

        chromosome = records[0]
        start = records[3]
        stop = records[4]
        strand = records[6]
        transcript_id = attributes.get('transcript_id', None)
        gene_name = attributes.get('gene_name', attributes['gene_id'])
        if record_type == 'exon':
            exon_starts.append(start)
            exon_stops.append(stop)

    write_refflat_record(
        outstream,
        transcript_id,
        gene_name,
        chromosome,
        strand,
        start,
        stop,
        exon_starts,
        exon_stops)


def write_refflat_record(outstream, transcript_id, gene_name,
                         chromosome, strand,
                         tx_start, tx_stop,
                         exon_starts, exon_stops):
    cds_start = str(min([int(x) for x in exon_starts]))
    cds_stop = str(max([int(x) for x in exon_stops]))

    outstream.write('\t'.join([
        transcript_id,
        gene_name,
        chromosome,
        strand,
        tx_start,
        tx_stop,
        cds_start,
        cds_stop,
        str(len(exon_starts)),
        ','.join(exon_starts) + ',',
        ','.join(exon_stops) + ',']))
    outstream.write(os.linesep)


if __name__ == '__main__':
    sys.exit(main())
