#!/usr/bin/python3
"""Make ribosomal RNA interval list for Picard CollectRnaSeqMetrics

The format appears to be a sam header followed by
chr start stop strand transcript_id
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

    try:
        with open(args.genome, 'rt') as instream:
            generate_sam_header_from_fasta(instream, outstream)
        with open(args.gtf, 'rt') as instream:
            generate_rrna_regions(instream, outstream)
    except Exception:
        if args.output:
            outstream.close()
            return 1
    return 0


def make_parser():
    parser = ArgumentParser()
    parser.add_argument('--genome', required=True,
                        help='Reference genome to use from reference lengths')
    parser.add_argument('--gtf', required=True,
                        help='GTF File to scan')
    parser.add_argument('-o', '--output', help='Destination filename')
    return parser


def generate_sam_header_from_fasta(instream, outstream):
    length = 0
    name = ""
    for line in instream:
        if line.startswith('>'):
            # Header
            if length > 0:
                write_header(outstream, name, length)
                length = 0
            name = line[1:].rstrip()
        else:
            length += len(line.rstrip())

    write_header(outstream, name, length)


def write_header(stream, name, length):
    stream.write("@SQ\tSN:{name}\tLN:{length}{eol}".format(
        name=name,
        length=length,
        eol=os.linesep,
    ))


def generate_rrna_regions(instream, outstream):
    for line in instream:
        records = line.split('\t')
        attributes = dict(parse_attributes(records[8], sep=' '))
        transcript_id = attributes['transcript_id']
        gene_type = attributes.get('gene_type', [None])
        if gene_type == 'rRNA':
            outstream.write("\t".join((records[0], records[3], records[4], records[6], transcript_id)))
            outstream.write(os.linesep)


if __name__ == '__main__':
    sys.exit(main())
