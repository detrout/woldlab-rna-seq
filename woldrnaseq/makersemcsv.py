#!/usr/bin/python3
import argparse
import logging

from woldrnaseq import models
from woldrnaseq import madqc

from woldrnaseq.common import (
    add_debug_arguments,
    configure_logging,
    get_seperator,
)

logger = logging.getLogger('RSEM CSV')


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    sep = get_seperator(args.sep)
    experiments = models.load_experiments(args.experiments, sep=sep)
    libraries = models.load_library_tables(args.libraries, sep=sep)

    output_sep = get_seperator(args.output_format)
    output_extension = {
        'TAB': '.tsv',
        ',': '.csv',
    }[args.output_format]

    if args.transcriptome:
        # isoforms
        load_quantifications = madqc.load_transcriptome_quantifications
        quantification_extension = '_isoform_' + args.quantification + output_extension
    else:
        # genes
        load_quantifications = madqc.load_genomic_quantifications
        quantification_extension = '_gene_' + args.quantification + output_extension

    for name in experiments:
        filename = name + quantification_extension
        replicates = experiments[name]
        logger.info("%s %s: %s",
                    name, args.quantification, ','.join(replicates))
        quantifications = load_quantifications(
            replicates, libraries, args.quantification)
        quantifications.to_csv(filename, sep=output_sep)


def make_parser():
    parser = argparse.ArgumentParser(
        description="Write combined CSV for a specific quantification.")
    parser.add_argument('-q', '--quantification',
                        choices=['FPKM', 'TPM', 'expected_count',
                                 'effective_length', 'length'],
                        default='FPKM',
                        help='which quantification value to use')
    parser.add_argument('-l', '--libraries', action='append', help='library information table')
    parser.add_argument('-e', '--experiments', action='append', help='experiment information table')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    parser.add_argument('--output-format', choices=['TAB', ','], default=',')
    parser.add_argument('--transcriptome', action='store_true', default=False,
                        help='Use RSEM transcriptome quantifications instead of genomic quantifications')
    add_debug_arguments(parser)
    return parser


if __name__ == '__main__':
    main()
