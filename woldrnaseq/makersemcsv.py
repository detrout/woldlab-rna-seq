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

    if args.add_names:
        if args.gtf_cache is None:
            parser.error('GTF-cache is needed to add names to the quantification file')
        else:
            logger.info('Loading GTF Cache %s', args.gtf_cache)
            annotation = models.load_gtf_cache(args.gtf_cache)
    else:
        annotation = None

    if args.transcriptome:
        # isoforms
        loader = IsoformRsemLoader(args.quantification, annotation)
    else:
        # genes
        loader = GeneRsemLoader(args.quantification, annotation)

    for i, experiment in experiments.iterrows():
        quantification = loader.load(experiment, libraries)
        loader.save(quantification, args.output_format)


def make_parser():
    parser = argparse.ArgumentParser(
        description="Write combined CSV for a specific quantification.")
    parser.add_argument('-q', '--quantification',
                        choices=['FPKM', 'TPM', 'expected_count',
                                 'effective_length', 'length'],
                        default='FPKM',
                        help='which quantification value to use')
    parser.add_argument('-l', '--libraries', action='append', required=True,
                        help='library information table')
    parser.add_argument('-e', '--experiments', action='append', required=True,
                        help='experiment information table')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    parser.add_argument('--output-format', choices=['TAB', ','], default=',')
    parser.add_argument('--transcriptome', action='store_true', default=False,
                        help='Use RSEM transcriptome quantifications instead of genomic quantifications')
    parser.add_argument('--gtf-cache', help='Specify name of GTF cache file')
    parser.add_argument('--add-names', action='store_true', default=False,
                        help='Add names to ouptut quantification file')
    add_debug_arguments(parser)
    return parser


class RsemLoader:
    def __init__(self, quantification_name, annotation):
        self.quantification_name = quantification_name
        self.annotation = annotation

    def save(self, quantifications, output_format, filename=None):
        output_sep = get_seperator(output_format)
        output_extension = {
            'TAB': '.tsv',
            ',': '.csv',
        }[output_format]

        if filename is None:
            filename = quantifications.name + output_extension

        quantifications.to_csv(filename, sep=output_sep)


class IsoformRsemLoader(RsemLoader):
    def load(self, experiment, libraries):
        replicates = experiment['replicates']
        logger.info("%s %s: %s",
                    experiment.name, self.quantification_name, ','.join(replicates))
        quantifications = madqc.load_transcriptome_quantifications(
            experiment, libraries, self.quantification_name)

        if self.annotation is not None:
            quantifications = models.lookup_gene_name_by_transcript_id(
                self.annotation, quantifications)

        quantifications.name = experiment.name + '_isoform_' + self.quantification_name
        return quantifications


class GeneRsemLoader(RsemLoader):
    def load(self, experiment, libraries):
        replicates = experiment['replicates']
        logger.info("%s %s: %s",
                    experiment.name, self.quantification_name, ','.join(replicates))
        quantifications = madqc.load_genomic_quantifications(
            experiment, libraries, self.quantification_name)

        if self.annotation is not None:
            quantifications = models.lookup_gene_name_by_gene_id(
                self.annotation, quantifications)

        quantifications.name = experiment.name + '_gene_' + self.quantification_name
        return quantifications


if __name__ == '__main__':
    main()
