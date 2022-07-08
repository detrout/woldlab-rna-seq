#!/usr/bin/python3
import argparse
import logging

from woldrnaseq import models
from woldrnaseq import madqc
from woldrnaseq.gtfcache import GTFCache

from woldrnaseq.common import (
    add_debug_arguments,
    configure_logging,
    get_seperator,
    sanitize_name,
)

logger = logging.getLogger('RSEM CSV')


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    sep = get_seperator(args.sep)
    experiments = models.load_experiments(args.experiments, sep=sep)
    libraries = models.load_library_tables(args.libraries, sep=sep)


    gtf_cache = None
    if args.add_names:
        if args.genome_dir is None:
            parser.error('genome-dir is needed to add names to the quantification file')
        else:
            gtf_cache = GTFCache(libraries, args.genome_dir)

    if len(args.quantification) > 0:
        quantification_list = args.quantification
    else:
        quantification_list = ['FPKM']

    if args.transcriptome:
        # isoforms
        RsemLoader = IsoformRsemLoader
    else:
        # genes
        RsemLoader = GeneRsemLoader

    for quantification in quantification_list:
        logger.info('Building expression matrix for %s', quantification)
        for i, experiment in experiments.iterrows():
            loader = RsemLoader(quantification, gtf_cache)
            matrix = loader.load(experiment, libraries)
            loader.save(matrix, args.output_format)


def make_parser():
    parser = argparse.ArgumentParser(
        description="Write combined CSV for a specific quantification produced by RSEM.")
    parser.add_argument('-q', '--quantification',
                        choices=['FPKM', 'TPM', 'expected_count',
                                 'effective_length', 'length'],
                        default=[],
                        action='append',
                        help='which quantification value to use')
    parser.add_argument('-l', '--libraries', action='append', required=True,
                        help='library metadata table')
    parser.add_argument('-e', '--experiments', action='append', required=True,
                        help='experiment metadata table')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    parser.add_argument('--output-format', choices=['TAB', ','], default=',')
    parser.add_argument('--transcriptome', action='store_true', default=False,
                        help='Use RSEM transcriptome quantifications instead of genomic quantifications')
    parser.add_argument('--genome-dir', help='Specify genome directory root')
    parser.add_argument('--add-names', action='store_true', default=False,
                        help='Add names to ouptut quantification file, requires --gtf-cache')
    add_debug_arguments(parser)
    return parser


class RsemLoader:
    """Help generate the combined quantification files

    Parameters
    ----------
    quantification_name: str
        the name of the quantification, such as, FPKM, TPM, effective_count.
    annotation: pandas.DataFrame
        A GTF cache file to help map gene_ids to gene_names.
    """
    def __init__(self, quantification_name, annotation=None):
        self.quantification_name = quantification_name
        self.annotation = annotation

    def save(self, quantifications, output_format, filename=None):
        output_sep = get_seperator(output_format)
        output_extension = {
            'TAB': '.tsv',
            ',': '.csv',
        }[output_format]

        if filename is None:
            filename = sanitize_name(quantifications.name) + output_extension

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
                self.annotation[replicates[0]], quantifications)

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
                self.annotation[replicates[0]], quantifications)

        quantifications.name = experiment.name + '_gene_' + self.quantification_name
        return quantifications


if __name__ == '__main__':
    main()
