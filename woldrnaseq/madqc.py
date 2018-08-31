import argparse
import collections
from glob import glob
import pandas
import numpy
import os
import scipy.stats
import sys
import logging

from woldrnaseq import models
from .common import get_seperator

logger = logging.getLogger('madQC')


def load_rsem_quantifications(experiment_files, index=None, column='FPKM'):
    """Load quantifications out of RSEM results into a pandas dataframe

    Columns will be library accession identifiers.

    Parameters
    ----------
    experiment_files: list of filenames
        list of paths to RSEM result files.
    index: list of column names
    column: str
        what RSEM column to load.

    Returns
    -------
    pandas.DataFrame
        gene id (rows) x library id (columns) table of the
        specified RSEM column
    """
    quantifications = []
    filenames = []
    for filename in experiment_files:
        try:
            table = pandas.read_csv(filename,
                                    sep='\t',
                                    index_col=0)
        except pandas.io.common.EmptyDataError as e:
            logger.error("Unable to read file {} empty data".format(filename))
            sys.exit(1)

        table = table[[column]]

        quantifications.append(table[column])
        _, name = os.path.split(filename)
        filenames.append(name)

    df = pandas.concat(quantifications, axis=1)
    if index:
        df.columns = index
    else:
        df.columns = filenames
    return df


def replicate_scores(table, rep1_name, rep2_name, Acutoff=0):
    """Compute correlations, MAD, and SD replicate comparison scores

    Parameters
    ----------
    table: pandas.DataFrame
       DataFrame of RSEM quantifications
    rep1_name: str
       name of one column to score (library name)
    rep2_name: str
       name of the other column to score (library name)
    Acutoff: float
       Require that the average of the two replicates is greater
       than this cutoff in log2 space

    Returns
    -------
    pandas.Series
       Containing all the scores between the two replicates
    """
    rep1 = table[rep1_name]
    rep2 = table[rep2_name]

    eitherzero = (rep1 == 0) | (rep2 == 0)
    replz1 = numpy.log2(rep1[eitherzero != True])
    replz2 = numpy.log2(rep2[eitherzero != True])

    M = replz1 - replz2
    A = (replz1 + replz2) / 2.0

    replz1_gt_Acutoff = replz1[A > Acutoff]
    replz2_gt_Acutoff = replz2[A > Acutoff]

    if len(replz1_gt_Acutoff) == 0:
        logger.warning("No data survived Acutoff filter in %s", rep1_name)
    if len(replz2_gt_Acutoff) == 0:
        logger.warning("No data survived Acutoff filter in %s", rep2_name)

    scores = pandas.Series({
        'total_rows': len(table),
        'passed_filter': len(replz1[A > Acutoff]),

        'naive_pearson': scipy.stats.pearsonr(rep1, rep2)[0],
        'naive_spearman': scipy.stats.spearmanr(rep1, rep2)[0],

        'rafa_pearson': scipy.stats.pearsonr(replz1_gt_Acutoff,
                                             replz2_gt_Acutoff)[0],
        'rafa_spearman': scipy.stats.spearmanr(replz1_gt_Acutoff,
                                               replz2_gt_Acutoff)[0],
        'MAD': numpy.round(numpy.median(numpy.abs(M)[A > Acutoff]) * 1.4826, 3),
        'SD': numpy.round(numpy.sqrt(numpy.mean(M[A > Acutoff] ** 2)), 3)
    },
        index=['total_rows', 'passed_filter',
               'naive_pearson', 'naive_spearman',
               'rafa_pearson', 'rafa_spearman',
               'MAD', 'SD']
    )
    return scores


def compute_all_vs_all_scores(fpkms, Acutoff=0):
    """Compute all the scores of note for a FPKM table.

    Parameters
    ----------
    fpkms: pandas.DataFrame
        Table of gene x library quantification scores

    Returns
    -------
    pandas.Panel
        Containing every replicate vs every other replicate by 
        the pandas.Series of scores from replicate_scores
    """
    all_scores = collections.OrderedDict()
    shape = (len(fpkms.columns), len(fpkms.columns))
    for i in range(0, len(fpkms.columns)):
        for j in range(0, len(fpkms.columns)):
            if i == j:
                continue

            rep1 = fpkms.columns[i]
            rep2 = fpkms.columns[j]
            scores = replicate_scores(fpkms, rep1, rep2, Acutoff)

            for name in scores.keys():
                if name not in all_scores:
                    a = numpy.empty(shape)
                    a.fill(numpy.nan)
                    all_scores[name] = pandas.DataFrame(
                        a,
                        index=fpkms.columns,
                        columns=fpkms.columns
                    )
                all_scores[name][rep1][rep2] = scores[name]
    return pandas.Panel(all_scores)


def load_genomic_quantifications(experiment, libraries, column='FPKM'):
    """Load all of the RSEM gene quantifications files for an experiment

    Parameters
    ----------
    experiment: pandas.Series
        row of information about an experiment
    libraries: pandas.DataFrame
        Table of library metadata
    columns: str
       name of quantification column to load

    Returns
    -------
    pandas.DataFrame
       Table of all the quantification folumn for all gene RSEM result files.
    """
    extension = '*_rsem.genes.results'
    return load_rsem_replicates(extension, experiment, libraries, column)


def load_transcriptome_quantifications(experiment, libraries, column='FPKM'):
    """Load all of the RSEM isoform quantifications files for an experiment

    Parameters
    ----------
    experiment: pandas.Series
        row of information about an experiment
    libraries: pandas.DataFrame
        Table of library metadata
    columns: str
       name of quantification column to load

    Returns
    -------
    pandas.DataFrame
       Table of all the quantification folumn for all isoform RSEM result files.
    """
    extension = '*_rsem.isoforms.results'
    return load_rsem_replicates(extension, experiment, libraries, column)


def load_rsem_replicates(extension, experiment, libraries, column):
    analysis_files = []
    library_ids = []

    assert isinstance(experiment, pandas.Series)
    for library_id in experiment.replicates:
        library_ids.append(library_id)
        library = libraries.loc[library_id]
        analysis_dir = library.analysis_dir
        pattern = os.path.join(analysis_dir, extension)
        result = glob(pattern)
        if len(result) != 1:
            raise RuntimeError(
                "Unexpected number of rsem gene result files {} using '{}'".format(
                    len(result),
                    pattern,
                )
            )
        elif not os.path.exists(result[0]):
            raise RuntimeError("{} doesn't exist".format(result[0]))
        else:
            analysis_files.append(result[0])
    quantifications = load_rsem_quantifications(
        analysis_files, index=library_ids, column=column
    )
    return quantifications


def create_quantification_cache(
        experiment, libraries, quantification_name,
        sep='\t'):
    score_filename = models.make_correlation_filename(experiment)
    quant_filename = models.make_quantification_filename(experiment,
                                                         quantification_name)

    quantifications = load_genomic_quantifications(
        experiment,
        libraries,
        quantification_name)
    if os.path.exists(quant_filename):
        os.unlink(quant_filename)

    store = pandas.HDFStore(quant_filename, complevel=9, complib='blosc')
    store.append('quantifications', quantifications)
    store.close()

    scores = compute_all_vs_all_scores(quantifications)
    if os.path.exists(score_filename):
        os.unlink(score_filename)

    store = pandas.HDFStore(score_filename)
    for key in scores:
        store.append(key, scores[key])
    store.close()


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    sep = get_seperator(args.sep)
    if args.experiments:
        experiments = models.load_experiments(args.experiments, sep=sep,
                                              analysis_root=args.root)
    else:
        if args.experiment_name is None:
            parser.error(
                "Please provide an experiment name. (Used as filename)")
        if len(args.replicates) == 0:
            parser.error(
                "Please provide list of replicates or experiment table")
        experiments = {args.experiment_name: args.replicates}

    if args.libraries is None:
        parser.error("Please provide library information tables")

    libraries = models.load_library_tables(args.libraries, sep=sep)

    for i, experiment in experiments.iterrows():
        logging.info('Processing:', experiment.name)
        create_quantification_cache(
            experiment,
            libraries,
            args.quantification,
            sep)


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--quantification', choices=['FPKM', 'TPM'],
                        default='FPKM',
                        help='which quantification value to use')
    # parser.add_argument('--model', choices=['gene', 'transcriptome'],
    #                     default='gene',
    #                     help="Quantify on genes or transcripts")

    parser.add_argument('-l', '--libraries', action='append',
                        help='library information table')
    parser.add_argument('-e', '--experiments', action='append',
                        help='experiments tables')
    parser.add_argument('-n', '--experiment-name', help='experiment name')
    parser.add_argument('replicates', nargs='*',
                        help='list of replicates for this experiment')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    parser.add_argument('--root', default=None,
                        help='write analysis files to this directory')
    return parser


if __name__ == '__main__':
    main()
