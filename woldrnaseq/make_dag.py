#!/usr/bin/python3
"""Master DAG generator.

This is the script for automating generating the mapping, quantification,
and some quality control metrics.
"""
from __future__ import absolute_import

import argparse
import datetime
import os
import logging
import pandas
from pkg_resources import resource_filename
from jinja2 import Environment, PackageLoader

from woldrnaseq import make_star_rsem_dag
from woldrnaseq import models
from woldrnaseq import __version__

from woldrnaseq.common import (
    add_default_path_arguments,
    add_debug_arguments,
    add_metadata_arguments,
    add_version_argument,
    configure_logging,
    find_fastqs,
    get_seperator,
    validate_path_args,
    validate_library_file_existance,
    validate_experiment_file_existance,
)

logger = logging.getLogger(__name__)


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    if not validate_path_args(args):
        parser.error('Please set required parameters')

    if not (validate_library_file_existance(args) and
            validate_experiment_file_existance(args)):
        parser.error('Fix path to files')

    sep = get_seperator(args.sep)
    library_filenames = args.libraries
    library_filenames.extend(args.other_libraries)

    libraries = models.load_library_tables(library_filenames, sep)
    read1 = dict(find_fastqs(libraries, 'read_1'))
    if 'read_2' in libraries.columns:
        read2 = dict(find_fastqs(libraries, 'read_2'))
    else:
        read2 = {}

    dags = generate_star_rsem_analysis(args, libraries, read1, read2)
    generate_combined_analysis(args, dags)

    return 0


def make_parser():
    parser = argparse.ArgumentParser(
        description="Generates a condor dagman script from a provided "
                    "library metadata file."
    )
    parser.add_argument('-o', '--output', default='run.dagman',
                        help='Name to save master dagman too')
    add_metadata_arguments(parser)
    parser.add_argument(
        'other_libraries', nargs='*',
        help="(Deprecated) specify a list of library metadata file names."
    )
    parser.add_argument(
        '-s', '--sep', choices=['TAB', ','], default='TAB',
        help="Specify the field separator character in the library metadata file"
    )
    add_default_path_arguments(parser)
    add_version_argument(parser)
    add_debug_arguments(parser)
    #parser.add_argument('--template', help='override default dagman template')

    return parser


def get_reference_prefix(libraries, library_id):
    """Get optional reference_prefix from model.

    Defaults to 'chr' if not present
    """
    prefix_name = 'reference_prefix'
    prefix_default = 'chr'
    if prefix_name not in libraries.columns:
        return prefix_default

    prefix = libraries.loc[library_id, prefix_name]
    if pandas.isnull(prefix):
        return prefix_default

    return prefix


def generate_star_rsem_analysis(args, libraries, read_1_fastqs, read_2_fastqs):
    """Generate dagmans for each analysis

    Parameters
    ----------
    args: Parsed Arguments from make_parser
    libraries: pandas.DataFrame containing library metadata
    read_1_fastqs: list of read 1 fastqs
    read_2_fastqs: list of read 2 fastqs

    Returns
    -------
    List of filenames to created per-analysis dagmans
    """
    dag = []
    for library_id in libraries.index:
        logger.debug("Creating script for %s", library_id)
        analysis = make_star_rsem_dag.AnalysisDAG()

        analysis.genome_dir = args.genome_dir
        analysis.star_dir = args.star_dir
        analysis.rsem_dir = args.rsem_dir
        analysis.georgi_dir = args.georgi_dir
        analysis.ucsc_tools_dir = args.ucsc_tools_dir

        analysis.genome = libraries.loc[library_id, 'genome']
        analysis.annotation = libraries.loc[library_id, 'annotation']
        analysis.sex = libraries.loc[library_id, 'sex']
        analysis.job_id = library_id
        analysis.analysis_dir = libraries.loc[library_id, 'analysis_dir']
        analysis.analysis_name = libraries.loc[library_id, 'analysis_name']
        analysis.read_1_fastqs = read_1_fastqs[library_id]
        analysis.read_2_fastqs = read_2_fastqs.get(library_id, [])
        analysis.stranded = libraries.loc[library_id, 'stranded']

        analysis.reference_prefix = get_reference_prefix(libraries, library_id)
        #if args.template:
        #    analysis.dagman_template = args.template

        if analysis.is_valid():
            target = os.path.join(analysis.analysis_dir, library_id + '.dagman')
            with open(target, 'wt') as outstream:
                outstream.write(str(analysis))
            dag.append({'library_id': library_id,
                        'subdag': target})
        else:
            raise ValueError("Unable to generate dagman script for {}".format(library_id))

    return dag


def generate_combined_analysis(args, dags):
    with open(args.output, 'wt') as outstream:
        env = Environment(loader=PackageLoader('woldrnaseq', 'templates'))
        template = env.get_template('full-encode.dagman')

        library_tsv = ' '.join(['-l '+l for l in args.libraries])
        experiment_tsv = ' '.join(['-e '+e for e in args.experiments])

        outstream.write(template.render(
            madqc=resource_filename(__name__, 'madqc.condor'),
            makersemcsv=resource_filename(__name__, 'makersemcsv.condor'),
            report=resource_filename(__name__, 'report.condor'),
            libraries=library_tsv,
            experiments=experiment_tsv,
            genome_dir=args.genome_dir,
            dags=dags,
            username=os.getlogin(),
            timestamp=datetime.datetime.now().isoformat(),
            woldrnaseq_version=__version__,
        ))


if __name__ == '__main__':
    main()
