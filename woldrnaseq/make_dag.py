#!/usr/bin/python3
"""Master DAG generator.

This is the script for automating generating the mapping, quantification,
and some quality control metrics.
"""
from __future__ import absolute_import

import argparse
import os
import logging
import pandas

from woldrnaseq import make_star_rsem_dag
from woldrnaseq import models
from woldrnaseq.version import get_git_version

from woldrnaseq.common import (
    add_default_path_arguments,
    add_debug_arguments,
    add_version_argument,
    configure_logging,
    find_fastqs,
    get_seperator,
    validate_args
)

logger = logging.getLogger(__name__)

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    if args.version:
        parser.exit(0, 'version: %s\n' % (get_git_version(),))

    if not validate_args(args):
        parser.error("Please set required parameters")

    sep = get_seperator(args.sep)
    library_filenames = args.libraries
    library_filenames.extend(args.other_libraries)

    libraries = models.load_library_tables(library_filenames, sep)
    read1 = dict(find_fastqs(libraries, 'read_1'))
    if 'read_2' in libraries.columns:
        read2 = dict(find_fastqs(libraries, 'read_2'))
    else:
        read2 = {}

    dag = generate_star_rsem_analysis(args, libraries, read1, read2)
    print(dag)
    
    return 0

def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sep', choices=['TAB',','], default='TAB')
    parser.add_argument('-l', '--libraries', action='append', default=[])
    parser.add_argument('other_libraries', nargs='*')
    add_default_path_arguments(parser)
    add_version_argument(parser)
    add_debug_arguments(parser)
    parser.add_argument('--template', help='override default dagman template')
    
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

        analysis.reference_prefix = get_reference_prefix(libraries, library_id)
        analysis.dagman_template = args.template

        if analysis.is_valid():
            dag.append(str(analysis))
        else:
            raise ValueError("Unable to generate dagman script for {}".format(library_id))

    return os.linesep.join(dag)

if __name__ == '__main__':
    main()
    
