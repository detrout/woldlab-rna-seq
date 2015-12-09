#!/usr/bin/python3
"""Master DAG generator.

This is the script for automating generating the mapping, quantification,
and some quality control metrics.
"""
import argparse
import collections
from glob import glob
import os
import pandas

import make_star_rsem_dag
import models

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    sep = models.get_seperator(args.sep)
    libraries = models.load_library_tables(args.libraries, sep)
    fastqs = dict(find_fastqs(libraries))

    dag = generate_star_rsem_analysis(args, libraries, fastqs)
    print(dag)
    
    return 0

def make_parser():
    defaults = make_star_rsem_dag.read_defaults()
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--sep', choices=['TAB',','], default='TAB')
    parser.add_argument('libraries', nargs='+')
    
    parser.add_argument('--condor-script-dir',
                        help="specify the directory that the condor scripts located in",
                        default=defaults['condor_script_dir'])
    parser.add_argument('--genome-dir',
                        help="specify the directory that has the genome indexes",
                        default=defaults['genome_dir'])
    return parser

def find_fastqs(table, fastq_column='read_1'):
    """Find fastqs for a library from a library table

    fastqs are a comma seperated glob pattern
    """
    if fastq_column in table.columns:
        for library_id in table.index:
            fastqs = find_fastqs_by_glob(
                table.loc[library_id, fastq_column].split(','))
            yield (library_id, list(fastqs))
    else:
        # eventually look up by library ID
        raise NotImplemented("Please specify fastq glob")
        
def find_fastqs_by_glob(fastq_globs):
    for fastq in fastq_globs:
        for filename in glob(fastq):
            if os.path.exists(filename):
                yield os.path.abspath(filename)
            else:
                logger.warn("Can't find fastq {}. skipping".format(filename))


def generate_star_rsem_analysis(args, libraries, read_1_fastqs):
    dag = []
    for library_id in libraries.index:
        read_1_files = read_1_fastqs[library_id]

        analysis = make_star_rsem_dag.AnalysisDAG()

        analysis.condor_script_dir = args.condor_script_dir
        analysis.genome_dir = args.genome_dir
    
        analysis.genome = libraries.loc[library_id, 'genome']
        analysis.annotation = libraries.loc[library_id, 'annotation']
        analysis.sex = libraries.loc[library_id, 'sex']
        analysis.job_id = library_id
        analysis.analysis_dir = libraries.loc[library_id, 'analysis_dir']
        analysis.fastqs = read_1_files

        if analysis.is_valid():
            dag.append(str(analysis))
        else:
            raise ValueError("Unable to generate dagman script for {}".format(library_id))

    return os.linesep.join(dag)

if __name__ == '__main__':
    main()
    
