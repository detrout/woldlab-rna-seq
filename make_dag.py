#!/usr/bin/python3

import argparse
import collections
from glob import glob
import os
import pandas

import make_star_rsem_dag

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    sep = get_seperator(args.sep)
    libraries = load_library_tables(args.libraries, sep)
    fastqs = list(find_fastqs(libraries))

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

def get_seperator(sep):
    if sep.lower() == 'tab':
        return '\t'
    elif sep == ',':
        return ','
    else:
        raise ValueError("Unrecognized seperator")

def load_library_tables(table_filenames, sep):
    tables = []
    for library_file in table_filenames:
        table = pandas.read_csv(library_file, sep)
        required_columns_present(table)
        tables.append(table)

    libraries = pandas.concat(tables)
    validate_library_ids(libraries)
    return libraries
    
def required_columns_present(table):
    missing = []
    for key in ['library_id', 'genome', 'sex', 'annotation', 'analysis_dir', 'fastqs']:
        if key not in table.columns:
            missing.append(key)
    if len(missing) != 0:
        raise ValueError("Required columns missing: {}".format(','.join(mising)))

def validate_library_ids(table):
    library_ids = collections.Counter()
    for library_id in table['library_id']:
        library_ids[library_id] += 1

    duplicates = []
    for library_id in library_ids:
        if library_ids[library_id] > 1:
            duplicates.append(library_id)

    if len(duplicates) > 0:
        raise ValueError("Duplicate library ids: {}".format(duplicates))


def find_fastqs(table):
    """Find fastqs for a library from a library table

    fastqs are a comma seperated glob pattern
    """
    if 'fastqs' in table.columns:
        for i in table.index:
            fastqs = find_fastqs_by_glob(table.loc[i, 'fastqs'].split(','))
            yield (table.loc[i, 'library_id'], fastqs)
    else:
        # eventually look up by library ID
        raise NotImplemented("Please specify fastq glob")
        
def find_fastqs_by_glob(fastq_globs):
    for fastq in fastq_globs:
        for filename in glob(fastq):
            if os.path.exists(filename):
                yield filename
            else:
                logger.warn("Can't find fastq {}. skipping".format(filename))


def generate_star_rsem_analysis(args, libraries, fastqs):
    dag = []
    for i in libraries.index:
        library_id = libraries.loc[i, 'library_id']
        fastq_library_id, filenames = fastqs[i]
        assert library_id == fastq_library_id

        analysis = make_star_rsem_dag.AnalysisDAG()

        analysis.condor_script_dir = args.condor_script_dir
        analysis.genome_dir = args.genome_dir
    
        analysis.genome = libraries.loc[i, 'genome']
        analysis.annotation = libraries.loc[i, 'annotation']
        analysis.sex = libraries.loc[i, 'sex']
        analysis.job_id = libraries.loc[i, 'library_id']
        analysis.analysis_dir = libraries.loc[i, 'analysis_dir']
        analysis.fastqs = filenames

        if analysis.is_valid():
            dag.append(str(analysis))
        else:
            raise ValueError("Unable to generate dagman script for {}".format(library_id))

    return os.linesep.join(dag)

if __name__ == '__main__':
    main()
    
