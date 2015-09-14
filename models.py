"""Load necessary information to run our analysis from text tables.
"""

import collections
from glob import glob
import os
import pandas

AnalysisFile = collections.namedtuple('AnalysisFile', ['library_id', 'filename'])

def read_line_from_stream(stream):
    """Read a line filtering out blank lines and comments.
    """
    for line in stream:
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            yield line


def load_library_tables(table_filenames, sep='\t'):
    """Load table describing libraries to be analyized

    Parameters:
      table_filenames - list of filenames to load
         the tables need to have a library_id, genome, sex, annotation,
         analysis_dir and fastqs columns.
      sep - separator character to use for table defaults to TAB

    Returns: a pandas dataframe containing required information

    Raises:
      ValueError if a required column is missing
      ValueError if there are duplicated library ids.
    """
    tables = []
    for library_file in table_filenames:
        table = pandas.read_csv(library_file, sep=sep,
                                index_col='library_id',
                                dtype={'library_id':str})
        required_library_columns_present(table)
        table.index = [str(x) for x in table.index]
        table.index.name = 'library id'
        tables.append(table)

    libraries = pandas.concat(tables)
    validate_library_ids(libraries)
    return libraries

def required_library_columns_present(table):
    """Verify that a table contains required columns
    """
    missing = []
    for key in ['genome', 'sex', 'annotation', 'analysis_dir', 'fastqs']:
        if key not in table.columns:
            missing.append(key)
    if len(missing) != 0:
        raise ValueError("Required columns missing: {}".format(','.join(missing)))

def validate_library_ids(table):
    """Validate that library ids are unique
    """
    library_ids = collections.Counter()
    for library_id in table.index:
        library_ids[library_id] += 1

    duplicates = []
    for library_id in library_ids:
        if library_ids[library_id] > 1:
            duplicates.append(library_id)

    if len(duplicates) > 0:
        raise ValueError("Duplicate library ids: {}".format(duplicates))

def load_experiments(experiment_filenames, libraries, sep='\t'):
    """Load table describing experiments

    Parameters:
      experiment_filenames - list of filenames to load
         the tables need to have a experiment name and list of
         library ids that are intended to be treated as related
         replicates
      sep - separator character to use for table defaults to TAB

    Returns a dictionary mapping experiment names to list of replicates
    """
    tables = []
    for experiment_filename in experiment_filenames:
        table = pandas.read_csv(experiment_filename, sep)
        tables.append(table)

    experiments = {}
    all_tables = pandas.concat(tables)
    for i in all_tables.index:
        experiment = all_tables.loc[i, 'experiment']
        replicates = all_tables.loc[i, 'replicates'].split(',')
        experiments[experiment] = replicates
    return experiments


def load_samstats(filename):
    floats = set(['Fraction Mapped', 'Complexity', 'Read Length, Average'])
    samstats = collections.OrderedDict()

    with open(filename, 'r') as instream:
        for line in read_line_from_stream(instream):
            name, value = line.split(':\t')
            if name in floats:
                samstats[name] = float(value)
            else:
                samstats[name] = int(value)

    return pandas.Series(samstats, index=samstats.keys())


def load_all_samstats(libraries):
    samstats = []
    library_ids = []
    analysis_files = find_library_analysis_file(libraries, '*.samstats')
    for library_id, filename in analysis_files:
        samstats.append(load_samstats(filename))
        library_ids.append(library_id)

    return pandas.DataFrame(samstats, index=library_ids)


def load_distribution(filename):
    distribution = collections.OrderedDict()
    with open(filename, 'rt') as instream:
        for line in read_line_from_stream(instream):
            name, value = line.split(':\t')
            distribution[name] = float(value)
    return pandas.Series(distribution, index=distribution.keys())


def load_all_distribution(libraries):
    distribution = []
    library_ids = []
    analysis_files = find_library_analysis_file(libraries, '*.sam_reads_genes')
    for library_id, filename in analysis_files:
        distribution.append(load_distribution(filename))
        library_ids.append(library_id)
    return pandas.DataFrame(distribution, index=library_ids)
    #return libraries.merge(distribution_df, left_index=True, right_index=True)


def load_coverage(filename, library_id=None):
    with open(filename, 'rt') as instream:
        coverage = []
        for line in instream:
            position, value = line.strip().split('\t')
            coverage.append(value)
    return pandas.DataFrame(coverage, columns=[library_id])


def load_all_coverage(libraries):
    coverage = []
    analysis_files = find_library_analysis_file(libraries, '*.coverage')
    for library_id, filename in analysis_files:
        coverage.append(load_coverage(filename, library_id))
    return pandas.concat(coverage, axis=1)


def find_library_analysis_file(libraries, extension):
    for library_id in libraries.index:
        analysis_dir = libraries.loc[library_id, 'analysis_dir']
        filenames = glob(os.path.join(analysis_dir, extension))
        if len(filenames) != 1:
            raise ("To many files {}".format(filenames))
        else:
            yield AnalysisFile(library_id, filenames[0])


def load_all_correlations(experiments):
    correlations = collections.OrderedDict()
    for name in experiments:
        replicates = experiments[name]
        correlation_filename = make_correlation_filename(name)
        if not os.path.exists(correlation_filename):
            raise RuntimeError(
                "Unable to open expected score file {}".format(
                    correlation_filename))
        store = pandas.HDFStore(correlation_filename)
        scores = {}
        for key in store.keys():
            key_name = key.replace('/','')
            scores[key_name] = store[key]
        store.close()
        correlations[name] = pandas.Panel(scores)
    return correlations


# really doesn't belong here
def get_seperator(sep):
    if sep.lower() == 'tab':
        return '\t'
    elif sep == ',':
        return ','
    else:
        raise ValueError("Unrecognized seperator")


def make_correlation_filename(experiment):
    return experiment + '_correlation.h5'
