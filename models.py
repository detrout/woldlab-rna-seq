"""Load necessary information to run our analysis from text tables.
"""

import collections
from glob import glob
import os
from functools import partial
import pandas

AnalysisFile = collections.namedtuple('AnalysisFile', ['library_id', 'filename'])

def prepend_path(path, file_name):
    return os.path.join(path, file_name)

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
         analysis_dir and fastq read_1 columns.
      sep - separator character to use for table defaults to TAB

    Returns: a pandas dataframe containing required information

    Raises:
      ValueError if a required column is missing
      ValueError if there are duplicated library ids.
    """
    tables = []
    for library_file in table_filenames:
        library_file = os.path.abspath(library_file)
        path, name = os.path.split(library_file)
        table = pandas.read_csv(library_file, sep=sep,
                                index_col='library_id',
                                dtype={'library_id':str,
                                       'analysis_dir': str})
        required_library_columns_present(table)
        table.index = [str(x) for x in table.index]
        table.index.name = 'library id'
        table['analysis_dir'] = table['analysis_dir'].apply(partial(prepend_path, path))
        table['analysis_name'] = table['analysis_dir'].apply(os.path.basename)
        tables.append(table)

    libraries = pandas.concat(tables)
    validate_library_ids(libraries)
    return libraries

def required_library_columns_present(table):
    """Verify that a table contains required columns
    """
    missing = []
    for key in ['genome', 'sex', 'annotation', 'analysis_dir', 'read_1']:
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

def load_experiments(experiment_filenames, sep='\t'):
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
        experiment_filename = os.path.abspath(experiment_filename)
        table = pandas.read_csv(experiment_filename, sep=sep)
        tables.append(table)

    experiments = {}
    all_tables = pandas.concat(tables)

    for i, row in all_tables.iterrows():
        replicates = row.replicates.split(',')
        experiments[row.experiment] = replicates
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
            coverage.append(float(value))
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
        if len(filenames) == 0:
            raise RuntimeError("No files found in {} for {}".format(
                analysis_dir, extension))
        elif len(filenames) > 1:
            raise RuntimeError("To many files {}".format(filenames))
        else:
            yield AnalysisFile(library_id, filenames[0])


def load_correlations(experiment):
    """Load correlation panel for an experiment
    """
    correlations = collections.OrderedDict()
    correlation_filename = make_correlation_filename(experiment)
    if not os.path.exists(correlation_filename):
        raise RuntimeError(
            "Unable to open expected score file {}".format(
                correlation_filename))
    store = pandas.HDFStore(correlation_filename)
    scores = {}
    for key in store.keys():
        key_name = normalize_hdf_key(key)
        scores[key_name] = store[key]
    store.close()
    return pandas.Panel(scores)

def load_quantifications(experiment, quantification_name='FPKM'):
    """Load quantifications for an experiment
    """
    quantification_filename = make_quantification_filename(
        experiment,
        quantification_name)
    store = pandas.HDFStore(quantification_filename)
    for key in store.keys():
        quantifications =  store[key]
    store.close()
    return quantifications


def normalize_hdf_key(key):
    return key.replace('/', '')

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


def make_quantification_filename(experiment, quantification='FPKM'):
    return experiment + '_' + quantification + '.h5'


def get_single_spike_cpc():
    """Load single-cell copies-per-cell
    """
    return get_bulk_spike_cpc() / 1000.0

def get_bulk_spike_cpc():
    """Return bulk copies per cell
    """
    cpc = [
        ('gSpikein_ERCC-00130',  3612000,),
        ('gSpikein_ERCC-00004',   903000,),
        ('gSpikein_ERCC-00136',   225750,),
        ('gSpikein_ERCC-00108',   112875,),
        ('gSpikein_ERCC-00116',    56438,),
        ('gSpikein_ERCC-00092',    28219,),
        ('gSpikein_ERCC-00095',    14109,),
        ('gSpikein_ERCC-00131',    14109,),
        ('gSpikein_ERCC-00062',     7055,),
        ('gSpikein_ERCC-00019',     3527,),
        ('gSpikein_ERCC-00144',     3527,),
        ('gSpikein_ERCC-00170',     1764,),
        ('gSpikein_ERCC-00154',      882,),
        ('gSpikein_ERCC-00085',      882,),
        ('gSpikein_ERCC-00028',      441,),
        ('gSpikein_ERCC-00033',      220,),
        ('gSpikein_ERCC-00134',      220,),
        ('gSpikein_ERCC-00147',      110,),
        ('gSpikein_ERCC-00097',       55,),
        ('gSpikein_ERCC-00156',       55,),
        ('gSpikein_ERCC-00123',       28,),
        ('gSpikein_ERCC-00017',       14,),
        ('gSpikein_ERCC-00083',        3,),
        ('gSpikein_ERCC-00096',  1806000,),
        ('gSpikein_ERCC-00171',   451500,),
        ('gSpikein_ERCC-00009',   112875,),
        ('gSpikein_ERCC-00042',    56438,),
        ('gSpikein_ERCC-00060',    28219,),
        ('gSpikein_ERCC-00035',    14109,),
        ('gSpikein_ERCC-00025',     7055,),
        ('gSpikein_ERCC-00051',     7055,),
        ('gSpikein_ERCC-00053',     3527,),
        ('gSpikein_ERCC-00148',     1764,),
        ('gSpikein_ERCC-00126',     1764,),
        ('gSpikein_ERCC-00034',      882,),
        ('gSpikein_ERCC-00150',      441,),
        ('gSpikein_ERCC-00067',      441,),
        ('gSpikein_ERCC-00031',      220,),
        ('gSpikein_ERCC-00109',      110,),
        ('gSpikein_ERCC-00073',      110,),
        ('gSpikein_ERCC-00158',       55,),
        ('gSpikein_ERCC-00104',       28,),
        ('gSpikein_ERCC-00142',       28,),
        ('gSpikein_ERCC-00138',       14,),
        ('gSpikein_ERCC-00117',        7,),
        ('gSpikein_ERCC-00075',        2,),
        ('gSpikein_ERCC-00074',  1806000,),
        ('gSpikein_ERCC-00113',   451500,),
        ('gSpikein_ERCC-00145',   112875,),
        ('gSpikein_ERCC-00111',    56438,),
        ('gSpikein_ERCC-00076',    28219,),
        ('gSpikein_ERCC-00044',    14109,),
        ('gSpikein_ERCC-00162',     7055,),
        ('gSpikein_ERCC-00071',     7055,),
        ('gSpikein_ERCC-00084',     3527,),
        ('gSpikein_ERCC-00099',     1764,),
        ('gSpikein_ERCC-00054',     1764,),
        ('gSpikein_ERCC-00157',      882,),
        ('gSpikein_ERCC-00143',      441,),
        ('gSpikein_ERCC-00039',      441,),
        ('gSpikein_ERCC-00058',      220,),
        ('gSpikein_ERCC-00120',      110,),
        ('gSpikein_ERCC-00040',      110,),
        ('gSpikein_ERCC-00164',       55,),
        ('gSpikein_ERCC-00024',       28,),
        ('gSpikein_ERCC-00016',       28,),
        ('gSpikein_ERCC-00012',       14,),
        ('gSpikein_ERCC-00098',        7,),
        ('gSpikein_ERCC-00057',        2,),
        ('gSpikein_ERCC-00002',  1806000,),
        ('gSpikein_ERCC-00046',   451500,),
        ('gSpikein_ERCC-00003',   112875,),
        ('gSpikein_ERCC-00043',    56438,),
        ('gSpikein_ERCC-00022',    28219,),
        ('gSpikein_ERCC-00112',    14109,),
        ('gSpikein_ERCC-00165',     7055,),
        ('gSpikein_ERCC-00079',     7055,),
        ('gSpikein_ERCC-00078',     3527,),
        ('gSpikein_ERCC-00163',     1764,),
        ('gSpikein_ERCC-00059',     1764,),
        ('gSpikein_ERCC-00160',      882,),
        ('gSpikein_ERCC-00014',      441,),
        ('gSpikein_ERCC-00077',      441,),
        ('gSpikein_ERCC-00069',      220,),
        ('gSpikein_ERCC-00137',      110,),
        ('gSpikein_ERCC-00013',      110,),
        ('gSpikein_ERCC-00168',       55,),
        ('gSpikein_ERCC-00041',       28,),
        ('gSpikein_ERCC-00081',       28,),
        ('gSpikein_ERCC-00086',       14,),
        ('gSpikein_ERCC-00061',        7,),
        ('gSpikein_ERCC-00048',        2,),
    ]
    s = pandas.Series(collections.OrderedDict(cpc), dtype=float)
    s.index.name = 'spike-in'
    s.name = 'cpc'
    return s
