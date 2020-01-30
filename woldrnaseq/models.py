"""Load necessary information to run our analysis from text tables.
"""
import collections
from glob import glob
import logging
import os
from functools import partial
import pandas

from woldrnaseq.common import validate_reference_type

logger = logging.getLogger(__name__)


def read_line_from_stream(stream):
    """Read a line filtering out blank lines and comments
    """
    for line in stream:
        line = line.strip()
        if len(line) > 0 and not line.startswith('#'):
            yield line


def load_experiments(experiment_filenames, sep='\t', analysis_root=None):
    """Load metadata table describing a set of experiments

    Parameters
    ----------
      experiment_filenames: list of str
         list of filenames to load the tables need to have a
         experiment name and list of library ids that are intended to
         be treated as related replicates
      sep: str
         separator character to use for table defaults to TAB

    Returns
    -------
    pandas.DataFrame
      A table of experiment names and a list of replicates
    """
    tables = []
    for experiment_filename in experiment_filenames:
        # Maybe I need a better test for handling remote urls
        if os.path.exists(experiment_filename):
            warn_if_spaces(experiment_filename)
            experiment_filename = os.path.abspath(experiment_filename)
        if analysis_root is None:
            analysis_cur_root, name = os.path.split(experiment_filename)
        else:
            analysis_cur_root = analysis_root

        table = pandas.read_csv(experiment_filename, sep=sep, comment='#',
                                skip_blank_lines=True,
                                header=0,
                                dtype={
                                    'experiment': str,
                                },
                                converters={
                                    'replicates': lambda x: x.split(',')
                                }
        )
        table.set_index('experiment', inplace=True)
        required_experiment_columns_present(table)
        table['analysis_dir'] = analysis_cur_root
        tables.append(table)

    return pandas.concat(tables)


def load_library_tables(table_filenames, sep='\t', analysis_root=None):
    """Load table describing libraries to be analyized

    Parameters
    ----------
    table_filenames: list of str
        a list of filenames to load the tables. The tables are
        required to have library_id, genome, sex, annotation,
        analysis_dir and read_1 column headings.

    sep: str
        separator character to use for table defaults to TAB
    analysis_root: str
         analysis_dirs are relative to this path, defaults to location
         of library.txt file

    Returns
    -------
    pandas.DataFrame
         Containing the columns from the input library and a inferred columns
         analysis_dir and analysis_name

    Raises
    ------
    ValueError
        if a required column is missing
    ValueError
        if there are duplicated library ids.
    """
    assert not isinstance(table_filenames, str)
    tables = []
    for library_file in table_filenames:
        # Maybe I need a better test for handling remote urls
        if os.path.exists(library_file):
            warn_if_spaces(library_file)
            library_file = os.path.abspath(library_file)
        if analysis_root is None:
            analysis_cur_root, name = os.path.split(library_file)
        else:
            analysis_cur_root = analysis_root

        table = pandas.read_csv(library_file, sep=sep,
                                index_col='library_id',
                                dtype={'library_id': str,
                                       'analysis_dir': str},
                                converters={
                                    'stranded': _normalize_stranded
                                },
                                comment='#',
                                skip_blank_lines=True,
                                )
        verify_library_columns(table)
        table.index = [str(x) for x in table.index]
        table.index.name = 'library id'
        table['analysis_dir'] = table['analysis_dir'].apply(partial(os.path.join, analysis_cur_root))
        table['analysis_name'] = table['analysis_dir'].apply(os.path.basename)
        if 'stranded' not in table.columns:
            table['stranded'] = 'unstranded'
        tables.append(table)

    libraries = pandas.concat(tables)
    validate_library_ids(libraries)
    return libraries


def _normalize_stranded(value):
    """Return standardized strand names from library control file

    This is intended to be used internally by load_library_tables

    Parameters
    ----------
    value: str
      string containing forward/+, reverse/-, or unstranded/blank

    Returns
    -------
      forward, reverse, or unstranded
    """
    if pandas.isnull(value) or value.lower() in ('unstranded'):
        return 'unstranded'
    elif value.lower() in ('forward', '+'):
        return 'forward'
    elif value.lower() in ('reverse', '-'):
        return 'reverse'
    else:
        raise ValueError("Unrecognized strand {}".format(value))


def genome_name_from_library(row):
    """Generate genome name triple from a library row

    Parameters
    ----------
    row: pandas.Series
        row of a library metadata table.

    Returns
    -------
    str
        Combined genome, annotation, and sex string.
    """
    return '-'.join([row.genome, row.annotation, row.sex])


def verify_library_columns(table):
    """Verify that a table contains required columns

    Parameters
    ----------
    table : pandas.DataFrame
        Verify that a library metadata table contains the required
        fields. Additionally we also check and warn if there are any
        unrecognized optional columns present.

    Raises
    ------
    ValueError : if a required column is missing
    """
    missing = []
    required = {'genome', 'sex', 'annotation', 'analysis_dir', 'read_1'}
    optional = {'read_2', 'reference_prefix', 'stranded'}
    known = required.union(optional)

    columns = set(table.columns)
    missing = required.difference(columns)
    if len(missing) != 0:
        raise ValueError("Required columns missing: {}".format(','.join(missing)))

    unknown = columns.difference(known)
    if len(unknown) > 0:
        logger.warning('Unrecognized columns present. Is this intended?: {}'.format(
            ','.join(unknown)))


def validate_library_ids(table):
    """Validate that library ids are unique

    Parameters
    ----------
    table : pandas.DataFrame
        A library metadata table

    Raises
    ------
    ValueError : if there are duplicate library_ids present
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


def required_experiment_columns_present(table):
    """Verify that a experiment table contains required columns
    """
    missing = set(('replicates',)).difference(table.columns)
    if table.index.name != 'experiment':
        missing.add('experiment')

    if len(missing) != 0:
        raise ValueError("Required columns missing: {}".format(','.join(missing)))


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


def parse_star_final_fieldname(name):
    name = name.strip()
    for sep in ['|', ':']:
        if name.endswith(sep):
            name = name[:-1]
    return name.strip()


def parse_star_final_percent(x):
    return float(x.replace('%', ''))


def load_star_counts(filename, column):
    """Load star count table

    Parameters
    ----------
    filename: str
        Path to a STAR ReadsPerGene.out.tab
    column: str
        Which column to return.
        U for unstranded
        + for first read
        - for second read

    Returns
    -------
    pandas.DataFrame
        Containing the annotation id and the requested column.
    """
    count_columns = {
        'U': 1,
        '+': 2,
        '-': 3,
    }

    if column not in count_columns:
        raise ValueError('Expected quantification column identifier {}'.format(list(count_columns)))

    data = pandas.read_csv(filename,
                           skiprows=4,
                           header=None,
                           index_col=0,
                           usecols=[0, count_columns[column]],
                           sep='\t')
    data.index.name = 'gene_id'
    data.columns = [column]
    return data


def load_all_star_counts(libraries, column):
    """Load STAR gene count tables

    Parameters
    ----------
    libraries: pandas.DataFrame
        Table of library metadata to load
    column: str
        Which column to return.
        U for unstranded
        + for first read
        - for second read
    """
    library_ids = []
    counts = []
    analysis_files = find_library_analysis_file(libraries, 'ReadsPerGene.out.tab')
    for library_id, filename in analysis_files:
        counts.append(load_star_counts(filename, column))
        library_ids.append(library_id)

    expression = pandas.concat(counts, axis=1)
    expression.columns = library_ids
    return expression


def load_star_final_log(filename):
    name_type = {
        'Started job on': str,
        'Started mapping on': str,
        'Finished on': str,
        'Uniquely mapped reads %': parse_star_final_percent,
        'Mismatch rate per base, %': parse_star_final_percent,
        'Deletion rate per base': parse_star_final_percent,
        'Insertion rate per base': parse_star_final_percent,
        '% of reads mapped to multiple loci': parse_star_final_percent,
        '% of reads mapped to too many loci': parse_star_final_percent,
        '% of reads unmapped: too many mismatches': parse_star_final_percent,
        '% of reads unmapped: too short': parse_star_final_percent,
        '% of reads unmapped: other': parse_star_final_percent,
        '% of chimeric reads': parse_star_final_percent,
    }
    prefix = ''
    index = []
    values = []
    with open(filename) as instream:
        for line in instream:
            fields = line.split('\t')
            if len(fields) == 0:
                # blank line
                continue
            elif len(fields) == 1:
                # header
                prefix = parse_star_final_fieldname(fields[0])
            else:
                name = parse_star_final_fieldname(fields[0])
                index.append((prefix, name))
                values.append(name_type.get(name, float)(fields[1].strip()))

    index = pandas.MultiIndex.from_tuples(index, names=['read_class', 'name'])
    return pandas.Series(values, index)


def load_all_star_final(libraries):
    final = []
    library_ids = []
    analysis_files = find_library_analysis_file(libraries, 'Log.final.out')
    for library_id, filename in analysis_files:
        data = load_star_final_log(filename)
        data.name = library_id
        final.append(data)
        library_ids.append(library_id)

    return pandas.DataFrame(final)


def load_distribution(filename):
    distribution = collections.OrderedDict()
    with open(filename, 'rt') as instream:
        try:
            for line in read_line_from_stream(instream):
                name, value = line.split(':\t')
                distribution[name] = float(value)
        except ValueError as e:
            raise ValueError('Unable to parse {}: {}'.format(filename, str(e)))
    return pandas.Series(distribution, index=distribution.keys())


def load_all_distribution(libraries):
    distribution = []
    library_ids = []
    analysis_files = find_library_analysis_file(libraries, '*.sam_reads_genes')
    for library_id, filename in analysis_files:
        distribution.append(load_distribution(filename))
        library_ids.append(library_id)
    return pandas.DataFrame(distribution, index=library_ids)


def load_coverage(filename, library_id=None):
    logger.info('loading {}'.format(filename))
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


def load_gene_coverage(filename, library_id, gene_normalization):
    """Load a geneList per gene coverage file
    """
    coverage = pandas.read_csv(filename, sep='\t', index_col=0, header=None)
    if gene_normalization == 'max':
        coverage = coverage.divide(coverage.max(axis=1), axis=0)
    coverage.name = library_id
    return coverage


def load_all_gene_coverage(libraries, gene_list, gene_normalization):
    if len(libraries) > 0:
        analysis_files = find_library_analysis_file(libraries, '*.geneList')
        for library_id, filename in analysis_files:
            yield load_gene_coverage(filename, library_id, gene_normalization)
    for filename in gene_list:
        library_id = os.path.split(filename.replace('.coverage.geneList', ''))[1]
        yield load_gene_coverage(filename, library_id, gene_normalization)


AnalysisFile = collections.namedtuple('AnalysisFile', ['library_id', 'filename'])
AnalysisFile.__doc__ = "Tuple of library_ids and the full path to a particular generated file type."


def find_library_analysis_file(libraries, extension):
    for library_id in libraries.index:
        analysis_dir = libraries.loc[library_id, 'analysis_dir']
        assert analysis_dir is not None
        filenames = glob(os.path.join(analysis_dir, extension))
        if len(filenames) == 0:
            logger.warning("No files found in {} for {}".format(
                analysis_dir, extension))
        elif len(filenames) > 1:
            raise RuntimeError("To many files {}".format(filenames))
        else:
            yield AnalysisFile(library_id, filenames[0])


def find_library_bam_file(library, reference_type='genome', analysis_root=None):
    """Generate the path to where the bam for this library is

    :param Series library: row from a library table DataFrame
    :param str reference_type: bam reference type, genome, transcriptome
    :param str analysis_root: root directory to be searching for track files
    :returns: path of bam file relative to analysis_root
    """
    validate_reference_type(reference_type)

    if 'genome' == reference_type:
        extension = '_genome.bam'
    else:
        extension = '_anno.bam'

    genome_triplet = genome_name_from_library(library)
    bam_name = library.analysis_name + '-' + genome_triplet + extension
    to_check = [
        os.path.join(library.analysis_dir, bam_name),
    ]
    if analysis_root is not None:
        to_check.append(os.path.join(analysis_root, bam_name))

    if 'genome' == reference_type:
        to_check.append(os.path.join(library.analysis_dir, 'Aligned.sortedByCoord.out.bam'))

    for pathname in to_check:
        if os.path.exists(pathname):
            if 'genome' == reference_type:
                bai = pathname + '.bai'
                if not os.path.exists(bai):
                    logger.warning('Missing index file for {}'.format(pathname))
            return pathname


def load_correlations(experiment):
    """Load correlation panel for an experiment
    """
    correlation_filename = make_correlation_filename(experiment)

    if not os.path.exists(correlation_filename):
        raise FileNotFoundError(
            "Unable to open expected score file {}".format(
                correlation_filename))
    store = pandas.HDFStore(correlation_filename)
    scores = {}
    for key in store.keys():
        key_name = normalize_hdf_key(key)
        scores[key_name] = store[key]
        logger.debug('Loading %s shape %s', key_name, scores[key_name].shape)
    store.close()
    return scores


def load_quantifications(experiment, quantification_name='FPKM'):
    """Load quantifications for an experiment
    """
    assert isinstance(experiment, pandas.Series)

    quantification_filename = make_quantification_filename(
        experiment,
        quantification_name)

    logger.debug("Opening quantification file: %s", quantification_filename)
    quantifications = None
    if os.path.exists(quantification_filename):
        store = pandas.HDFStore(quantification_filename)
        if len(store.keys()) == 0:
            logger.error("Quantification cache file %s is empty",
                         quantification_filename)
        else:
            for key in store.keys():
                quantifications = store[key]
        store.close()
        return quantifications
    else:
        logger.error("Quantification cache file %s not available",
                     quantification_filename)

        return None


def load_gtf_cache(filename):
    """Load GTF cache file produced by gff2table

    :Paramters:
      - filename: (str) path to .h5 gtf cache file
    :Returns: DataFrame of annotations
    """
    if os.path.exists(filename):
        store = pandas.HDFStore(filename, 'r')
        assert len(store.keys()) == 1
        annotation = store[store.keys()[0]]
        store.close()
        return annotation
    else:
        raise FileNotFoundError('Unable to load file {}'.format(filename))


def lookup_gene_name_by_gene_id(annotation, table):
    """Add gene names using gene_id
    """
    return lookup_gene_name_by_id(annotation, table, 'gene_id')


def lookup_gene_name_by_transcript_id(annotation, table):
    """Add gene names using transcript_id
    """
    return lookup_gene_name_by_id(annotation, table, 'transcript_id')


def lookup_gene_name_by_id(annotation, table, column):
    """Add gene names using specified column
    """
    column_id_to_type = {
        'gene_id': 'gene',
        'transcript_id': 'transcript',
    }
    annotation_type = column_id_to_type.get(column)
    if annotation_type is None:
        raise ValueError("Unsupported column type %s" % (column,))

    annotation = annotation.set_index(column)
    annotation_filter = annotation['type'] == annotation_type
    annotated_table = table.merge(
        annotation[['gene_name']][annotation_filter],
        left_index=True,
        right_index=True,
        how='left')
    # we should have the same number of rows as the initial table
    assert annotated_table.shape[0] == table.shape[0]
    preferred_order = ['gene_name'] + list(table.columns)
    return annotated_table[preferred_order].fillna('')


def normalize_hdf_key(key):
    return key.replace('/', '')


def make_correlation_filename(experiment, reference_type='genome'):
    validate_reference_type(reference_type)
    assert isinstance(experiment, pandas.Series)
    components = [experiment.name]
    if reference_type == 'transcriptome':
        components.append(reference_type)
    components.append('correlation')
    name = '_'.join(components) + '.h5'
    return os.path.join(experiment['analysis_dir'], name)


def make_quantification_filename(
        experiment,
        quantification='FPKM',
        reference_type='gene'):
    validate_reference_type
    assert isinstance(experiment, pandas.Series)
    components = [experiment.name]
    if reference_type == 'transcriptome':
        components.append(reference_type)

    components.append(quantification)

    name = '_'.join(components) + '.h5'
    return os.path.join(experiment['analysis_dir'], name)


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


def warn_if_spaces(filename):
    with open(filename, 'rt') as instream:
        line = instream.readline()
        if ' ' in line:
            logger.warning("There are spaces in the header line, is this intentional?")
