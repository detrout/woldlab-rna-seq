"""Common utilities shared by several components
"""
import configparser
from glob import glob
import logging
import os
from pathlib import Path
import subprocess
from collections import abc
import numpy

from . import __version__

logger = logging.getLogger(__name__)


def add_metadata_arguments(parser):
    parser.add_argument('-l', '--libraries', action='append',
                        help='library information table')
    parser.add_argument('-e', '--experiments', action='append',
                        help='experiments tables')


def add_default_path_arguments(parser):
    """Add arguments to allow overriding location of dependencies
    """
    defaults = read_defaults()
    group = parser.add_argument_group(
        title='Paths',
        description='Paths to required external programs. Defaults is set from .htsworkflow.ini')
    group.add_argument(
        '--genome-dir',
        help="specify the directory that has the genome indexes",
        default=defaults['genome_dir'])
    group.add_argument(
        '--star-dir',
        default=defaults['star_dir'],
        help='Specify the directory where STAR is installed')
    group.add_argument(
        '--rsem-dir',
        default=defaults['rsem_dir'],
        help='Specify the directory where rsem is installed')
    group.add_argument(
        '--georgi-dir',
        default=defaults['georgi_dir'],
        help='Specify the directory where georgi scripts are installed')
    group.add_argument(
        '--ucsc-tools-dir',
        default=defaults['ucsc_tools_dir'],
        help='Specify the directory where the UCSC tools are installed')
    return parser


def add_debug_arguments(parser):
    """Add arguments for tuning logging
    """
    group = parser.add_argument_group('Verbosity', 'Set logging level')
    group.add_argument(
        '-v', '--verbose', default=False, action='store_true',
        help="Display informational level log messages"
    )
    group.add_argument(
        '-d', '--debug', default=False, action='store_true',
        help="Display debugging level log messages"
    )
    return group


def add_separator_argument(parser):
    """Add argument for setting table separator comma or tab
    """
    parser.add_argument(
        '-s', '--sep', choices=['TAB', ','], default='TAB',
        help="Specify the field separator character in the library metadata file"
    )
    return parser


def add_version_argument(parser):
    """Add Version argument
    """
    parser.add_argument(
        '--version', action='version', version=__version__
    )
    return parser


def configure_logging(args, **kwargs):
    """run logging.basicConfig based on common verbosity command line arguments
    """
    if args.debug:
        logging.basicConfig(level=logging.DEBUG, **kwargs)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO, **kwargs)
    else:
        logging.basicConfig(level=logging.WARNING, **kwargs)


def get_seperator(sep):
    """Parse seperator argument
    """
    if sep.lower() == 'tab':
        return '\t'
    elif sep == ',':
        return ','
    else:
        raise ValueError("Unrecognized seperator")


def read_defaults():
    """Read defaults from configuration file
    """
    defaults = {
        'genome_dir': None,
        'star_dir': None,
        'rsem_dir': None,
        'georgi_dir': None,
        'ucsc_tools_dir': None,
    }
    config = configparser.ConfigParser()
    config.read([
        Path('~/.htsworkflow.ini').expanduser(),
        Path('/etc/htsworkflow.ini')
    ])

    if config.has_section('analysis'):
        analysis = config['analysis']
        for name in defaults.keys():
            defaults[name] = normalize_path(analysis.get(name))

    return defaults


def normalize_path(path):
    """Make sure defined paths end with a /

    Don't do anything to empty or None paths.
    """
    if path is None:
        return None
    elif len(path) == 0:
        return path
    else:
        return os.path.join(os.path.expanduser(path), '')


def get_fixed_range(limits):
    """Return minimum and maxmium limits
    """
    # this wont work if we have the larger value on the bottom
    assert numpy.all([a[0] <= a[1] for a in limits]), 'Constraint error disordered list {}'.format(limits)
    a = min((a[0] for a in limits))
    b = max((a[1] for a in limits))
    return (a, b)


def save_fixed_height(plots):
    """reformat plot heights based on largest height
    """

    limits = []
    for name in plots:
        for a in plots[name].axes:
            limits.append(a.get_ylim())

    fixed = get_fixed_range(limits)

    for png_name in plots:
        f = plots[png_name]
        for a in f.axes:
            a.set_ylim(fixed)
        print('saving {}'.format(png_name))
        f.savefig(png_name, bbox_inches='tight')


def validate_path_args(args):
    """Warn if path arguments aren't set

    Parameters
    ----------
    args : argparse.Namespace
        parsed arguments as generated by argparse.parse_args

    Returns
    -------
    bool
        True if we have required arguments
    """
    can_continue = True
    if args.genome_dir is None:
        logger.error("Need path to genome indexes")
        can_continue = False

    if args.star_dir is None:
        logger.warning("Path to STAR not provided, assuming its on the PATH")

    if args.rsem_dir is None:
        logger.warning("Path to rsem-calculate-expression not provided, assuming its on the PATH")

    if args.georgi_dir is None:
        logger.error('Path to "GeorgiScripts" python scripts not provided.')
        can_continue = False

    return can_continue


def validate_experiment_file_existance(args):
    """Warn if experiment metadata files are missing

    Parameters
    ----------
    args : argparse.Namespace
        parsed arguments as generated by argparse.parse_args

    Returns
    -------
    bool
        True if we have required arguments
    """
    can_continue = True

    for exp in args.experiments:
        if not os.path.exists(exp):
            logger.error('{} does not exist'.format(exp))
            can_continue = False

    return can_continue


def validate_library_file_existance(args):
    """Warn if library metadata files are missing

    Parameters
    ----------
    args : argparse.Namespace
        parsed arguments as generated by argparse.parse_args

    Returns
    -------
    bool
        True if we have required arguments
    """
    can_continue = True

    for lib in args.libraries:
        if not os.path.exists(lib):
            logger.error('{} does not exist'.format(lib))
            can_continue = False

    for exp in args.experiments:
        if not os.path.exists(exp):
            logger.error('{} does not exist'.format(exp))
            can_continue = False

    return can_continue


def validate_reference_type(reference_type):
    if reference_type not in ['genome', 'transcriptome']:
        raise ValueError(
            'Reference type must be genome or transcriptome not {}'.format(
                reference_type))


def find_fastqs(table, fastq_column):
    """Find fastqs for a library from a library table

    fastqs are a comma seperated glob pattern
    """
    assert fastq_column in ['read_1', 'read_2'], "Unrecognized fastq column name {}".format(fastq_column)
    if fastq_column in table.columns:
        for library_id in table.index:
            fastq_field = table.loc[library_id, fastq_column]
            if isinstance(fastq_field, abc.Sequence):
                fastqs = find_fastqs_by_glob(fastq_field)
                yield (library_id, list(fastqs))
    else:
        # eventually look up by library ID
        raise NotImplementedError("Please specify fastq glob")


def find_fastqs_for_library(table, library_id, fastq_column):
    """Find fastqs for a library from a library table

    fastqs are a comma seperated glob pattern
    """
    assert fastq_column in ['read_1', 'read_2'], "Unrecognized fastq column name {}".format(fastq_column)
    if fastq_column in table.columns and library_id in table.index:
        fastq_field = table.loc[library_id, fastq_column]
        if isinstance(fastq_field, abc.Sequence):
            fastqs = find_fastqs_by_glob(fastq_field)
            yield (library_id, list(fastqs))
    else:
        # eventually look up by library ID
        raise NotImplementedError("Please specify fastq glob {} {}".format(fastq_column in table.columns, library_id in table.index))


def find_fastqs_by_glob(fastq_globs):
    """Generate a file names from the provided list of names

    Parameters
    ----------
    fastq_globs: list of str
        A list of filenames that potentally contain the file name
        pattern match characters '*' and '?'

    Returns
    -------
    generator of str
        A list of file names with all of the sorted patterns expanded
    """
    for fastq in fastq_globs:
        fastq_list = sorted(glob(str(os.path.expanduser(fastq))))
        if len(fastq_list) == 0:
            logger.warning("No fastqs matched: %s", fastq)
        for filename in fastq_list:
            if os.path.exists(filename):
                yield os.path.abspath(filename)
            else:
                logger.warning("Can't find fastq {}. skipping".format(filename))


def add_trailing_slash(path):
    """Add a trailing slash if not present

    Parameters
    ----------
    path : str
        A string representing a path

    Returns
    -------
    str
        A new string with a trailing / if not previously present.
    """
    if path[-1] != '/':
        path += '/'
    return path


def get_username():
    """Return username

    Return a useful username even if we are running under HT-Condor.

    Returns
    -------
    str : username
    """
    batch_system = os.environ.get('BATCH_SYSTEM')
    if batch_system == 'HTCondor':
        return os.environ.get('USER', '*Unknown user*')

    return os.getlogin()


def get_star_version(star_dir=None):
    if star_dir is not None:
        star_cmd = Path(star_dir) / "STAR"
    else:
        star_cmd = "STAR"
    star = subprocess.run([star_cmd, "--version"], stdout=subprocess.PIPE)
    return star.stdout.decode("utf-8").rstrip()


def get_rsem_version(rsem_dir=None):
    if rsem_dir is not None:
        rsem_cmd = Path(rsem_dir) / "rsem-calculate-expression"
    else:
        rsem_cmd = "rsem-calculate-expression"
    rsem = subprocess.run([rsem_cmd, "--version"], stdout=subprocess.PIPE)
    stdout = rsem.stdout.decode("utf-8").split()
    return stdout[-1]


def sanitize_name(name):
    """Make names more browser and filename safe"""

    return name.replace("/", "_per_").replace(" ", "_").replace(":", "")
