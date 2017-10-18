"""Common utilities shared by several components
"""
import configparser
from glob import glob
import logging
import os

import pandas
import numpy

logger = logging.getLogger(__name__)

def add_default_path_arguments(parser):
    """Add arguments to allow overriding location of dependencies
    """
    defaults = read_defaults()
    group = parser.add_argument_group('Paths', 'paths to needed external resources')
    group.add_argument('--genome-dir',
                        help="specify the directory that has the genome indexes",
                        default=defaults['genome_dir'])
    group.add_argument('--star-dir',
                        default=defaults['star_dir'],
                        help='Specify the directory where STAR is installed')
    group.add_argument('--rsem-dir',
                        default=defaults['rsem_dir'],
                        help='Specify the directory where rsem is installed')
    group.add_argument('--georgi-dir',
                        default=defaults['georgi_dir'],
                        help='Specify the directory where georgi scripts are installed')
    group.add_argument('--ucsc-tools-dir',
                       default=defaults['ucsc_tools_dir'],
                       help='Specify the directory where the UCSC tools are installed')
    return parser


def add_debug_arguments(parser):
    """Add arguments for tuning logging
    """
    group = parser.add_argument_group('Verbosity', 'Set logging level')
    parser.add_argument('-v', '--verbose', default=False, action='store_true')
    parser.add_argument('-d', '--debug', default=False, action='store_true')
    return parser


def add_version_argument(parser):
    """Add Version argument
    """
    parser.add_argument('--version', default=False, action='store_true')
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
    config.read([os.path.expanduser('~/.htsworkflow.ini'),
                 '/etc/htsworkflow.ini'])

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
        return os.path.join(path, '')


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


def validate_args(args):
    """Warn if path arguments aren't set
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


def find_fastqs(table, fastq_column):
    """Find fastqs for a library from a library table

    fastqs are a comma seperated glob pattern
    """
    if fastq_column in table.columns:
        for library_id in table.index:
            fastq_field = table.loc[library_id, fastq_column]
            if not pandas.isnull(fastq_field):
                fastqs = find_fastqs_by_glob(fastq_field.split(','))
                yield (library_id, list(fastqs))
    else:
        # eventually look up by library ID
        raise NotImplemented("Please specify fastq glob")


def find_fastqs_by_glob(fastq_globs):
    for fastq in fastq_globs:
        fastq_list = sorted(glob(fastq))
        if len(fastq_list) == 0:
            logger.warn("No fastqs matched: %s", fastq)
        for filename in fastq_list:
            if os.path.exists(filename):
                yield os.path.abspath(filename)
            else:
                logger.warn("Can't find fastq {}. skipping".format(filename))


