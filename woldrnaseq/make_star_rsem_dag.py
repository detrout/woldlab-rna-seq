"""STAR & RSEM based RNA-Seq pipeline.
"""
from __future__ import print_function

from argparse import ArgumentParser
import configparser
import os
import logging
from pkg_resources import resource_filename
from jinja2 import Environment, PackageLoader

logger = logging.getLogger(__name__)

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    analysis = AnalysisDAG()

    analysis.genome_dir = args.genome_dir
    analysis.star_dir = args.star_dir
    analysis.rsem_dir = args.rsem_dir
    analysis.georgi_dir = args.georgi_dir
    analysis.genome = args.genome
    analysis.annotation = args.annotation
    analysis.sex = args.sex
    analysis.job_id = args.library_id
    analysis.analysis_dir = args.analysis_dir
    analysis.read_1_fastqs = args.read1
    analysis.read_2_fastqs = args.read2

    if analysis.is_valid():
        print(str(analysis))

def make_parser():
    parser = ArgumentParser()
    parser.add_argument('-g', '--genome')
    parser.add_argument('-a', '--annotation')
    parser.add_argument('-s', '--sex')
    parser.add_argument('-l', '--library-id')
    parser.add_argument('--analysis-dir', help='target dir to store analysis')
    parser.add_argument('--read1', nargs='+', help='path to read 1 fastqs')
    parser.add_argument('--read2', nargs='*', default=[],
                        help='path to read 2 fastqs')

    add_default_path_arguments(parser)
    add_debug_arguments(parser)

    return parser

def add_default_path_arguments(parser):
    """Add arguments to allow overriding location of dependencies
    """
    defaults = read_defaults()
    parser.add_argument('--genome-dir',
                        help="specify the directory that has the genome indexes",
                        default=defaults['genome_dir'])
    parser.add_argument('--star-dir',
                        default=defaults['star_dir'],
                        help='Specify the directory where STAR is installed')
    parser.add_argument('--rsem-dir',
                        default=defaults['rsem_dir'],
                        help='Specify the directory where rsem is installed')
    parser.add_argument('--georgi-dir',
                        default=defaults['georgi_dir'],
                        help='Specify the directory where georgi scripts are installed')
    return parser

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

def add_debug_arguments(parser):
    """Add arguments for tuning logging
    """
    parser.add_argument('-v', '--verbose', default=False, action='store_true')
    parser.add_argument('-d', '--debug', default=False, action='store_true')
    return parser

def read_defaults():
    defaults = {
        'genome_dir': None,
        'star_dir': None,
        'rsem_dir': None,
        'georgi_dir': None,
    }
    config = configparser.ConfigParser()
    config.read([os.path.expanduser('~/.htsworkflow.ini'),
                 '/etc/htsworkflow.ini'])

    if config.has_section('analysis'):
        analysis = config['analysis']
        for name in ['genome_dir', 'star_dir', 'rsem_dir', 'georgi_dir']:
            defaults[name] = normalize_path(analysis.get(name))

    return defaults

def normalize_path(path):
    if path is None:
        return None
    elif len(path) == 0:
        return path
    else:
        return os.path.join(path, '')

class AnalysisDAG:
    def __init__(self):
        self.genome_dir = None
        self.star_dir = None
        self.rsem_dir = None
        self.georgi_dir = None
        self.job_id = None
        self.genome = None
        self.annotation = None
        self.sex = None
        self.analysis_dir = None
        self.read_1_fastqs = []
        self.read_2_fastqs = []

    def is_valid(self):
        for key in self.__dict__:
            if key == 'read_1_fastqs' and len(self.read_1_fastqs) == 0:
                raise ValueError("Read 1 fastqs are required")
            elif getattr(self, key) is None:
                raise ValueError("{} is not set".format(key))
        return True
    
    def __str__(self):
        env = Environment(loader=PackageLoader('woldrnaseq', 'templates'))
        template = env.get_template('star_rsem.dagman')

        return template.render(
            align_star=resource_filename(__name__, 'align-star.condor'),
            sort_samtools=resource_filename(__name__, 'sort-samtools.condor'),
            quant_rsem=resource_filename(__name__, 'quant-rsem.condor'),
            index_samtools=resource_filename(__name__, 'index-samtools.condor'),
            qc_samstats=resource_filename(__name__, 'qc-samstats.condor'),
            bedgraph_star=resource_filename(__name__, 'bedgraph-star.condor'),
            qc_coverage=resource_filename(__name__, 'qc-coverage.condor'),
            qc_distribution=resource_filename(__name__, 'qc-distribution.condor'),
            bedgraph2bigwig=resource_filename(__name__, 'bedgraph2bigwig.condor'),
            sort_samtools_sh=resource_filename(__name__, 'sort-samtools.sh'),
            genome_dir=self.genome_dir,
            star_dir=self.star_dir,
            rsem_dir=self.rsem_dir,
            georgi_dir=self.georgi_dir,
            job_id=self.job_id,
            genome=self.genome,
            annotation=self.annotation,
            sex=self.sex,
            analysis_dir=self.analysis_dir,
            read_1_fastqs=",".join(self.read_1_fastqs),
            read_2_fastqs=",".join(self.read_2_fastqs),
        )


if __name__ == "__main__":
    main()
