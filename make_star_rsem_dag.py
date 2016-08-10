#!/usr/bin/python3
"""STAR & RSEM based RNA-Seq pipeline.
"""
from __future__ import print_function

from argparse import ArgumentParser
import configparser
import os
import logging

logger = logging.getLogger(__name__)

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    analysis = AnalysisDAG()

    analysis.condor_script_dir = args.condor_script_dir
    analysis.genome_dir = args.genome_dir
    analysis.star_dir = args.star_dir
    analysis.rsem_dir = args.rsem_dir
    analysis.georgi_dir = args.georgi_dir
    analysis.genome = args.genome
    analysis.annotation = args.annotation
    analysis.sex = args.sex
    analysis.job_id = args.library_id
    analysis.analysis_dir = args.analysis_dir
    if args.analysis_name is None:
        analysis.analysis_name = os.path.basename(args.analysis_dir)
    analysis.fastqs = args.fastqs

    if analysis.is_valid():
        print(str(analysis))

def make_parser():
    defaults = read_defaults()
    
    parser = ArgumentParser()
    parser.add_argument('--condor-script-dir',
                        help="specify the directory that the condor scripts located in",
                        default=defaults['condor_script_dir'])
    parser.add_argument('--genome-dir',
                        help="specify the directory that has the genome indexes",
                        default=defaults['genome_dir'])
    
    parser.add_argument('-g', '--genome')
    parser.add_argument('-a', '--annotation')
    parser.add_argument('-s', '--sex')
    parser.add_argument('-l', '--library-id')
    parser.add_argument('--analysis-dir', help='target dir to store analysis')
    parser.add_argument('--analysis-name', help='name to store analysis')
    parser.add_argument('fastqs', nargs='+', help='path to fastqs')

    return parser

def add_default_path_arguments(parser):
    """Add arguments to allow overriding location of dependencies
    """
    defaults = read_defaults()
    parser.add_argument('--condor-script-dir',
                        help="specify the directory that the condor scripts located in",
                        default=defaults['condor_script_dir'])
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

    if args.condor_script_dir is None:
        logger.error("Need path to condor templates")
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
        'condor_script_dir': None,
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
        for name in ['condor_script_dir', 'genome_dir',
                     'star_dir', 'rsem_dir', 'georgi_dir']:
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
    template = """JOB {job_id}_align-star-se {condor_script_dir}align-star-se.condor
JOB {job_id}_sort-samtools {condor_script_dir}sort-samtools.condor
JOB {job_id}_quant-rsem {condor_script_dir}quant-rsem.condor
JOB {job_id}_index-samtools {condor_script_dir}index-samtools.condor
JOB {job_id}_qc-samstats {condor_script_dir}qc-samstats.condor
JOB {job_id}_bedgraph-star {condor_script_dir}bedgraph-star.condor
JOB {job_id}_qc-coverage {condor_script_dir}qc-coverage.condor
JOB {job_id}_qc-distribution {condor_script_dir}qc-distribution.condor
JOB {job_id}_bedgraph2bigwig {condor_script_dir}bedgraph2bigwig.condor

PARENT {job_id}_align-star-se  CHILD {job_id}_sort-samtools
PARENT {job_id}_align-star-se  CHILD {job_id}_index-samtools
PARENT {job_id}_align-star-se  CHILD {job_id}_bedgraph-star
PARENT {job_id}_index-samtools CHILD {job_id}_qc-samstats
PARENT {job_id}_index-samtools CHILD {job_id}_qc-distribution
PARENT {job_id}_sort-samtools  CHILD {job_id}_quant-rsem
PARENT {job_id}_bedgraph-star  CHILD {job_id}_qc-coverage
PARENT {job_id}_bedgraph-star  CHILD {job_id}_bedgraph2bigwig

# VARS {job_id}_align-star-se   analysis_name="{analysis_name}"
VARS {job_id}_sort-samtools   analysis_name="{analysis_name}"
VARS {job_id}_quant-rsem      analysis_name="{analysis_name}"
# VARS {job_id}_index-samtools  analysis_name="{analysis_name}"
VARS {job_id}_qc-samstats     analysis_name="{analysis_name}"
VARS {job_id}_bedgraph-star   analysis_name="{analysis_name}"
VARS {job_id}_qc-coverage     analysis_name="{analysis_name}"
VARS {job_id}_qc-distribution analysis_name="{analysis_name}"
VARS {job_id}_bedgraph2bigwig analysis_name="{analysis_name}"

VARS {job_id}_align-star-se   curdir="{analysis_dir}"
VARS {job_id}_sort-samtools   curdir="{analysis_dir}"
VARS {job_id}_quant-rsem      curdir="{analysis_dir}"
VARS {job_id}_index-samtools  curdir="{analysis_dir}"
VARS {job_id}_qc-samstats     curdir="{analysis_dir}"
VARS {job_id}_bedgraph-star   curdir="{analysis_dir}"
VARS {job_id}_qc-coverage     curdir="{analysis_dir}"
VARS {job_id}_qc-distribution curdir="{analysis_dir}"
VARS {job_id}_bedgraph2bigwig curdir="{analysis_dir}"

VARS {job_id}_align-star-se   genome_root="{genome_dir}"
VARS {job_id}_sort-samtools   genome_root="{genome_dir}"
VARS {job_id}_quant-rsem      genome_root="{genome_dir}"
VARS {job_id}_index-samtools  genome_root="{genome_dir}"
VARS {job_id}_qc-samstats     genome_root="{genome_dir}"
VARS {job_id}_bedgraph-star   genome_root="{genome_dir}"
VARS {job_id}_qc-coverage     genome_root="{genome_dir}"
VARS {job_id}_qc-distribution genome_root="{genome_dir}"
VARS {job_id}_bedgraph2bigwig genome_root="{genome_dir}"

VARS {job_id}_align-star-se   star_dir="{star_dir}"
VARS {job_id}_bedgraph-star   star_dir="{star_dir}"

VARS {job_id}_quant-rsem      rsem_dir="{rsem_dir}"

VARS {job_id}_qc-samstats     georgi_dir="{georgi_dir}"
VARS {job_id}_qc-coverage     georgi_dir="{georgi_dir}"
VARS {job_id}_qc-distribution georgi_dir="{georgi_dir}"

VARS {job_id}_align-star-se   genome="{genome}" annotation="{annotation}" sex="{sex}" 
VARS {job_id}_sort-samtools   genome="{genome}" annotation="{annotation}" sex="{sex}" 
VARS {job_id}_quant-rsem      genome="{genome}" annotation="{annotation}" sex="{sex}"
VARS {job_id}_index-samtools  genome="{genome}" annotation="{annotation}" sex="{sex}" 
VARS {job_id}_qc-samstats     genome="{genome}" annotation="{annotation}" sex="{sex}" 
VARS {job_id}_bedgraph-star   genome="{genome}" annotation="{annotation}" sex="{sex}"
VARS {job_id}_qc-coverage     genome="{genome}" annotation="{annotation}" sex="{sex}" 
VARS {job_id}_qc-distribution genome="{genome}" annotation="{annotation}" sex="{sex}" 
VARS {job_id}_bedgraph2bigwig genome="{genome}" annotation="{annotation}" sex="{sex}" 
VARS {job_id}_align-star-se read1="{fastqs}"

"""

    def __init__(self):
        self.condor_script_dir = None
        self.genome_dir = None
        self.star_dir = None
        self.rsem_dir = None
        self.georgi_dir = None
        self.job_id = None
        self.genome = None
        self.annotation = None
        self.sex = None
        self.analysis_dir = None
        self.fastqs = None

    def is_valid(self):
        for key in self.__dict__:
            if getattr(self, key) is None:
                raise ValueError("{} is not set".format(key))
        return True
    
    def __str__(self):
        return self.template.format(
            condor_script_dir=self.condor_script_dir,
            genome_dir=self.genome_dir,
            star_dir=self.star_dir,
            rsem_dir=self.rsem_dir,
            georgi_dir=self.georgi_dir,
            job_id=self.job_id,
            genome=self.genome,
            annotation=self.annotation,
            sex=self.sex,
            analysis_dir=self.analysis_dir,
            analysis_name=self.analysis_name,
            fastqs=",".join(self.fastqs),
        )


if __name__ == "__main__":
    main()
