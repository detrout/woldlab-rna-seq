#!/usr/bin/python3
from __future__ import print_function

from argparse import ArgumentParser
import configparser
import os


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    analysis = AnalysisDAG()

    analysis.condor_script_dir = args.condor_script_dir
    analysis.genome_dir = args.genome_dir
    analysis.genome = args.genome
    analysis.annotation = args.annotation
    analysis.sex = args.sex
    analysis.job_id = args.library_id
    analysis.analysis_dir = args.analysis_dir
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
    parser.add_argument('fastqs', nargs='+', help='path to fastqs')

    return parser

def read_defaults():
    defaults = {
        'condor_script_dir': None,
        'genome_dir': None
    }
    config = configparser.ConfigParser()
    config.read([os.path.expanduser('~/.htsworkflow.ini'), '/etc/htsworkflow.ini'])
    if config.has_section('analysis'):
        analysis = config['analysis']
        defaults['condor_script_dir'] = analysis['condor_script_dir']
        defaults['genome_dir'] =analysis['genome_dir']
    return defaults

class AnalysisDAG:
    template = """JOB {job_id}_align-star-se {condor_script_dir}/align-star-se.condor
JOB {job_id}_sort-samtools {condor_script_dir}/sort-samtools.condor
JOB {job_id}_quant-rsem {condor_script_dir}/quant-rsem.condor
JOB {job_id}_index-samtools {condor_script_dir}/index-samtools.condor
JOB {job_id}_qc-samstats {condor_script_dir}/qc-samstats.condor
JOB {job_id}_bedgraph-star {condor_script_dir}/bedgraph-star.condor
JOB {job_id}_qc-coverage {condor_script_dir}/qc-coverage.condor
JOB {job_id}_qc-distribution {condor_script_dir}/qc-distribution.condor
JOB {job_id}_bedgraph2bigwig {condor_script_dir}/bedgraph2bigwig.condor

PARENT {job_id}_align-star-se  CHILD {job_id}_sort-samtools
PARENT {job_id}_align-star-se  CHILD {job_id}_index-samtools
PARENT {job_id}_align-star-se  CHILD {job_id}_bedgraph-star
PARENT {job_id}_index-samtools CHILD {job_id}_qc-samstats
PARENT {job_id}_index-samtools CHILD {job_id}_qc-distribution
PARENT {job_id}_sort-samtools  CHILD {job_id}_quant-rsem
PARENT {job_id}_bedgraph-star  CHILD {job_id}_qc-coverage
PARENT {job_id}_bedgraph-star  CHILD {job_id}_bedgraph2bigwig

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
            job_id=self.job_id,
            genome=self.genome,
            annotation=self.annotation,
            sex=self.sex,
            analysis_dir=self.analysis_dir,
            fastqs=",".join(self.fastqs),
        )


if __name__ == "__main__":
    main()