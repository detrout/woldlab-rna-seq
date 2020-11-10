"""STAR & RSEM based RNA-Seq pipeline.
"""
from __future__ import absolute_import, print_function

from argparse import ArgumentParser
import datetime
import itertools
import logging
import os
from pkg_resources import resource_filename
import stat
from jinja2 import Environment, PackageLoader

from woldrnaseq import __version__
from .common import (
    add_default_path_arguments,
    add_debug_arguments,
    add_version_argument,
    configure_logging,
)

logger = logging.getLogger(__name__)


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    analysis = AnalysisDAG()

    analysis.genome_dir = args.genome_dir
    analysis.star_dir = args.star_dir
    analysis.rsem_dir = args.rsem_dir
    analysis.georgi_dir = args.georgi_dir
    analysis.ucsc_tools_dir = args.ucsc_tools_dir
    analysis.genome = args.genome
    analysis.annotation = args.annotation
    analysis.sex = args.sex
    analysis.job_id = args.library_id
    analysis.analysis_dir = args.analysis_dir
    analysis.analysis_name = args.analysis_name
    analysis.read_1_fastqs = args.read1
    analysis.read_2_fastqs = args.read2
    analysis.splice_template = args.splice_template

    if analysis.is_valid():
        print(str(analysis))


def make_parser():
    parser = ArgumentParser()
    parser.add_argument('-g', '--genome')
    parser.add_argument('-a', '--annotation')
    parser.add_argument('-s', '--sex')
    parser.add_argument('-l', '--library-id')
    parser.add_argument('--analysis-dir', help='target dir to store analysis')
    parser.add_argument('--analysis-name', help='name to store analysis')
    parser.add_argument('--read1', nargs='+', help='path to read 1 fastqs')
    parser.add_argument('--read2', nargs='*', default=[],
                        help='path to read 2 fastqs')
    parser.add_argument('--splice-template', default='star_rsem.dagman',
                        help='Override splice dagman template')

    add_default_path_arguments(parser)
    add_version_argument(parser)
    add_debug_arguments(parser)

    return parser


class AnalysisDAG:
    def __init__(self):
        self.genome_dir = None
        self.star_dir = None
        self.rsem_dir = None
        self.georgi_dir = None
        self.ucsc_tools_dir = None
        self.job_id = None
        self.genome = None
        self.annotation = None
        self.sex = None
        self.analysis_dir = None
        self._analysis_name = None
        self.read_1_fastqs = []
        self.read_2_fastqs = []
        self.stranded = 'unstranded'
        self.reference_prefix = 'chr'
        self.splice_template = 'star_rsem.dagman'

    @property
    def fastq_size(self):
        filesize = 0
        for filename in itertools.chain(self.read_1_fastqs, self.read_2_fastqs):
            filesize += os.stat(filename)[stat.ST_SIZE]
        return filesize

    def is_valid(self):
        for key in self.__dict__:
            if key == 'read_1_fastqs' and len(self.read_1_fastqs) == 0:
                raise ValueError("Read 1 fastqs are required for library {}".format(self.job_id))
            elif key == '_analysis_name':
                # analysis name will default to analysis dir
                pass
            elif key == 'analysis_dir':
                if self.analysis_dir is None:
                    raise ValueError("analysis_dir is not set")
                if not os.path.exists(self.analysis_dir):
                    logger.warning("Analysis dir %s doesn't exist", self.analysis_dir)
            elif getattr(self, key) is None:
                raise ValueError("{} is not set".format(key))
        return True

    @property
    def analysis_name(self):
        if self._analysis_name is not None:
            return self._analysis_name

        if self.analysis_dir is not None:
            self._analysis_name = os.path.basename(self.analysis_dir)

        return self._analysis_name

    @analysis_name.setter
    def analysis_name(self, value):
        self._analysis_name = value

    def __str__(self):
        env = Environment(loader=PackageLoader('woldrnaseq', 'templates'))
        template = env.get_template(self.splice_template)

        rsem_paired_argument = '--paired-end' if len(self.read_2_fastqs) > 0 else ''
        if self.stranded == 'forward':
            rsem_strand_probability = '--forward-prob 1'
        elif self.stranded == 'reverse':
            rsem_strand_probability = '--forward-prob 0'
        else:
            rsem_strand_probability = '--forward-prob 0.5'

        return template.render(
            align_star=resource_filename(__name__, 'align-star.condor'),
            pre_star=resource_filename(__name__, 'pre_star'),
            post_star=resource_filename(__name__, 'post_star'),
            sort_samtools=resource_filename(__name__, 'sort-samtools.condor'),
            quant_rsem=resource_filename(__name__, 'quant-rsem.condor'),
            index_samtools=resource_filename(__name__, 'index-samtools.condor'),
            qc_samstats=resource_filename(__name__, 'qc-samstats.condor'),
            bedgraph_star=resource_filename(__name__, 'bedgraph-star.condor'),
            qc_coverage=resource_filename(__name__, 'qc-coverage.condor'),
            qc_distribution=resource_filename(__name__, 'qc-distribution.condor'),
            bedgraph2bigwig=resource_filename(__name__, 'bedgraph2bigwig.condor'),
            sort_samtools_sh=resource_filename(__name__, 'sort-samtools.sh'),
            bedgraph_bedsort_sh=resource_filename(__name__, 'bedsort.sh'),
            picard_markdup=resource_filename(__name__, 'picard-markdup.condor'),
            picard_multi_metrics=resource_filename(__name__, 'picard-multi-metrics.condor'),
            rrna_premap=resource_filename(__name__, 'rrna-premap.condor'),
            rrna_premap_sh=resource_filename(__name__, 'rrna-premap.sh'),
            genome_dir=self.genome_dir,
            star_dir=self.star_dir,
            rsem_dir=self.rsem_dir,
            georgi_dir=self.georgi_dir,
            ucsc_tools_dir=self.ucsc_tools_dir,
            job_id=self.job_id,
            genome=self.genome,
            annotation=self.annotation,
            sex=self.sex,
            analysis_dir=self.analysis_dir,
            analysis_name=self.analysis_name,
            read_1_fastqs=",".join(self.read_1_fastqs),
            read_2_fastqs=",".join(self.read_2_fastqs),
            star_request_memory_megabytes=int(self.fastq_size/(1024**2) * 3),
            star_request_memory_bytes=int(self.fastq_size * 3),
            star_request_disk_kilobytes=int(self.fastq_size/1024 * 4),
            reference_prefix=self.reference_prefix,
            rsem_paired_argument=rsem_paired_argument,
            rsem_strand_probability=rsem_strand_probability,
            rsem_request_disk=int(self.fastq_size/1024 * 7),
            username=os.getlogin(),
            timestamp=datetime.datetime.now().isoformat(),
            woldrnaseq_version=__version__,
        )


if __name__ == "__main__":
    main()
