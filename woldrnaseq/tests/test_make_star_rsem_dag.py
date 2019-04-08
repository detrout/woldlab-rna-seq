from unittest import TestCase, main
import logging
import os
import tempfile
from pkg_resources import resource_filename
from woldrnaseq import make_star_rsem_dag
from woldrnaseq import common

class TestMakeDag(TestCase):
    def setUp(self):
        self.home_fd, self.home_pathname = tempfile.mkstemp(
            suffix='wold-rna-seq-home')
        self.old_home_variable = os.environ.get('HOME', tempfile.gettempdir())
        os.environ['HOME'] = self.home_pathname

    def tearDown(self):
        if os.path.isdir(self.home_pathname):
            os.rmdir(self.home_pathname)
        os.environ['HOME'] = self.old_home_variable
        
    def test_arguments_se(self):
        parser = make_star_rsem_dag.make_parser()

        args = parser.parse_args(['-g', 'genome', '-a', 'annotation','--read1', 'a.fastq.gz','b.fastq.gz'])
        with self.assertLogs(common.logger, logging.INFO) as log:
            common.validate_args(args)
        self.assertEqual(len(log.output), 4)

        args = parser.parse_args(['-g', 'genome', '-a', 'annotation',
                                  '--georgi-dir', '~/scripts',
                                  '--genome-dir', '~/genome',
                                  '--read1', 'a.fastq.gz'])
        with self.assertLogs(common.logger, logging.INFO) as log:
            common.validate_args(args)
        self.assertEqual(len(log.output), 2)

    def test_arguments_pe(self):
        parser = make_star_rsem_dag.make_parser()

        args = parser.parse_args(['-g', 'genome', '-a', 'annotation','--read1', 'a.fastq.gz','b.fastq.gz', '--read2', 'a2.fastq.gz','b2.fastq.gz'])
        with self.assertLogs(common.logger, logging.INFO) as log:
            common.validate_args(args)
        self.assertEqual(len(log.output), 4)

        args = parser.parse_args(['-g', 'genome', '-a', 'annotation',
                                  '--georgi-dir', '~/scripts',
                                  '--genome-dir', '~/genome',
                                  '--read1', 'a.fastq.gz',
                                  '--read2', 'b.fastq.gz',
        ])
        with self.assertLogs(common.logger, logging.INFO) as log:
            common.validate_args(args)
        self.assertEqual(len(log.output), 2)

    def test_analysis(self):
        analysis = make_star_rsem_dag.AnalysisDAG()

        analysis.genome_dir = '~/genome_dir'
        analysis.star_dir = '/usr/bin'
        analysis.rsem_dir = '/usr/bin' # should this actually be required?
        analysis.georgi_dir = '~/GeorgiScripts'
        analysis.ucsc_tools_dir = '~/x86_64-304'
        analysis.genome = 'mm10'
        analysis.annotation = 'M4'
        analysis.sex = 'male'
        analysis.job_id = '12345'
        analysis.analysis_dir = '/tmp/12345'
        analysis.read_1_fastqs = ['*.fastqs']

        self.assertTrue(analysis.is_valid())

        dag = str(analysis)

        samtools = resource_filename('woldrnaseq', 'sort-samtools.condor')
        condor_dir = os.path.dirname(samtools)
        for line in str(analysis).split(os.linesep):
            if line.startswith("JOB"):
                self.assertIn(condor_dir, line)
        print(dag)


if __name__ == '__main__':
    main()
