from unittest import TestCase, main
from unittest.mock import patch
import logging
import os
import re
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
            common.validate_path_args(args)
        self.assertEqual(len(log.output), 4)

        args = parser.parse_args(['-g', 'genome', '-a', 'annotation',
                                  '--georgi-dir', '~/scripts',
                                  '--genome-dir', '~/genome',
                                  '--read1', 'a.fastq.gz'])
        with self.assertLogs(common.logger, logging.INFO) as log:
            common.validate_path_args(args)
        self.assertEqual(len(log.output), 2)

    def test_arguments_pe(self):
        parser = make_star_rsem_dag.make_parser()

        args = parser.parse_args([
            '-g', 'genome',
            '-a', 'annotation',
            '--library-id', '12345',
            '--read1', 'a.fastq.gz', 'b.fastq.gz',
            '--read2', 'a2.fastq.gz', 'b2.fastq.gz'
        ])
        with self.assertLogs(common.logger, logging.INFO) as log:
            common.validate_path_args(args)
        self.assertEqual(len(log.output), 4)

        args = parser.parse_args(['-g', 'genome', '-a', 'annotation',
                                  '--georgi-dir', '~/scripts',
                                  '--genome-dir', '~/genome',
                                  '--read1', 'a.fastq.gz',
                                  '--read2', 'b.fastq.gz',
        ])
        with self.assertLogs(common.logger, logging.INFO) as log:
            common.validate_path_args(args)
        self.assertEqual(len(log.output), 2)

    def test_analysis(self):
        total_size = 1024 ** 3
        with patch('woldrnaseq.make_star_rsem_dag.AnalysisDAG.fastq_size', return_value=total_size):
            analysis = make_star_rsem_dag.AnalysisDAG()

            read_1 = ['22157_GTCATCGT_L001_R1_001.fastq.gz', '22157_GTCATCGT_L001_R1_002.fastq.gz']
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
            analysis.read_1_fastqs = read_1

            self.assertTrue(analysis.is_valid())

            dag = str(analysis)

            samtools = resource_filename('woldrnaseq', 'sort-samtools.condor')
            condor_dir = os.path.dirname(samtools)
            for line in str(analysis).split(os.linesep):
                if line.startswith("JOB"):
                    self.assertIn(condor_dir, line)
                match = re.search('request_disk="(?P<request>[0-9]+)"', line)
                if match:
                    request = match.group(request)
                    self.assertAlmostEqual(float(request), total_size/1024*4)
            print(dag)


    def test_fastq_size(self):
        analysis = make_star_rsem_dag.AnalysisDAG()

        read_1 = ['22157_GTCATCGT_L001_R1_001.fastq.gz', '22157_GTCATCGT_L001_R1_002.fastq.gz']
        read_2 = ['22157_GTCATCGT_L001_R2_001.fastq.gz', '22157_GTCATCGT_L001_R2_002.fastq.gz']
        sizes = {
            '22157_GTCATCGT_L001_R1_001.fastq.gz':   5_000,
            '22157_GTCATCGT_L001_R1_002.fastq.gz':  10_000,
            '22157_GTCATCGT_L001_R2_001.fastq.gz':  20_000,
            '22157_GTCATCGT_L001_R2_002.fastq.gz': 100_000,
        }
        def fake_stat(path):
            if path.endswith('.fastq.gz'):
                return os.stat_result((
                    755,  # st_mode,
                    100,  # st_ino,
                    0,    # st_dev,
                    1,    # st_nlink,
                    1000, # st_uid,
                    1000, # st_guid,
                    sizes[path], # st_size,
                    1596478646,  # st_atime,
                    1596478646,  # st_mtime,
                    1596478646,  # st_ctime
                ))
            else:
                return os.stat(path)
        total_size = sum([sizes[k] for k in sizes])

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
        analysis.read_1_fastqs = read_1
        analysis.read_2_fastqs = read_2

        with patch('os.stat', autospec=True, side_effect=fake_stat):
            self.assertEqual(analysis.fastq_size, total_size)


if __name__ == '__main__':
    main()
