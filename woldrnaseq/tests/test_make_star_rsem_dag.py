from unittest import TestCase, main
import logging
import os
from pkg_resources import resource_filename
from woldrnaseq import make_star_rsem_dag

class TestMakeDag(TestCase):
    def test_arguments(self):
        parser = make_star_rsem_dag.make_parser()

        args = parser.parse_args(['-g', 'genome', '-a', 'annotation','a.fastq.gz','b.fastq.gz'])
        with self.assertLogs(make_star_rsem_dag.logger, logging.INFO) as log:
            make_star_rsem_dag.validate_args(args)
        self.assertEqual(len(log.output), 4)

        args = parser.parse_args(['-g', 'genome', '-a', 'annotation',
                                  '--georgi-dir', '~/scripts',
                                  '--genome-dir', '~/genome',
                                  'a.fastq.gz'])
        with self.assertLogs(make_star_rsem_dag.logger, logging.INFO) as log:
            make_star_rsem_dag.validate_args(args)
        self.assertEqual(len(log.output), 2)

    def test_analysis(self):
        analysis = make_star_rsem_dag.AnalysisDAG()

        analysis.genome_dir = '~/genome_dir'
        analysis.star_dir = '/usr/bin'
        analysis.rsem_dir = '/usr/bin' # should this actually be required?
        analysis.georgi_dir = '~/GeorgiScripts'
        analysis.genome = 'mm10'
        analysis.annotation = 'M4'
        analysis.sex = 'male'
        analysis.job_id = '12345'
        analysis.analysis_dir = '/tmp/12345'
        analysis.fastqs = ['*.fastqs']

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
