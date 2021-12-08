from contextlib import contextmanager
import shutil
import os
from tempfile import TemporaryDirectory
from unittest import TestCase, main
from pkg_resources import resource_filename
import pandas
import pytest
from woldrnaseq import models


def count_valid_records(filename):
    total_lines = sum([1 for l in open(filename, 'rt')])
    comment_lines = sum([1 for l in open(filename, 'rt') if l.strip().startswith('#')])
    blank_lines = sum([1 for l in open(filename, 'rt') if len(l.strip()) == 0])
    header_lines = 1
    return total_lines - (comment_lines + header_lines + blank_lines)


class TestModel(TestCase):
    def test_read_line_from_stream(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        with open(mm10tsv) as instream:
            lines = list(models.read_line_from_stream(instream))
        mm10 = models.load_library_tables([mm10tsv])
        # add one to mm10 dataframe because the header is not counted in len()
        self.assertEqual(len(lines), len(mm10) + 1)

    def test_verify_library_columns_missing_required(self):
        df = pandas.DataFrame({'library_id': ['1'],
                               'genome': ['mm10'],
                               'sex': ['male'],
                               'annotation': ['M4'],
                               'analysis_dir': ['1/']})
        df.set_index('library_id', inplace=True)

        self.assertRaises(ValueError, models.verify_library_columns, df)

    def test_verify_library_columns_optional_present(self):
        df = pandas.DataFrame({'library_id': ['1'],
                               'genome': ['mm10'],
                               'sex': ['male'],
                               'annotation': ['M4'],
                               'analysis_dir': ['1/'],
                               'stranded': ['forward'],
                               'read_1': ['R1.fastq.gz'],
                               'read_2': ['R2.fastq.gz']})
        df.set_index('library_id', inplace=True)

        models.verify_library_columns(df)

    def test_verify_library_columns_misspelled_optional(self):
        df = pandas.DataFrame({'library_id': ['1'],
                               'genome': ['mm10'],
                               'sex': ['male'],
                               'annotation': ['M4'],
                               'analysis_dir': ['1/'],
                               'read_1': ['R1.fastq.gz'],
                               'read2': ['R2.fastq.gz']})
        df.set_index('library_id', inplace=True)

        with self.assertLogs(models.__name__, level='WARNING') as log_manager:
            models.verify_library_columns(df)
        self.assertEqual(
            log_manager.output,
            ['WARNING:{}:Unrecognized columns present. Is this intended?: read2'.format(
                models.__name__)])

    def test_invalid_library(self):
        tsvname = resource_filename(__name__, 'library-invalid.tsv')
        self.assertRaises(ValueError, models.load_library_tables, [tsvname])

    def test_load_library(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        hg38tsv = resource_filename(__name__, 'library-hg38-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])
        self.assertEqual(len(mm10), count_valid_records(mm10tsv))
        hg38 = models.load_library_tables([hg38tsv])
        both = models.load_library_tables([mm10tsv, hg38tsv])
        self.assertEqual(len(mm10) + len(hg38), len(both))

    def test_load_duplicate_libraries(self):
        hg38tsv = resource_filename(__name__, 'library-hg38-se.tsv')
        self.assertRaises(ValueError, models.load_library_tables, [hg38tsv, hg38tsv])

    def test_load_stranded_library(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-stranded.tsv')
        mm10 = models.load_library_tables([mm10tsv])
        expected = ['forward', 'reverse', 'unstranded', 'forward', 'reverse', 'unstranded']
        for strand, (library_id, row) in zip(expected, mm10.iterrows()):
            self.assertEqual(strand, row.stranded)

    def test_load_stranded_invalid_library(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-stranded-invalid.tsv')
        self.assertRaises(ValueError, models.load_library_tables, [mm10tsv])

    def test_load_invalid_experiment(self):
        invalid = resource_filename(__name__, 'experiments-invalid.tsv')
        self.assertRaises(ValueError, models.load_experiments, [invalid])

    def test_load_numeric_experiment(self):
        filename = resource_filename(__name__, 'experiments-numeric.tsv')
        experiment = models.load_experiments([filename])
        for name in experiment.index:
            self.assertIsInstance(name, str)

    def test_load_experiment(self):
        mm10tsv = resource_filename(__name__, 'experiments-mm10.tsv')
        hg38tsv = resource_filename(__name__, 'experiments-hg38.tsv')
        mm10 = models.load_experiments([mm10tsv])
        self.assertIn('replicates', mm10.columns)
        self.assertEqual(len(mm10), count_valid_records(mm10tsv))
        hg38 = models.load_experiments([hg38tsv])
        both = models.load_experiments([mm10tsv, hg38tsv])
        self.assertEqual(len(mm10) + len(hg38), len(both))

        self.assertEqual(mm10.loc['expm']['replicates'],
                         ['12307', '12308'])

    def test_load_experiments_analysis_root(self):
        with TemporaryDirectory() as analysis_dir:
            with chdir(analysis_dir):
                mm10tsv = resource_filename(__name__, 'experiments-mm10.tsv')
                tmpname = os.path.join(analysis_dir, 'experiments-mm10.tsv')
                shutil.copy(mm10tsv, tmpname)
                analysis_root = os.path.dirname(mm10tsv)
                mm10 = models.load_experiments([mm10tsv])
                mm10tmp = models.load_experiments([tmpname],
                                                  analysis_root=analysis_root)
                for i in mm10['analysis_dir'].index:
                    self.assertEqual(mm10['analysis_dir'][i],
                                     mm10tmp['analysis_dir'][i])

    def test_load_library_analysis_root(self):
        with TemporaryDirectory() as analysis_dir:
            with chdir(analysis_dir):
                mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
                tmpname = os.path.join(analysis_dir, 'library-mm10-se.tsv')
                shutil.copy(mm10tsv, tmpname)
                analysis_root = os.path.dirname(mm10tsv)
                mm10 = models.load_library_tables([mm10tsv])
                mm10tmp = models.load_library_tables([tmpname],
                                                     analysis_root=analysis_root)
                for i in mm10['analysis_dir'].index:
                    self.assertEqual(mm10['analysis_dir'][i],
                                     mm10tmp['analysis_dir'][i])

    def test_load_find_library_analysis_file(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])

        cwd_files = list(models.find_library_analysis_file(mm10, '*.coverage'))
        self.assertGreaterEqual(len(cwd_files), 1)
        for f in cwd_files:
            self.assertTrue(isinstance(f, models.AnalysisFile))

        with TemporaryDirectory() as analysis_dir:
            with chdir(analysis_dir):
                mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
                tmpname = os.path.join(analysis_dir, 'library-mm10-se.tsv')
                shutil.copy(mm10tsv, tmpname)
                analysis_root = os.path.dirname(mm10tsv)
                mm10 = models.load_library_tables([tmpname],
                                                  analysis_root=analysis_root)

                abs_files = list(models.find_library_analysis_file(mm10, '*.coverage'))
                self.assertGreaterEqual(len(abs_files), 1)
                for f in abs_files:
                    self.assertTrue(isinstance(f, models.AnalysisFile))

        self.assertEqual(len(cwd_files), len(abs_files))
        self.assertEqual(cwd_files[0].filename, abs_files[0].filename)

    def test_load_all_coverage(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])
        coverage = models.load_all_coverage(mm10)
        self.assertEqual(coverage.shape, (100, 1))

    def test_load_flagstat_ENCFF284YOU(self):
        filename = resource_filename(__name__, 'ENCFF284YOU_flagstat.txt')
        flagstat = models.load_flagstat(filename)
        expected = {
            'duplicates': 0,
            'total': 51759944,
            'mapped_qc_failed': 3779026,
            'mapped': 49892660,
            'mapped_pct': '96.39%',
            'total_qc_failed': 4818581,
            'duplicates_qc_failed': 0,
        }
        for k in expected:
            self.assertEqual(flagstat[k], expected[k])
        self.assertEqual(sorted(expected), sorted(flagstat))

    def test_load_flagstat_ENCFF704SWF(self):
        filename = resource_filename(__name__, 'ENCFF704SWF_flagstat.txt')
        flagstat = models.load_flagstat(filename)
        expected = {
            'paired_properly': 192704068,
            'paired_properly_pct': '90.01%',
            'paired': 214087832,
            'with_itself_qc_failed': 0,
            'with_itself': 192704068,
            'paired_properly_qc_failed': 0,
            'read2_qc_failed': 0,
            'total_qc_failed': 0,
            'mapped_qc_failed': 0,
            'duplicates': 0,
            'diff_chroms': 0,
            'diff_chroms_qc_failed': 0,
            'singletons_pct': '0.00%',
            'singletons': 0,
            'mapped': 192704068,
            'read1': 107043916,
            'total': 214087832,
            'duplicates_qc_failed': 0,
            'read1_qc_failed': 0,
            'mapped_pct': '90.01%',
            'singletons_qc_failed': 0,
            'read2': 107043916,
        }
        for k in expected:
            self.assertEqual(
                flagstat[k], expected[k],
                '{} did not match'.format(k),
            )
        self.assertEqual(sorted(expected), sorted(flagstat))

    def test_load_all_samstats(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])
        samstats = models.load_all_samstats(mm10)
        self.assertEqual(samstats.shape, (1, 10))
        self.assertEqual(samstats.index[0], '12304')

    def test_load_one_star_count(self):
        expected = {
            'U': {
                'ENSG00000223972.5': 0,
                'ENSG00000239945.1': 9,
            },
            '+': {
                'ENSG00000223972.5': 0,
                'ENSG00000239945.1': 7,
            },
            '-': {
                'ENSG00000223972.5': 0,
                'ENSG00000239945.1': 2,
            },
        }
        lib12304 = resource_filename(__name__, 'lib_12304/ReadsPerGene.out.tab')

        for column in expected:
            star_count = models.load_star_counts(lib12304, column)
            self.assertEqual(star_count.shape, (11, 1))
            for gene in expected[column]:
                self.assertEqual(star_count.loc[gene][column], expected[column][gene])

    def test_load_all_star_counts(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])
        scores = models.load_all_star_counts(mm10, '+')
        self.assertEqual(scores.shape, (11, 2))
        self.assertEqual(scores.index.name, 'gene_id')
        self.assertEqual(list(scores.columns), ['12304', '12305'])

    def test_load_star_final_log(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])
        final_out = models.load_all_star_final(mm10)
        self.assertEqual(final_out.shape, (1, 27))
        self.assertEqual(final_out.index[0], '12304')

    @pytest.mark.xfail
    def test_load_star_solo_quality_metric(self):
        raise NotImplementedError("Need to wait until I have a star release")

    def test_load_all_distribution(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])
        distribution = models.load_all_distribution(mm10)
        self.assertEqual(distribution.shape, (1, 3))
        self.assertEqual(distribution.index[0], '12304')

    def test_warn_if_spaces(self):
        mm10invalid = resource_filename(__name__, 'library-invalid.tsv')

        with self.assertLogs(models.logger, level='WARN') as log:
            models.warn_if_spaces(mm10invalid)
        self.assertEqual(log.output,
                         ['WARNING:woldrnaseq.models:There are spaces in the header line, is this intentional?'])

    def test_make_correlation_filename_default(self):
        """Does make_correlation_filename work with default analysis_root
        """
        results = {
            'genome': 'expf_correlation.h5',
            'transcriptome': 'expf_transcriptome_correlation.h5',
        }
        for reference_type in results:
            mm10tsv = resource_filename(__name__, 'experiments-mm10.tsv')
            path, _ = os.path.split(mm10tsv)
            mm10 = models.load_experiments([mm10tsv])

            filename = models.make_correlation_filename(
                mm10.iloc[0],
                reference_type=reference_type,
            )
            expected = os.path.join(path, results[reference_type])
            self.assertEqual(filename, expected)

    def test_make_correlation_filename_other(self):
        """Does make_correlation_filename work with an alternate analysis_root
        """
        results = {
            'genome': 'expf_correlation.h5',
            'transcriptome': 'expf_transcriptome_correlation.h5',
        }
        for reference_type in results:
            mm10tsv = resource_filename(__name__, 'experiments-mm10.tsv')
            path = '/tmp'
            mm10 = models.load_experiments([mm10tsv], analysis_root=path)

            filename = models.make_correlation_filename(
                mm10.iloc[0],
                reference_type=reference_type,
            )
            expected = os.path.join(path, results[reference_type])
            self.assertEqual(filename, expected)

    def test_make_quantification_filename_default(self):
        """Does make_quantification_filename work with default analysis_root
        """
        results = {
            'genome': 'expf_FPKM.h5',
            'transcriptome': 'expf_transcriptome_FPKM.h5',
        }
        for reference_type in results:
            mm10tsv = resource_filename(__name__, 'experiments-mm10.tsv')
            path, _ = os.path.split(mm10tsv)
            mm10 = models.load_experiments([mm10tsv])

            filename = models.make_quantification_filename(
                mm10.iloc[0],
                reference_type=reference_type,
            )
            expected = os.path.join(path, results[reference_type])
            self.assertEqual(filename, expected)

    def test_make_quantification_filename_other(self):
        """Does make_quantification_filename work with an alternate analysis_root
        """
        results = {
            'genome': 'expf_FPKM.h5',
            'transcriptome': 'expf_transcriptome_FPKM.h5',
        }
        for reference_type in results:
            mm10tsv = resource_filename(__name__, 'experiments-mm10.tsv')
            path = '/tmp'
            mm10 = models.load_experiments([mm10tsv], analysis_root=path)

            filename = models.make_quantification_filename(
                mm10.iloc[0],
                reference_type=reference_type,
            )
            expected = os.path.join(path, results[reference_type])
            self.assertEqual(filename, expected)

    def test_genome_name_from_library_series(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])

        self.assertEqual(models.genome_name_from_library(mm10.loc['12304']), 'mm10-M4-female')
        self.assertEqual(models.genome_name_from_library(mm10.loc['12309']), 'mm10-M4-male')

    def test_genome_name_from_library_dict(self):
        d = {
            'genome': 'mm10',
            'annotation': 'M21_minimal',
            'sex': 'male',
        }

        self.assertEqual(models.genome_name_from_library(d), 'mm10-M21_minimal-male')
        self.assertRaises(ValueError, models.genome_name_from_library, 10)


@contextmanager
def chdir(d):
    """Change dir in a context manager"""
    olddir = os.getcwd()
    os.chdir(d)
    yield olddir
    os.chdir(olddir)


if __name__ == '__main__':
    main()
