from contextlib import contextmanager
import shutil
import os
from tempfile import TemporaryDirectory
from unittest import TestCase, main
from pkg_resources import resource_filename
from woldrnaseq import models

def count_valid_records(filename):
    total_lines = sum([ 1 for l in open(filename, 'rt') ])
    comment_lines = sum([ 1 for l in open(filename, 'rt') if l.strip().startswith('#')])
    blank_lines = sum([ 1 for l in open(filename, 'rt') if len(l.strip()) == 0 ])
    header_lines = 1
    return total_lines - (comment_lines + header_lines + blank_lines)

class TestModel(TestCase):
    def test_invalid_library(self):
        tsvname = resource_filename(__name__, 'library-invalid.tsv')
        self.assertRaises(ValueError, models.load_library_tables, tsvname)

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

    def test_load_invalid_experiment(self):
        invalid = resource_filename(__name__, 'experiments-invalid.tsv')
        self.assertRaises(ValueError, models.load_experiments, [invalid])

    def test_load_experiment(self):
        mm10tsv = resource_filename(__name__, 'experiments-mm10.tsv')
        hg38tsv = resource_filename(__name__, 'experiments-hg38.tsv')
        mm10 = models.load_experiments([mm10tsv])
        self.assertEqual(len(mm10), count_valid_records(mm10tsv))
        hg38 = models.load_experiments([hg38tsv])
        both = models.load_experiments([mm10tsv, hg38tsv])
        self.assertEqual(len(mm10) + len(hg38), len(both))

    def test_load_library_analysis_root(self):
        with TemporaryDirectory() as analysis_dir:
            print(analysis_dir)
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


@contextmanager
def chdir(d):
    """Change dir in a context manager"""
    olddir = os.getcwd()
    os.chdir(d)
    yield olddir
    os.chdir(olddir)

if __name__ == '__main__':
    main()
