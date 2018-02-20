import os
from pkg_resources import resource_filename
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

import pandas

from woldrnaseq import models
from woldrnaseq import gtfcache


class TestGTFCache(TestCase):
    def setUp(self):
        self.mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        self.mm10 = models.load_library_tables([self.mm10tsv])

        self.female = pandas.DataFrame()
        self.female.index.name = 'female'
        
        self.male = pandas.DataFrame()
        self.male.index.name = 'male'

    def test_gtf_cache_mapping(self):
        """Does GTFCache properly behave as a mapping?
        """
        with TemporaryDirectory('gtfcache') as tempdir:
            cache = gtfcache.GTFCache(self.mm10, tempdir)

            # test __len__
            uncommented_libraries = 6
            self.assertEqual(len(cache), uncommented_libraries)

            # test __iter__
            for l, r in zip(self.mm10.index, cache):
                self.assertEqual(l, r)

            cache._gtf_cache['mm10-M4-female'] = self.female
            self.assertEqual(cache['12304'].index.name, 'female')

            cache._gtf_cache['mm10-M4-male'] = self.male
            self.assertEqual(cache['12307'].index.name, 'male')

    def test_gtf_cache_filename(self):
        cache = gtfcache.GTFCache(self.mm10, '/tmp')

        filename = cache._get_gtf_cache_filename('genome')
        self.assertEqual(filename, os.path.join('/tmp', 'genome', 'genome.h5'))

    def test_get_gtf_cache(self):
        with TemporaryDirectory('gtfcache') as tempdir:
            with patch('woldrnaseq.gtfcache.load_gtf_cache') as load_gtf_cache:
                cache = gtfcache.GTFCache(self.mm10, tempdir)

                cache._gtf_cache['mm10-M4-female'] = self.female
                self.assertEqual(cache['12304'].index.name, 'female')
                self.assertFalse(load_gtf_cache.called)

                load_gtf_cache.return_value = self.male
                df = cache['12307']
                load_gtf_cache.called_with('mm10-M4-male')
                self.assertEqual(df.index.name, 'male')
