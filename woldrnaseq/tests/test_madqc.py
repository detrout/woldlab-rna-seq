import os
import tempfile
from pkg_resources import resource_filename
from unittest import TestCase

from woldrnaseq import madqc
from woldrnaseq import models

class TestMadqc(TestCase):
    def setUp(self):
        self.exp_tsv = resource_filename(__name__, 'experiments-mm10.tsv')
        self.lib_tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        self.libraries = models.load_library_tables([self.lib_tsv])
        self.experiments = models.load_experiments([self.exp_tsv])
        
    def test_load_genomic_quantifications(self):
        quantifications = madqc.load_genomic_quantifications(
            self.experiments.ix[0],
            self.libraries,
            'FPKM'
        )

    def test_create_quantification_cache_samedir(self):
        quant = 'FPKM'
        score_filename = models.make_correlation_filename(self.experiments.ix[0])
        quant_filename = models.make_quantification_filename(self.experiments.ix[0],
                                                             quant)

        assert not os.path.exists(score_filename)
        assert not os.path.exists(quant_filename)
        cache = madqc.create_quantification_cache(
            self.experiments.ix[0],
            self.libraries,
            quant)

        self.assertTrue(os.path.exists(score_filename))
        os.remove(score_filename)
        self.assertTrue(os.path.exists(quant_filename))
        os.remove(quant_filename)

    def test_create_quantification_cache_tempdir(self):
        with tempfile.TemporaryDirectory() as tempdir:
            temp_experiments = models.load_experiments([self.exp_tsv],
                                                       analysis_root=tempdir)
            quant = 'FPKM'
            score_filename = models.make_correlation_filename(
                temp_experiments.ix[0])
            quant_filename = models.make_quantification_filename(
                temp_experiments.ix[0], quant)

            print(temp_experiments)
            print(tempdir, score_filename)
            self.assertTrue(score_filename.startswith(tempdir))
            self.assertTrue(quant_filename.startswith(tempdir))
            
            assert not os.path.exists(score_filename)
            assert not os.path.exists(quant_filename)
            cache = madqc.create_quantification_cache(
                temp_experiments.ix[0],
                self.libraries,
                quant)

            self.assertTrue(os.path.exists(score_filename))
            os.remove(score_filename)
            self.assertTrue(os.path.exists(quant_filename))
            os.remove(quant_filename)

            
