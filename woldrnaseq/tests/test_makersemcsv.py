from pkg_resources import resource_filename
from unittest import TestCase
from unittest.mock import MagicMock, patch

from woldrnaseq import makersemcsv
from woldrnaseq import models


class TestMakeRSEMCSV(TestCase):
    def setUp(self):
        self.exp_tsv = resource_filename(__name__, 'experiments-mm10.tsv')
        self.lib_tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        self.libraries = models.load_library_tables([self.lib_tsv])
        self.experiments = models.load_experiments([self.exp_tsv])
        
    def test_gene_loader_no_annotation(self):
        for quantification_name in ['FPKM', 'TPM', 'expected_count']:
            loader = makersemcsv.GeneRsemLoader(quantification_name, None)
            q = loader.load(self.experiments.loc['expf'], self.libraries)

            self.assertEqual(q.name, 'expf_gene_{}'.format(quantification_name))
            self.assertEqual(q.shape, (144, 2))

    def test_isoform_loader_no_annotation(self):
        for quantification_name in ['FPKM', 'TPM', 'expected_count']:
            loader = makersemcsv.IsoformRsemLoader(quantification_name, None)
            q = loader.load(self.experiments.loc['expf'], self.libraries)

            self.assertEqual(q.name, 'expf_isoform_{}'.format(quantification_name))
            self.assertEqual(q.shape, (144, 2))

    def test_loader_save(self):
        for quantification_name in ['FPKM', 'TPM', 'expected_count']:
            loader = makersemcsv.RsemLoader(quantification_name, None)

            q = MagicMock()
            q.name = 'expf_isoform_{}'.format(quantification_name)

            loader.save(q, 'TAB')
            q.to_csv.assert_called_with(q.name + '.tsv', sep='\t')

            loader.save(q, ',')
            q.to_csv.assert_called_with(q.name + '.csv', sep=',')        

    def test_makersemcsv_invalid_commandline(self):
        # arguments are required
        with self.assertRaises(SystemExit):
            makersemcsv.main([])

    def test_makersemcsv_commandline(self):            
        with patch.object(makersemcsv.GeneRsemLoader, 'load') as load:
            with patch.object(makersemcsv.GeneRsemLoader, 'save') as save:
                q = MagicMock()
                q.name = 'quantification'
                load.return_value = q
                makersemcsv.main(['-e', self.exp_tsv, '-l', self.lib_tsv])

                load.assert_called()
                self.assertEqual(len(load.call_args), 2)
                args, kwargs = load.call_args
                self.assertTrue(self.experiments.ix[-1].equals(args[0]))
                self.assertTrue(self.libraries.equals(args[1]))

                save.assert_called()
                self.assertEqual(len(save.call_args), 2)
                args, kwargs = save.call_args
                self.assertEqual(args[0], q)
                self.assertEqual(args[1], ',')
