import unittest
from pkg_resources import resource_filename
from io import StringIO
import numpy

from woldrnaseq.gff2table import (
    AttributesParser,
    GFFParser,
    parse_score,
    parse_strand,
    parse_phase,
)


class TestAttributeParser(unittest.TestCase):
    def test_quoting_missing(self):
        p = AttributesParser(' ')
        self.assertRaises(ValueError, list, p.tokenize('a "b'))

    def test_embedded_sep(self):
        p = AttributesParser(' ')
        t = list(p.tokenize('a "b;c"; d 1'))
        self.assertEqual(t, ['a', ' ', '"b;c"', ';', 'd', ' ', '1'])

    def test_empty(self):
        p = AttributesParser()
        attributes = p('')
        self.assertEqual(attributes, 0)

    def test_gtf_extra_spaces(self):
        p = AttributesParser(' ')
        t = list(p.tokenize(' a 1;    b   "b"; '))
        self.assertEqual(t, ['a', ' ', '1', ';', 'b', ' ', '"b"', ';'])

    def test_gff_extra_spaces(self):
        p = AttributesParser('=')
        t = list(p.tokenize(' a=1;    b=  "b"; '))
        self.assertEqual(t, ['a', '=', '1', ';', 'b', '=', '"b"', ';'])

    @unittest.skip('We should suppress extra semicolons')
    def test_gff_extra_semicolons(self):
        p = AttributesParser('=')
        t = list(p.tokenize(' a=1;; b="b"; '))
        self.assertEqual(t, ['a', '=', '1', ';', 'b', '=', '"b"', ';'])
        
    def test_simple_gtf(self):
        p = AttributesParser()
        attributes = p('gene_id "gene1"')
        self.assertEqual(attributes, 1)
        self.assertIn('gene_id', p.terms)
        self.assertEqual(p.terms['gene_id'], {0: 'gene1'})
        self.assertIsInstance(p.terms['gene_id'][0], str)

    def test_numeric_gtf(self):
        p = AttributesParser()
        attributes = p('exon_number 1')
        self.assertEqual(attributes, 1)
        self.assertIn('exon_number', p.terms)
        self.assertEqual(p.terms['exon_number'], {0: 1})
        self.assertIsInstance(p.terms['exon_number'][0], int)
        
    def test_multi_value_gtf(self):
        p = AttributesParser(' ')
        attributes = p(
            'gene_id "gene1"; transcript_id "transcript1"; exon_number 1')
        self.assertEqual(attributes, 3)
        line_no = 0
        for name, expected, attribute_type in [
                ('gene_id', 'gene1', str),
                ('transcript_id', 'transcript1', str),
                ('exon_number', 1, int)]:
            
            self.assertIn(name, p.terms)
            self.assertEqual(p.terms[name], {line_no: expected})
            self.assertIsInstance(p.terms[name][line_no], attribute_type)

    def test_multi_value_gff(self):
        p = AttributesParser('=')
        attributes = p('gene_id=gene1;transcript_id=transcript1;exon_number=1')
        self.assertEqual(attributes, 3)
        line_no = 0
        for name, expected, attribute_type in [
                ('gene_id', 'gene1', str),
                ('transcript_id', 'transcript1', str),
                ('exon_number', 1, int)]:

            self.assertIn(name, p.terms)
            self.assertEqual(p.terms[name], {line_no: expected})
            self.assertIsInstance(p.terms[name][line_no], attribute_type)

    def test_gff_value_with_spaces(self):
        p = AttributesParser('=')
        t = list(p.tokenize('mol_type=genomic DNA'))
        self.assertEqual(t, ['mol_type', '=', 'genomic DNA'])

    def test_grcm38(self):
        p = AttributesParser('=')
        attributes = p('ID=id0;Dbxref=taxon:10090;Name=1;chromosome=1;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=C57BL/6J')
        self.assertEqual(attributes, 8)
        line_no = 0
        for name, expected, attribute_type in [
                ('ID', 'id0', str),
                ('Dbxref', 'taxon:10090', str),
                ('Name', 1, int),
                ('chromosome1', 1, int),
                ('gbkey', 'Src', str),
                ('genome', 'chromosome', str),
                ('mol_type', 'genomic DNA', str),
                ]:

            self.assertIn(name, p.terms)
            self.assertEqual(p.terms[name], {line_no: expected})
            self.assertIsInstance(p.terms[name][line_no], attribute_type)

    def test_ignore_space(self):
        p = AttributesParser(' ', ignore=['chromosome', 'junk'])
        record = 'chromosome 1;ID id0;junk drawer'

        t = list(p.tokenize(record))
        self.assertEqual(t, [
            'chromosome', ' ', '1', ';',
            'ID', ' ', 'id0', ';',
            'junk', ' ', 'drawer'
        ])
        attributes = p(record)
        self.assertEqual(attributes, 1)
        self.assertIn('ID', p.terms)
        self.assertEqual(p.terms['ID'][0], 'id0')
        self.assertNotIn('chromosome', p.terms)
        self.assertEqual(p.reserved['chromosome'], 0)

    def test_ignore_equal(self):
        p = AttributesParser('=', ignore=['chromosome', 'junk'])
        record = 'chromosome=1;ID=id0;junk=drawer'

        t = list(p.tokenize(record))
        self.assertEqual(t, [
            'chromosome', '=', '1', ';',
            'ID', '=', 'id0', ';',
            'junk', '=', 'drawer'
        ])
        attributes = p(record)
        self.assertEqual(attributes, 1)
        self.assertIn('ID', p.terms)
        self.assertEqual(p.terms['ID'][0], 'id0')
        self.assertNotIn('chromosome', p.terms)
        self.assertEqual(p.reserved['chromosome'], 0)


class TestParseHelpers(unittest.TestCase):
    def test_parse_score(self):
        self.assertIs(parse_score('.'), numpy.nan)
        self.assertEqual(parse_score('1000'), 1000)
        self.assertEqual(parse_score('0'), 0)
        self.assertEqual(parse_score('0.5'), 0.5)

    def test_parse_strand(self):
        self.assertEqual(parse_strand('+'), 1)
        self.assertIs(parse_strand('.'), numpy.nan)
        self.assertEqual(parse_strand('-'), -1)
        self.assertIs(parse_strand('junk'), numpy.nan)

    def test_parse_phase(self):
        # should we test for
        self.assertIs(parse_phase('.'), numpy.nan)
        for i in ['0', '1', '2']:
            self.assertEqual(parse_phase(i), int(i))
        self.assertRaises(ValueError, parse_phase, '5')


class TestGFFParser(unittest.TestCase):
    def test_parse_simple_gtf(self):
        text = StringIO('''chr1	source	exon	1	1000	.	+	.	junk "drawer"''')
        p = GFFParser(' ')
        p.read_gff(text)
        self.assertEqual(p.gtf.shape, (1, 9))
        self.assertIn('chromosome', p.gtf.columns)
        self.assertIn('junk', p.gtf.columns)

    def test_parse_collision_gtf(self):
        """Make sure that the mandatory first 8 columns keep the correct name
        """
        text = StringIO('''chr1	source	exon	1	1000	.	+	.	chromosome "drawer"''')
        p = GFFParser(' ')
        p.read_gff(text)
        self.assertEqual(p.gtf.shape, (1, 9))
        self.assertEqual('chromosome', p.gtf.columns[0])
        self.assertEqual('chromosome1', p.gtf.columns[8])


