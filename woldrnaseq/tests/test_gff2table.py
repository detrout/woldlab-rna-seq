from unittest import TestCase, main
from pkg_resources import resource_filename
from io import StringIO

from woldrnaseq.gff2table import (
    AttributesParser,
    GFFParser,
)

class TestAttributeParser(TestCase):
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
                ('chromosome', 1, int),
                ('gbkey', 'Src', str),
                ('genome', 'chromosome', str),
                ('mol_type', 'genomic DNA', str),
                ]:
            
            self.assertIn(name, p.terms)
            self.assertEqual(p.terms[name], {line_no: expected})
            self.assertIsInstance(p.terms[name][line_no], attribute_type)

    def test_ignore(self):
        p = AttributesParser(' ', ignore=['chromosome', 'junk'])
        attributes = p('chromsome=1;ID=id0;junk=drawer')
        self.assertEqual(attributes, 1)
        self.assertIn('ID', p.terms)
        self.assertEqual(p.terms['ID'][0], 'id0')
        self.assertNotIn('chromosome', p.terms)
        self.assertEqual(p.reserved['chromosome'], 0)


class TestAttributeParser(TestCase):
    def test_parse_simple_gtf(self):
        text = StringIO('''chr1	source	exon	1	1000	.	+	.	junk "drawer"''')
        p = GFFParser(' ')
        p.read_gff(text)
        self.assertEqual(p.gtf.shape, (1, 9))
        self.assertIn('chromosome', p.gtf.columns)
        self.assertIn('junk', p.gtf.columns)

    def test_parse_collision_gtf(self):
        text = StringIO('''chr1	source	exon	1	1000	.	+	.	chromosome "drawer"''')
        p = GFFParser(' ')
        p.read_gff(text)
        self.assertEqual(p.gtf.shape, (1, 9))
        self.assertIn('chromosome', p.gtf.columns)
        self.assertIn('chromosome1', p.gtf.columns)

