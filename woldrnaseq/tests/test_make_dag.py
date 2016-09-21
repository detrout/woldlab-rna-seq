from unittest import TestCase, main
from pkg_resources import resource_filename
from woldrnaseq import models, make_dag
import pandas

class TestMakeDag(TestCase):
    def test_reference_prefix(self):
        spurtsv = resource_filename(__name__, 'library-spur-se.tsv')
        spur = models.load_library_tables([spurtsv])

        self.assertEqual(make_dag.get_reference_prefix(spur, '12304'), 'scaffold')
        self.assertEqual(make_dag.get_reference_prefix(spur, '12307'), 'chr')

    def test_reference_prefix_missing(self):
        mm10tsv = resource_filename(__name__, 'library-mm10-se.tsv')
        mm10 = models.load_library_tables([mm10tsv])

        self.assertEqual(make_dag.get_reference_prefix(mm10, '12304'), 'chr')


if __name__ == '__main__':
    main()
