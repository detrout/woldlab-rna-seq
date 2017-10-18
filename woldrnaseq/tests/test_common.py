from unittest import TestCase, main
from pkg_resources import resource_filename

from woldrnaseq.common import (
    get_seperator,
    normalize_path,
)

class TestCommonUtilities(TestCase):
    def test_get_seperator(self):
        self.assertEqual(get_seperator('tab'), '\t')
        self.assertEqual(get_seperator('TaB'), '\t')
        self.assertEqual(get_seperator(','), ',')
        self.assertRaises(ValueError, get_seperator, '|')


    def test_normalize_path(self):
        self.assertEqual(normalize_path(None), None)
        self.assertEqual(normalize_path(''), '')
        self.assertEqual(normalize_path('/a/b'), '/a/b/')
        self.assertEqual(normalize_path('/a/b/'), '/a/b/')

if __name__ == '__main__':
    main()
