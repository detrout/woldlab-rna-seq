from unittest import TestCase, main
from pkg_resources import resource_filename

import numpy

from woldrnaseq.common import (
    get_seperator,
    normalize_path,
    get_fixed_range,
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

    def test_get_fixed_range(self):
        d = [(3,4), (1,2)]
        self.assertEqual(get_fixed_range(d), (1,4))

        d = [(4,3), (1,2)]
        self.assertRaises(AssertionError, get_fixed_range, d)

        d = [(3,4), (1,2), (numpy.nan, 5)]
        self.assertRaises(AssertionError, get_fixed_range, d)

if __name__ == '__main__':
    main()
