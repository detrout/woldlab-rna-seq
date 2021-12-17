from unittest import TestCase

from woldrnaseq.downloader import (
    make_short_fastq_name
)


class TestDownloader(TestCase):
    def test_make_short_fastq_name(self):
        self.assertEqual(
            make_short_fastq_name('22160_GAACGAAG_L002_R1_005.fastq.gz'),
            '22160_GAACGAAG_L002_R1.fastq.gz')

        self.assertRaises(
            ValueError,
            make_short_fastq_name,
            '22160_GAACGAAG_R1.fastq.gz')
