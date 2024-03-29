from unittest import TestCase

from .samsimulator import (
    tokenize_cigar,
    make_sam_flags,
    SAMGenerator,
)


def create_reads():
    sam = SAMGenerator()
    sam.add_paired_read0("r0001", "chr1", [0, 20], seq=["A" * 5, "T" * 10]),
    sam.add_paired_read0("r0002", "chr1", [5, 26], seq=["A" * 10, "C" * 10]),
    sam.add_paired_read0("r0003", "chr1", [10, 31], seq=["A" * 10, "C" * 10]),
    sam.add_paired_read0(
        "r0004", "chr1", [10, 31], seq=["A" * 10, "C" * 10], is_reverse=[True, True]
    ),
    sam.add_paired_read0(
        "r0005",
        "chr1",
        [451, 551],
        seq=["A" * 60, "A" * 60],
        is_reverse=[True, False],
    ),
    sam.add_paired_read0("r0006", "chr1", [500, 601], seq=["A" * 50, "T" * 50]),
    sam.add_paired_read0("r0007", "chr1", [740, 841], seq=["A" * 70, "A" * 70]),
    sam.add_paired_read0("r0008", "chr1", [840, 941], seq=["A" * 10, "A" * 10]),
    sam.add_paired_read0(
        "r0009",
        "chr1",
        [901, 1201],
        seq=["A" * 200, "A" * 200],
        cigar=["100M100N100M", "150M50S"],
    ),
    sam.add_paired_read0("r0010", "chr1", [1500, 1601], seq=["A" * 50, "A" * 50]),
    sam.add_paired_read0("r0011", "chr1", [2000, 2101], seq=["A" * 50, "A" * 50]),
    sam.add_paired_read0("r0012", "chr1", [3000, 3201], seq=["A" * 100, "T" * 100]),
    sam.add_paired_read0(
        "r0013", "*", [0, 0], flag=[4, 4], seq=["A" * 50, "A" * 50]
    ),

    return sam


class TestTokenizeCigar(TestCase):
    def test_tokenize(self):
        self.assertEqual(tokenize_cigar("100M"), [(100, "M")])
        self.assertEqual(tokenize_cigar("10M20S10M"), [(10, "M"), (20, "S"), (10, "M")])


class TestSamFlags(TestCase):
    def test_make_sam_flags(self):
        self.assertEqual(make_sam_flags(), 2)
        self.assertEqual(make_sam_flags(unmapped=True), 6)
        self.assertEqual(make_sam_flags(aligned=False, unmapped=True), 4)


class TestSamGenerator(TestCase):
    def test_make_single_sam(self):
        sam = SAMGenerator()
        sam.add_single_read0("r0001", "chr1", 0, seq="A" * 5),
        sam.add_single_read0("r0002", "chr1", 5, seq="A" * 10),
        sam.add_single_read0("r0003", "chr1", 10, seq="A" * 10),
        sam.add_single_read0("r0004", "chr1", 10, seq="A" * 10, is_reverse=True),
        sam.add_single_read0("r0005", "chr1", 450, seq="A" * 60, is_reverse=True),
        sam.add_single_read0("r0006", "chr1", 500, seq="A" * 50),
        sam.add_single_read0("r0007", "chr1", 740, seq="A" * 70),
        sam.add_single_read0("r0008", "chr1", 840, seq="A" * 10),
        sam.add_single_read0("r0009", "chr1", 900, seq="A" * 200, cigar="100M100N100M"),
        sam.add_single_read0("r0010", "chr1", 1500, seq="A" * 50),
        sam.add_single_read0("r0011", "chr1", 2000, seq="A" * 50),
        sam.add_single_read0("r0012", "chr1", 3000, seq="A" * 100),
        sam.add_single_read0("r0013", "*", 0, flag=4, seq="A" * 50),

        with sam.to_sam() as sam:
            # self.assertEqual(bam.references, ("chr1",))
            # self.assertEqual(bam.lengths,  (3101,))

            for read in sam:
                if read.query_name == "r001":
                    self.assertEqual(read.get_reference_position(), 1)
                    self.assertEqual(read.is_reverse, False)
                elif read.query_name == "r0004":
                    self.assertTrue(read.is_reverse)
                elif read.query_name == "r009":
                    self.assertEqual(read.cigar, [(100, "M"), (100, "N"), (100, "M")])
                elif read.query_name == "r0013":
                    self.assertEqual(read.is_unmapped, True)

    def test_make_single_sam(self):
        sam = create_reads()
        with sam.to_sam() as stream:
            # self.assertEqual(bam.references, ("chr1",))
            # self.assertEqual(bam.lengths,  (3101,))

            for read in stream:
                if read.query_name == "r001":
                    self.assertEqual(read.get_reference_position(), 1)
                    self.assertEqual(read.is_reverse, False)
                elif read.query_name == "r0004":
                    self.assertTrue(read.is_reverse)
                elif read.query_name == "r009":
                    self.assertEqual(read.cigar, [(100, "M"), (100, "N"), (100, "M")])
                elif read.query_name == "r0013":
                    self.assertEqual(read.is_unmapped, True)

    def test_make_single_bam(self):
        sam = create_reads()
        with sam.to_bam() as stream:
            # self.assertEqual(bam.references, ("chr1",))
            # self.assertEqual(bam.lengths,  (3101,))

            for read in stream:
                if read.query_name == "r001":
                    self.assertEqual(read.get_reference_position(), 1)
                    self.assertEqual(read.is_reverse, False)
                elif read.query_name == "r0004":
                    self.assertTrue(read.is_reverse)
                elif read.query_name == "r009":
                    self.assertEqual(read.cigar, [(100, "M"), (100, "N"), (100, "M")])
                elif read.query_name == "r0013":
                    self.assertEqual(read.is_unmapped, True)
                    
