import pandas
from unittest import TestCase

from .samsimulator import SAMGenerator
from ..compute_read_distribution import (
    build_gene_locations,
    count_exonic_genomic_reads_for_reference,
    ReferenceDistribution,
)


# def count_exonic_genomic_reads(bam, gtf_name, stranded="unstranded"):
# def count_exonic_genomic_reads_for_reference(


class TestBuildGeneLocations(TestCase):
    def setUp(self):
        # assume input GTF file is 1-based
        self.gtf = pandas.DataFrame(
            {
                "chromosome": ["chr1"] * 3,
                "type": ["gene", "transcript", "exon"],
                "start": [100, 100, 105],
                "stop": [200, 200, 195],
                "strand": [1, 1, 1],
                "gene_id": ["Gene-1", "Gene-1", "Gene-1"],
            }
        )

    def test_build_gene_locations_simple_forward(self):
        gene_plus, gene_minus = build_gene_locations(self.gtf, "test")
        print(gene_plus)
        assert len(gene_minus) == 0
        assert list(range(99, 200)) == list(gene_plus)
        for i in range(99, 200):
            if i >= 104 and i < 195:
                assert gene_plus[i] == "E", "{} {} != E".format(i, gene_plus[i])
            else:
                assert gene_plus[i] == "G", "{} {} != G".format(i, gene_plus[i])

    def test_build_gene_locations_simple_reverse(self):
        self.gtf["strand"] = [-1] * len(self.gtf["start"])
        gene_plus, gene_minus = build_gene_locations(self.gtf, "test")
        assert len(gene_plus) == 0
        assert list(range(99, 200)) == list(gene_minus)
        for i in range(100, 200):
            if i >= 104 and i < 195:
                assert gene_minus[i] == "E", "{} {} != E".format(i, gene_minus[i])
            else:
                assert gene_minus[i] == "G", "{} {} != G".format(i, gene_minus[i])

    def test_build_gene_locations_simple_reverse_unstranded(self):
        self.gtf["strand"] = [-1] * len(self.gtf["start"])
        gene_plus, gene_minus = build_gene_locations(
            self.gtf, "test", check_strand=False
        )
        assert len(gene_minus) == 0
        assert list(range(99, 200)) == list(gene_plus)
        for i in range(100, 200):
            if i >= 104 and i < 195:
                assert gene_plus[i] == "E", "{} {} != E".format(i, gene_plus[i])
            else:
                assert gene_plus[i] == "G", "{} {} != G".format(i, gene_plus[i])


class TestCountExonicGenes(TestCase):
    def setUp(self):
        # Gtf files are one based
        start = [101, 101, 106, 151, 181]
        stop = [201, 201, 126, 176, 196]
        gtf = pandas.DataFrame(
            {
                "chromosome": ["chr1"] * len(start),
                "type": ["gene", "transcript", "exon", "exon", "exon"],
                "start": start,
                "stop": stop,
                "strand": [1] * len(start),
                "gene_id": ["Gene-1"] * len(start),
            }
        )

        self.gene_plus, self.gene_minus = build_gene_locations(gtf, "chr1")

    def test_one_exon_forward_read(self):
        sam = SAMGenerator()
        chromosome = "chr1"
        sam.add_single_read0("r0001", chromosome, 105, seq="A" * 25)

        expected_data = [
            ("forward", ReferenceDistribution(1.0, 0, 0)),
            ("unstranded", ReferenceDistribution(1.0, 0, 0)),
            ("reverse", ReferenceDistribution(0, 0, 1.0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)

    def test_one_intron_forward_read(self):
        sam = SAMGenerator()
        sam.add_single_read0("r0001", "chr1", 125, seq="A" * 20)

        expected_data = [
            ("forward", ReferenceDistribution(0, 1.0, 0)),
            ("unstranded", ReferenceDistribution(0, 1.0, 0)),
            ("reverse", ReferenceDistribution(0, 0, 1.0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)

    def test_one_genome_forward_read(self):
        sam = SAMGenerator()
        sam.add_single_read0("r0001", "chr1", 0, seq="A" * 50)

        expected_data = [
            ("forward", ReferenceDistribution(0, 0, 1.0)),
            ("unstranded", ReferenceDistribution(0, 0, 1.0)),
            ("reverse", ReferenceDistribution(0, 0, 1.0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)

    def test_one_exon_reverse_read(self):
        sam = SAMGenerator()
        chromosome = "chr1"
        sam.add_single_read0("r0001", chromosome, 105, seq="A" * 25, is_reverse=True)

        expected_data = [
            ("forward", ReferenceDistribution(0, 0, 1.0)),
            ("unstranded", ReferenceDistribution(1.0, 0, 0)),
            ("reverse", ReferenceDistribution(1.0, 0, 0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)

    def test_one_intron_reverse_read(self):
        sam = SAMGenerator()
        sam.add_single_read0("r0001", "chr1", 125, seq="A" * 20, is_reverse=True)

        expected_data = [
            ("forward", ReferenceDistribution(0, 0, 1.0)),
            ("unstranded", ReferenceDistribution(0, 1.0, 0)),
            ("reverse", ReferenceDistribution(0, 1.0, 0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)

    def test_one_genome_reverse_read(self):
        sam = SAMGenerator()
        sam.add_single_read0("r0001", "chr1", 0, seq="A" * 50, is_reverse=True)

        expected_data = [
            ("forward", (0, 0, 1.0)),
            ("unstranded", (0, 0, 1.0)),
            ("reverse", (0, 0, 1.0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)

    def test_one_partial_exon_forward_read(self):
        sam = SAMGenerator()
        chromosome = "chr1"
        sam.add_single_read0("r0001", chromosome, 101, seq="A" * 10)

        expected_data = [
            ("forward", (1.0, 0, 0)),
            ("unstranded", (1.0, 0, 0)),
            ("reverse", (0, 0, 1.0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)

    def test_one_partial_intron_on_exon_forward_read(self):
        sam = SAMGenerator()
        chromosome = "chr1"
        sam.add_single_read0("r0001", chromosome, 99, seq="A" * 10)

        expected_data = [
            ("forward", (0, 1.0, 0)),
            ("unstranded", (0, 1.0, 0)),
            ("reverse", (0, 0, 1.0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)

    def test_one_partial_gene_forward_read(self):
        sam = SAMGenerator()
        chromosome = "chr1"
        # this is 6 bases over gene, 4 bases over exon
        # 3 _, 5 G, 1 E
        sam.add_single_read0("r0001", chromosome, 95, seq="A" * 10)

        expected_data = [
            ("forward", (0, 1.0, 0)),
            ("unstranded", (0, 1.0, 0)),
            ("reverse", (0, 0, 1.0)),
        ]

        for stranded, expected in expected_data:
            with sam.to_bam() as bam:
                counts = count_exonic_genomic_reads_for_reference(
                    bam, "chr1", self.gene_plus, self.gene_minus, stranded
                )

                self.assertEqual(counts, expected)
