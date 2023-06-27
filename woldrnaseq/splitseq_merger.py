"""Read a STAR Solo.out matrix and merge barcodes into logical cells

Split-seq is different from many of the other UMI tagging methods as
they use two different barcodes to represent a single cell, one
barcode is optimized to detect polyA messages, while the other is
random primed.

However that means to generate comparable data we need to merge the
two related barcodes into a single logical cell.
"""
import csv
import os
from pathlib import Path
import shutil
import scipy
from scipy.io.mmio import MMFile


splitseq_decoder = {
    "v1": {
        'ACTCGTAA': '_0', 'CTGCTTTG': '_0',
        'AAACGATA': '_1', 'CATGATCA': '_1',
        'TTACCTCG': '_2', 'GGGTAGCG': '_2',
        'GCCTGCAA': '_3', 'CCGAGAAA': '_3',
        'TGGTATAC': '_4', 'ACGGACTC': '_4',
        'CGTTCGAG': '_5', 'ACTTACGA': '_5',
        'TCTATTAC': '_6', 'TATTTAAG': '_6',
        'ATAAGCTC': '_7', 'ACCGTACG': '_7',
        'ATTCATGG': '_8', 'TATAGTCG': '_8',
        'ATCCGCGA': '_9', 'TGGGCATC': '_9',
        'ATCGCATA': '_10', 'TACCTAGA': '_10',
        'CCGTTCTA': '_11', 'GCTGCATG': '_11',
        'TGGCGCGC': '_12', 'GTCATATG': '_12',
        'TGTCTGAA': '_13', 'ATATTGGC': '_13',
        'CTGTCCCG': '_14', 'CTAAGGGA': '_14',
        'AATTTCTC': '_15', 'TCGTTTCG': '_15',
        'CGCGACTA': '_16', 'GAATAATG': '_16',
        'TGAAGCAA': '_17', 'ACTGCGCA': '_17',
        'TTATTCTG': '_18', 'GCTTATAG': '_18',
        'GTTCAACA': '_19', 'ATCATGCA': '_19',
        'ACGCCGGC': '_20', 'ACGTTAAC': '_20',
        'TTGTCTTA': '_21', 'CCATCTTG': '_21',
        'TACGGTTA': '_22', 'CATAGCTA': '_22',
        'TTGGGAGA': '_23', 'GAGGTTGA': '_23',
        'TGCTTGGG': '_24', 'GCACTGAC': '_24',
        'TAAATATC': '_25', 'TTCATCGC': '_25',
        'CACAATTG': '_26', 'GAAATTAG': '_26',
        'GTGCTAGC': '_27', 'AGGATTAA': '_27',
        'CGCCCGGA': '_28', 'AATAGAAC': '_28',
        'GCTCGCGG': '_29', 'TCTTAATC': '_29',
        'CTTTGGTC': '_30', 'TAATACGC': '_30',
        'TTCCGATC': '_31', 'GTTTGTGA': '_31',
        'TTCGCTAC': '_32', 'CGAACGTC': '_32',
        'AGCGAAAC': '_33', 'GGTTCTTC': '_33',
        'AAATAGCA': '_34', 'GCAAATTC': '_34',
        'CGTCTAGG': '_35', 'GCTATGCG': '_35',
        'GCCGTGTA': '_36', 'CTACCCTA': '_36',
        'CGCTTAAA': '_37', 'GTGGGTTC': '_37',
        'GACCTTTC': '_38', 'GTCCGTAG': '_38',
        'GGTGGAGC': '_39', 'TGCGATCG': '_39',
        'TACTCGAA': '_40', 'TATCCGGG': '_40',
        'CATTTGGA': '_41', 'AGGTAATA': '_41',
        'GAGCACAA': '_42', 'CGTGGTTG': '_42',
        'GTCGCGCG': '_43', 'GACAAAGC': '_43',
        'GTTACGTA': '_44', 'GGGCGATG': '_44',
        'CTATTTCA': '_45', 'ATCTATAA': '_45',
        'ACTATATA': '_46', 'GCCCATGA': '_46',
        'TCACTTTA': '_47', 'CTGAAAGG': '_47'
    },
    "v2": {
        'CATTCCTA': '_0', 'CATCATCC': '_0',
        'CTTCATCA': '_1', 'CTGCTTTG': '_1',
        'CCTATATC': '_2', 'CTAAGGGA': '_2',
        'ACATTTAC': '_3', 'GCTTATAG': '_3',
        'ACTTAGCT': '_4', 'TCTGATCC': '_4',
        'CCAATTCT': '_5', 'TCTCTTGG': '_5',
        'GCCTATCT': '_6', 'CAATTTCC': '_6',
        'ATGCTGCT': '_7', 'AGTCTCTT': '_7',
        'CATTTACA': '_8', 'TGCTGCTC': '_8',
        'ACTCGTAA': '_9', 'GTATTTCC': '_9',
        'CCTTTGCA': '_10', 'TTCCTGTG': '_10',
        'ACTCCTGC': '_11', 'GCTGCTTC': '_11',
        'ATTTGGCA': '_12', 'TATGTGTC': '_12',
        'TTATTCTG': '_13', 'CAATTCTC': '_13',
        'TCATGCTC': '_14', 'TGGTCTCC': '_14',
        'CATACTTC': '_15', 'GCTCTTTA': '_15',
        'CCGTTCTA': '_16', 'GCTGCATG': '_16',
        'GCTTCATA': '_17', 'ACTCATTT': '_17',
        'CTCTGTGC': '_18', 'AGTCTTGG': '_18',
        'CCCTTATA': '_19', 'GGTTCTTC': '_19',
        'ACTGCTCT': '_20', 'TCATGTTG': '_20',
        'CTCTAATC': '_21', 'ATTTTGCC': '_21',
        'ACCCTTGC': '_22', 'CTTCTGTA': '_22',
        'ATCTTAGG': '_23', 'GTCCATCT': '_23',
        'CATGTCTC': '_24', 'GCTATCTC': '_24',
        'TCATTGCA': '_25', 'TAGTTTCC': '_25',
        'ACACCTTT': '_26', 'TCCATTAT': '_26',
        'AATTTCTC': '_27', 'AGGATTAA': '_27',
        'ATTCATGG': '_28', 'AATCTTTC': '_28',
        'ACTTTACC': '_29', 'GTCATATG': '_29',
        'CTTCTAAC': '_30', 'GTGCTTCC': '_30',
        'CTATTTCA': '_31', 'ATGTGTTG': '_31',
        'TCTCATGC': '_32', 'CCATCTTG': '_32',
        'ATCCTTAC': '_33', 'TACTGTCT': '_33',
        'TAAATATC': '_34', 'TTCATCGC': '_34',
        'TTACCTGC': '_35', 'ACTGTGGG': '_35',
        'CACTTTCA': '_36', 'TCTGTGCC': '_36',
        'CACCTTTA': '_37', 'TCAATCTC': '_37',
        'CTGACTTC': '_38', 'GTCCTCTG': '_38',
        'CATTTGGA': '_39', 'TTACATTC': '_39',
        'GCTCTACT': '_40', 'ATTCTGTC': '_40',
        'GTTACGTA': '_41', 'TGTGTATG': '_41',
        'CCTGTTGC': '_42', 'TCCATTTG': '_42',
        'CTATCATC': '_43', 'TTAGCTTC': '_43',
        'GCTATCAT': '_44', 'GTGCTTGA': '_44',
        'ACATTCAT': '_45', 'GTTTGTGA': '_45',
        'TTCGCTAC': '_46', 'GAAATTAG': '_46',
        'CATTCTAC': '_47', 'GCAAATTC': '_47',
        'CACTTATC': '_48', 'GAGGTTGA': '_48',
        'ATAAGCTC': '_49', 'CCTGTCTG': '_49',
        'TCATCCTG': '_50', 'GTGGGTTC': '_50',
        'CCTGGTAT': '_51', 'TTTGCATC': '_51',
        'TGGTATAC': '_52', 'AGGTAATA': '_52',
        'TTGGGAGA': '_53', 'GTGCCTTC': '_53',
        'ACTTCATC': '_54', 'ATGTTTCC': '_54',
        'TCTCTAGC': '_55', 'CTTAATTC': '_55',
        'ATGCCCTT': '_56', 'TCTGGCTC': '_56',
        'CCCAATTT': '_57', 'CATCATTT': '_57',
        'ACTATATA': '_58', 'GTTGTCTC': '_58',
        'CTCTATAC': '_59', 'ATCTTCTG': '_59',
        'CTGTCTCA': '_60', 'TGTTTGCC': '_60',
        'GACCTTTC': '_61', 'TTCTGTCA': '_61',
        'GATTTGGC': '_62', 'ACGGACTC': '_62',
        'CGTCTAGG': '_63', 'TTTGGTCA': '_63',
        'TACTCGAA': '_64', 'TATCCGGG': '_64',
        'CAGCCTTT': '_65', 'TGTCATTC': '_65',
        'CCTCATTA': '_66', 'ATTCTCTG': '_66',
        'CTTATACC': '_67', 'TGGCTTCC': '_67',
        'TCTATTAC': '_68', 'TTGTTGCC': '_68',
        'CCTGCATT': '_69', 'GTCATCTC': '_69',
        'CAATCCTT': '_70', 'TTGCTCAT': '_70',
        'TTGTCTTA': '_71', 'CTGTCTGC': '_71',
        'TCACTTTA': '_72', 'TATATTCC': '_72',
        'TGCTTGGG': '_73', 'ATATTGGC': '_73',
        'CGCTCATT': '_74', 'GTGTCCTC': '_74',
        'GCCTCTAT': '_75', 'ATCTTCAT': '_75',
        'GAGCACAA': '_76', 'CGTGGTTG': '_76',
        'CTCTTAAC': '_77', 'TTGCATCC': '_77',
        'TCTAGGCT': '_78', 'TCTTAATC': '_78',
        'AATTCTGC': '_79', 'TGCATTTC': '_79',
        'CATTCTCA': '_80', 'GATGTTTC': '_80',
        'ACTTGCCT': '_81', 'ATCTTGTC': '_81',
        'ATCATTGC': '_82', 'TCATATTC': '_82',
        'GTTCAACA': '_83', 'TGGCCTCT': '_83',
        'CCATTTGC': '_84', 'CGTTGTCT': '_84',
        'GACTTTGC': '_85', 'TCTTGTCA': '_85',
        'ATTGGCTC': '_86', 'TATTCCTG': '_86',
        'GTGCTAGC': '_87', 'TCCATGTC': '_87',
        'CTTTCAAC': '_88', 'TTGTCATC': '_88',
        'ACTATTGC': '_89', 'ATTTCCTG': '_89',
        'ACTGGCTT': '_90', 'GTGTCTCC': '_90',
        'ATTAGGCT': '_91', 'GTGTGTGT': '_91',
        'GCCTTTCA': '_92', 'TATGCTTC': '_92',
        'ATTCTAGG': '_93', 'ATGGTGTT': '_93',
        'CCTTACAT': '_94', 'GAATAATG': '_94',
        'ACATTTGG': '_95', 'CCTCTGTG': '_95',
    },
    "v2-t": {
        'CATTCCTA': 'CATTCCTA', 'CATCATCC': 'CATTCCTA',
        'CTTCATCA': 'CTTCATCA', 'CTGCTTTG': 'CTTCATCA',
        'CCTATATC': 'CCTATATC', 'CTAAGGGA': 'CCTATATC',
        'ACATTTAC': 'ACATTTAC', 'GCTTATAG': 'ACATTTAC',
        'ACTTAGCT': 'ACTTAGCT', 'TCTGATCC': 'ACTTAGCT',
        'CCAATTCT': 'CCAATTCT', 'TCTCTTGG': 'CCAATTCT',
        'GCCTATCT': 'GCCTATCT', 'CAATTTCC': 'GCCTATCT',
        'ATGCTGCT': 'ATGCTGCT', 'AGTCTCTT': 'ATGCTGCT',
        'CATTTACA': 'CATTTACA', 'TGCTGCTC': 'CATTTACA',
        'ACTCGTAA': 'ACTCGTAA', 'GTATTTCC': 'ACTCGTAA',
        'CCTTTGCA': 'CCTTTGCA', 'TTCCTGTG': 'CCTTTGCA',
        'ACTCCTGC': 'ACTCCTGC', 'GCTGCTTC': 'ACTCCTGC',
        'ATTTGGCA': 'ATTTGGCA', 'TATGTGTC': 'ATTTGGCA',
        'TTATTCTG': 'TTATTCTG', 'CAATTCTC': 'TTATTCTG',
        'TCATGCTC': 'TCATGCTC', 'TGGTCTCC': 'TCATGCTC',
        'CATACTTC': 'CATACTTC', 'GCTCTTTA': 'CATACTTC',
        'CCGTTCTA': 'CCGTTCTA', 'GCTGCATG': 'CCGTTCTA',
        'GCTTCATA': 'GCTTCATA', 'ACTCATTT': 'GCTTCATA',
        'CTCTGTGC': 'CTCTGTGC', 'AGTCTTGG': 'CTCTGTGC',
        'CCCTTATA': 'CCCTTATA', 'GGTTCTTC': 'CCCTTATA',
        'ACTGCTCT': 'ACTGCTCT', 'TCATGTTG': 'ACTGCTCT',
        'CTCTAATC': 'CTCTAATC', 'ATTTTGCC': 'CTCTAATC',
        'ACCCTTGC': 'ACCCTTGC', 'CTTCTGTA': 'ACCCTTGC',
        'ATCTTAGG': 'ATCTTAGG', 'GTCCATCT': 'ATCTTAGG',
        'CATGTCTC': 'CATGTCTC', 'GCTATCTC': 'CATGTCTC',
        'TCATTGCA': 'TCATTGCA', 'TAGTTTCC': 'TCATTGCA',
        'ACACCTTT': 'ACACCTTT', 'TCCATTAT': 'ACACCTTT',
        'AATTTCTC': 'AATTTCTC', 'AGGATTAA': 'AATTTCTC',
        'ATTCATGG': 'ATTCATGG', 'AATCTTTC': 'ATTCATGG',
        'ACTTTACC': 'ACTTTACC', 'GTCATATG': 'ACTTTACC',
        'CTTCTAAC': 'CTTCTAAC', 'GTGCTTCC': 'CTTCTAAC',
        'CTATTTCA': 'CTATTTCA', 'ATGTGTTG': 'CTATTTCA',
        'TCTCATGC': 'TCTCATGC', 'CCATCTTG': 'TCTCATGC',
        'ATCCTTAC': 'ATCCTTAC', 'TACTGTCT': 'ATCCTTAC',
        'TAAATATC': 'TAAATATC', 'TTCATCGC': 'TAAATATC',
        'TTACCTGC': 'TTACCTGC', 'ACTGTGGG': 'TTACCTGC',
        'CACTTTCA': 'CACTTTCA', 'TCTGTGCC': 'CACTTTCA',
        'CACCTTTA': 'CACCTTTA', 'TCAATCTC': 'CACCTTTA',
        'CTGACTTC': 'CTGACTTC', 'GTCCTCTG': 'CTGACTTC',
        'CATTTGGA': 'CATTTGGA', 'TTACATTC': 'CATTTGGA',
        'GCTCTACT': 'GCTCTACT', 'ATTCTGTC': 'GCTCTACT',
        'GTTACGTA': 'GTTACGTA', 'TGTGTATG': 'GTTACGTA',
        'CCTGTTGC': 'CCTGTTGC', 'TCCATTTG': 'CCTGTTGC',
        'CTATCATC': 'CTATCATC', 'TTAGCTTC': 'CTATCATC',
        'GCTATCAT': 'GCTATCAT', 'GTGCTTGA': 'GCTATCAT',
        'ACATTCAT': 'ACATTCAT', 'GTTTGTGA': 'ACATTCAT',
        'TTCGCTAC': 'TTCGCTAC', 'GAAATTAG': 'TTCGCTAC',
        'CATTCTAC': 'CATTCTAC', 'GCAAATTC': 'CATTCTAC',
        'CACTTATC': 'CACTTATC', 'GAGGTTGA': 'CACTTATC',
        'ATAAGCTC': 'ATAAGCTC', 'CCTGTCTG': 'ATAAGCTC',
        'TCATCCTG': 'TCATCCTG', 'GTGGGTTC': 'TCATCCTG',
        'CCTGGTAT': 'CCTGGTAT', 'TTTGCATC': 'CCTGGTAT',
        'TGGTATAC': 'TGGTATAC', 'AGGTAATA': 'TGGTATAC',
        'TTGGGAGA': 'TTGGGAGA', 'GTGCCTTC': 'TTGGGAGA',
        'ACTTCATC': 'ACTTCATC', 'ATGTTTCC': 'ACTTCATC',
        'TCTCTAGC': 'TCTCTAGC', 'CTTAATTC': 'TCTCTAGC',
        'ATGCCCTT': 'ATGCCCTT', 'TCTGGCTC': 'ATGCCCTT',
        'CCCAATTT': 'CCCAATTT', 'CATCATTT': 'CCCAATTT',
        'ACTATATA': 'ACTATATA', 'GTTGTCTC': 'ACTATATA',
        'CTCTATAC': 'CTCTATAC', 'ATCTTCTG': 'CTCTATAC',
        'CTGTCTCA': 'CTGTCTCA', 'TGTTTGCC': 'CTGTCTCA',
        'GACCTTTC': 'GACCTTTC', 'TTCTGTCA': 'GACCTTTC',
        'GATTTGGC': 'GATTTGGC', 'ACGGACTC': 'GATTTGGC',
        'CGTCTAGG': 'CGTCTAGG', 'TTTGGTCA': 'CGTCTAGG',
        'TACTCGAA': 'TACTCGAA', 'TATCCGGG': 'TACTCGAA',
        'CAGCCTTT': 'CAGCCTTT', 'TGTCATTC': 'CAGCCTTT',
        'CCTCATTA': 'CCTCATTA', 'ATTCTCTG': 'CCTCATTA',
        'CTTATACC': 'CTTATACC', 'TGGCTTCC': 'CTTATACC',
        'TCTATTAC': 'TCTATTAC', 'TTGTTGCC': 'TCTATTAC',
        'CCTGCATT': 'CCTGCATT', 'GTCATCTC': 'CCTGCATT',
        'CAATCCTT': 'CAATCCTT', 'TTGCTCAT': 'CAATCCTT',
        'TTGTCTTA': 'TTGTCTTA', 'CTGTCTGC': 'TTGTCTTA',
        'TCACTTTA': 'TCACTTTA', 'TATATTCC': 'TCACTTTA',
        'TGCTTGGG': 'TGCTTGGG', 'ATATTGGC': 'TGCTTGGG',
        'CGCTCATT': 'CGCTCATT', 'GTGTCCTC': 'CGCTCATT',
        'GCCTCTAT': 'GCCTCTAT', 'ATCTTCAT': 'GCCTCTAT',
        'GAGCACAA': 'GAGCACAA', 'CGTGGTTG': 'GAGCACAA',
        'CTCTTAAC': 'CTCTTAAC', 'TTGCATCC': 'CTCTTAAC',
        'TCTAGGCT': 'TCTAGGCT', 'TCTTAATC': 'TCTAGGCT',
        'AATTCTGC': 'AATTCTGC', 'TGCATTTC': 'AATTCTGC',
        'CATTCTCA': 'CATTCTCA', 'GATGTTTC': 'CATTCTCA',
        'ACTTGCCT': 'ACTTGCCT', 'ATCTTGTC': 'ACTTGCCT',
        'ATCATTGC': 'ATCATTGC', 'TCATATTC': 'ATCATTGC',
        'GTTCAACA': 'GTTCAACA', 'TGGCCTCT': 'GTTCAACA',
        'CCATTTGC': 'CCATTTGC', 'CGTTGTCT': 'CCATTTGC',
        'GACTTTGC': 'GACTTTGC', 'TCTTGTCA': 'GACTTTGC',
        'ATTGGCTC': 'ATTGGCTC', 'TATTCCTG': 'ATTGGCTC',
        'GTGCTAGC': 'GTGCTAGC', 'TCCATGTC': 'GTGCTAGC',
        'CTTTCAAC': 'CTTTCAAC', 'TTGTCATC': 'CTTTCAAC',
        'ACTATTGC': 'ACTATTGC', 'ATTTCCTG': 'ACTATTGC',
        'ACTGGCTT': 'ACTGGCTT', 'GTGTCTCC': 'ACTGGCTT',
        'ATTAGGCT': 'ATTAGGCT', 'GTGTGTGT': 'ATTAGGCT',
        'GCCTTTCA': 'GCCTTTCA', 'TATGCTTC': 'GCCTTTCA',
        'ATTCTAGG': 'ATTCTAGG', 'ATGGTGTT': 'ATTCTAGG',
        'CCTTACAT': 'CCTTACAT', 'GAATAATG': 'CCTTACAT',
        'ACATTTGG': 'ACATTTGG', 'CCTCTGTG': 'ACATTTGG',
    }
}


class SimpleMMWriter(MMFile):
    # The default scipy writer defaults to more zeros
    # than we need.
    @staticmethod
    def _field_template(field, precision):
        return {
            MMFile.FIELD_REAL: '%.{}g\n'.format(precision),
            MMFile.FIELD_INTEGER: '%i\n',
            MMFile.FIELD_UNSIGNED: '%u\n',
            MMFile.FIELD_COMPLEX: '%.{p}e %%.{p}e\n'.format(p=precision),
        }.get(field, None)


def read_parse_cell_barcode_lineno_map(stream):
    barcodes = {}
    reader = csv.reader(stream, delimiter="\t")
    for i, line in enumerate(reader):
        barcodes[line[0]] = i + 1

    return barcodes


def compute_parse_map(raw_barcodes, mapper):
    raw_to_collapsed_mapping = {}
    combined_indexes = {}
    for barcode in raw_barcodes:
        fragments = barcode.split("_")
        fragments[2] = str(mapper[fragments[2]])
        combined_barcode = "{}{}{}".format(*fragments)

        combined_index = combined_indexes.setdefault(combined_barcode, len(combined_indexes)+1)
        raw_index = raw_barcodes[barcode]
        raw_to_collapsed_mapping[raw_index] = combined_index
    return (combined_indexes, raw_to_collapsed_mapping)


def _parse_mmread(matrix_filename, merged_barcodes, merged_mapping):
    header = True
    matrix = None
    with open(matrix_filename, "rt") as instream:
        dtype = None
        for line in instream:
            if line.startswith("%%MatrixMarket"):
                mmid, matrix_word, fmt, field, symmetry = \
                    [part.strip() for part in line.split()]
                if not matrix_word.lower() == 'matrix':
                    raise ValueError("Problem reading file header")

                dtype = MMFile.DTYPES_BY_FIELD[field.lower()]
            elif line.startswith("%"):
                pass
            elif header:
                # After the comment comes the one header line
                total_features, total_cells, total_counts = [
                    int(x) for x in line.rstrip().split()
                ]
                matrix = scipy.sparse.dok_matrix((total_features, len(merged_barcodes)), dtype=dtype)
                header = False
            else:
                # row, column, count
                feature_index, cell_index, count = line.rstrip().split()
                feature_index = int(feature_index)
                cell_index = int(cell_index)
                count = float(count)

                new_cell_index = merged_mapping[cell_index]
                matrix[feature_index-1, new_cell_index-1] = matrix.get((feature_index-1, new_cell_index-1), 0) + count

    return matrix


def load_merged_parse_matrix(filename, splitseq_version="v1"):
    """Load STAR Solo matrix that needs to have split-seq cell barcodes merged

    It does assume that there will be a barcodes.tsv file next to the
    matrix file.
    """
    filename = Path(filename)
    cell_barcode_filename = filename.parent / "barcodes.tsv"

    with open(cell_barcode_filename, "rt") as instream:
        cell_barcodes = read_parse_cell_barcode_lineno_map(instream)
    merged_barcodes, merged_mapping = compute_parse_map(
        cell_barcodes, splitseq_decoder[splitseq_version])

    matrix = _parse_mmread(filename, merged_barcodes, merged_mapping)

    assert matrix.shape[1] == len(merged_barcodes), "The number of columns doesn't match the number of barcodes"
    return matrix, list(merged_barcodes.keys())


def write_merged_splitseq_matrix(matrix_filename, destination, library_id=None, splitseq_version="v1"):
    matrix_filename = Path(matrix_filename)
    destination = Path(destination)
    if not destination.exists():
        destination.mkdir()

    if not destination.is_dir():
        IOError("We can only write results to a directory {}".format(destination))

    source = matrix_filename.parent
    if not source.is_dir():
        IOError("Source {} needs to be a directory".format(source))

    matrix, cells = load_merged_parse_matrix(matrix_filename, splitseq_version=splitseq_version)

    with open(destination / "barcodes.tsv", "wt") as outstream:
        for barcode in cells:
            if library_id is not None:
                outstream.write(f"{library_id}-{barcode}\n")
            else:
                outstream.write(f"{barcode}\n")

    feature_filename = source / "features.tsv"
    shutil.copy(feature_filename, destination / "features.tsv")

    destination_matrix = str(destination / matrix_filename.name)
    SimpleMMWriter().write(
        destination_matrix, matrix, comment='', field=None, precision=None, symmetry=None)


def archive_split_seq_solo(
    solo_root,
    config,
    quantification="GeneFull",
    multiread="Unique",
    matrix="raw",
    *,
    destination=None,
    splitseq_version="v1",
):
    """Archive a STAR solo directory that contains a split-seq experiment

    In theory this shouldn't be needed as the current pipeline handle
    merging the two split-seq barcodes into one logical cell behind
    the scenes.

    But if one runs STAR by hand this might be useful to archive split
    seq result files.

    Parameters
    ----------
    solo_root : path to STAR's Solo.out directory where
        a file named {quantification}_{multiread}_{matrix}.tar.gz will be
        written.
    config : dictionary of configuration options
    quantification : Which counting method to use "Gene", "GeneFull",
        "GeneFull_Ex50pAS"
    multiread : which STAR EM processing level to use "Unique", "EM"
    matrix : which matrix to read either "raw" or "filtered"
    destination : what directory to write the archive to, defaults to
        solo_root/..

    """
    import gzip
    from io import BytesIO, StringIO
    import tarfile
    from mex_gene_archive.manifest import (
        compute_md5sums,
        create_metadata,
        write_manifest,
    )
    from mex_gene_archive.starsolo import (
        make_archive_root_name,
        make_list_of_archive_files,
        make_output_type_term,
        make_tar_archive_name,
        MULTIREAD_NAME,
        parse_star_log_out,
        update_tarinfo,
        validate_star_solo_out_arguments,
    )

    validate_star_solo_out_arguments(quantification, multiread, matrix)
    matrix_filename = solo_root / quantification / matrix / MULTIREAD_NAME[multiread]
    X, cells = load_merged_parse_matrix(matrix_filename)

    archive_files = {}
    for pathname in make_list_of_archive_files(solo_root, quantification, multiread, matrix):
        name = pathname.name
        if name == "barcodes.tsv":
            value = BytesIO(os.linesep.encode("ascii").join([x.encode("ascii") for x in cells]))
        elif name == "features.tsv":
            value = pathname
        else:
            value = BytesIO()
            SimpleMMWriter().write(value, X, comment='', field=None, precision=None, symmetry=None)
            value.seek(0)
        archive_files[pathname] = value

    config['output_type'] = make_output_type_term(quantification, multiread, matrix)
    config.update(parse_star_log_out(solo_root / ".." / "Log.out"))
    md5s = compute_md5sums(archive_files.values())
    manifest = create_metadata(config, md5s)
    manifest_buffer = BytesIO(
        write_manifest(StringIO(), manifest).getvalue().encode("utf-8")
    )
    archive_root = make_archive_root_name(solo_root, quantification, multiread, matrix)
    manifest_filename = (archive_root / "manifest.tsv").relative_to(solo_root)

    tar_name = make_tar_archive_name(
        solo_root, quantification, multiread, matrix, destination)
    with gzip.GzipFile(tar_name, "wb", mtime=0) as gzipstream:
        with tarfile.open(mode="w", fileobj=gzipstream, format=tarfile.PAX_FORMAT) as archive:
            info = tarfile.TarInfo(str(manifest_filename))
            update_tarinfo(info, fileobj=manifest_buffer)
            archive.addfile(info, manifest_buffer)
            for filename in archive_files:
                fileobj = archive_files[filename]
                info = tarfile.TarInfo(str(filename.relative_to(solo_root)))
                if isinstance(fileobj, str) or isinstance(fileobj, Path):
                    update_tarinfo(info, filename)
                    with open(filename, "rb") as instream:
                        archive.addfile(info, instream)
                else:
                    update_tarinfo(info, fileobj=fileobj)
                    archive.addfile(info, fileobj)

    return tar_name
