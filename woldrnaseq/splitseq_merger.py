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
        'ACTCGTAA': '0', 'CTGCTTTG': '0',
        'AAACGATA': '1', 'CATGATCA': '1',
        'TTACCTCG': '2', 'GGGTAGCG': '2',
        'GCCTGCAA': '3', 'CCGAGAAA': '3',
        'TGGTATAC': '4', 'ACGGACTC': '4',
        'CGTTCGAG': '5', 'ACTTACGA': '5',
        'TCTATTAC': '6', 'TATTTAAG': '6',
        'ATAAGCTC': '7', 'ACCGTACG': '7',
        'ATTCATGG': '8', 'TATAGTCG': '8',
        'ATCCGCGA': '9', 'TGGGCATC': '9',
        'ATCGCATA': '10', 'TACCTAGA': '10',
        'CCGTTCTA': '11', 'GCTGCATG': '11',
        'TGGCGCGC': '12', 'GTCATATG': '12',
        'TGTCTGAA': '13', 'ATATTGGC': '13',
        'CTGTCCCG': '14', 'CTAAGGGA': '14',
        'AATTTCTC': '15', 'TCGTTTCG': '15',
        'CGCGACTA': '16', 'GAATAATG': '16',
        'TGAAGCAA': '17', 'ACTGCGCA': '17',
        'TTATTCTG': '18', 'GCTTATAG': '18',
        'GTTCAACA': '19', 'ATCATGCA': '19',
        'ACGCCGGC': '20', 'ACGTTAAC': '20',
        'TTGTCTTA': '21', 'CCATCTTG': '21',
        'TACGGTTA': '22', 'CATAGCTA': '22',
        'TTGGGAGA': '23', 'GAGGTTGA': '23',
        'TGCTTGGG': '24', 'GCACTGAC': '24',
        'TAAATATC': '25', 'TTCATCGC': '25',
        'CACAATTG': '26', 'GAAATTAG': '26',
        'GTGCTAGC': '27', 'AGGATTAA': '27',
        'CGCCCGGA': '28', 'AATAGAAC': '28',
        'GCTCGCGG': '29', 'TCTTAATC': '29',
        'CTTTGGTC': '30', 'TAATACGC': '30',
        'TTCCGATC': '31', 'GTTTGTGA': '31',
        'TTCGCTAC': '32', 'CGAACGTC': '32',
        'AGCGAAAC': '33', 'GGTTCTTC': '33',
        'AAATAGCA': '34', 'GCAAATTC': '34',
        'CGTCTAGG': '35', 'GCTATGCG': '35',
        'GCCGTGTA': '36', 'CTACCCTA': '36',
        'CGCTTAAA': '37', 'GTGGGTTC': '37',
        'GACCTTTC': '38', 'GTCCGTAG': '38',
        'GGTGGAGC': '39', 'TGCGATCG': '39',
        'TACTCGAA': '40', 'TATCCGGG': '40',
        'CATTTGGA': '41', 'AGGTAATA': '41',
        'GAGCACAA': '42', 'CGTGGTTG': '42',
        'GTCGCGCG': '43', 'GACAAAGC': '43',
        'GTTACGTA': '44', 'GGGCGATG': '44',
        'CTATTTCA': '45', 'ATCTATAA': '45',
        'ACTATATA': '46', 'GCCCATGA': '46',
        'TCACTTTA': '47', 'CTGAAAGG': '47'
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
        combined_barcode = "{}{}_{}".format(*fragments)

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


def write_merged_splitseq_matrix(matrix_filename, destination):
    matrix_filename = Path(matrix_filename)
    destination = Path(destination)
    if not destination.exists():
        destination.mkdir()

    if not destination.is_dir():
        IOError("We can only write results to a directory {}".format(destination))

    source = matrix_filename.parent
    if not source.is_dir():
        IOError("Source {} needs to be a directory".format(source))

    matrix, cells = load_merged_parse_matrix(matrix_filename)

    with open(destination / "barcodes.tsv", "wt") as outstream:
        for barcode in cells:
            outstream.write(barcode)
            outstream.write(os.linesep)

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
