#!/usr/bin/env python3
"""Convert a gff/gtf file to an indexed HDF5 file.

It takes about 11 minutes to run on the ENCODE GENCODE v19 file,
but afterwords I can query for attributes by gene name in 10-15 ms.
"""
import argparse
import os
import collections
import logging

import pandas
import numpy
import time
from xopen import xopen

logger = logging.getLogger('gff2table')


def tokenize(cell, sep):
    """Generator returning tokens for gff/gtf attributes
    """
    SEP = 0
    STRING = 1
    QUOTED_STRING = 2

    state = SEP
    current = []
    for character in cell:
        if state == SEP:
            if character.isspace():
                continue
            elif character == sep:
                state = SEP
                yield character
            elif character == '"':
                current.append(character)
                state = QUOTED_STRING
            elif character == ';':
                state = SEP
                yield character
            else:
                current.append(character)
                state = STRING
        elif state == STRING:
            if character in (sep, ';'):
                state = SEP
                yield ''.join(current)
                current = []
                yield character
            else:
                current.append(character)
        elif state == QUOTED_STRING:
            if character == '"':
                state = SEP
                current.append('"')
                yield ''.join(current)
                current = []
            else:
                current.append(character)

    if len(current) > 0:
        if state == QUOTED_STRING:
            raise ValueError('Unclosed quote')
        else:
            yield ''.join(current)


def parse_attributes(cell, sep, ignore={}, reserved={}):
    tokens = tokenize(cell, sep)
    for term in tokens:
        name = term

        sep = next(tokens)

        value = next(tokens)

        if name not in ignore:
            if name in reserved:
                suffix = reserved[name] + 1
                reserved[name] = suffix
                name = name + str(suffix)

            if value in ('"NULL"', 'NULL', 'nan'):
                value = None
            elif value[0] == '"':
                value = value[1:-1]
            elif value[0].isdigit():
                value = int(value)

            yield name, value

        try:
            next(tokens)
        except StopIteration:
            break


class AttributesParser:
    def __init__(self, sep=' ', ignore=None, max_terms=50):
        self.index = 0
        self.sep = sep
        self.terms = collections.OrderedDict()
        self.max_string = collections.OrderedDict()
        self.max_terms = max_terms
        self.ignore = ignore if ignore is not None else []
        self.reserved = {
            'chromosome': 0,
            'source': 0,
            'type': 0,
            'start': 0,
            'stop': 0,
            'score': 0,
            'strand': 0,
            'frame': 0,
        }

    def tokenize(self, cell):
        yield from tokenize(cell, self.sep)

    def __call__(self, cell):
        attributes_count = 0
        for name, value in parse_attributes(cell, self.sep, self.ignore, self.reserved):
            if isinstance(value, str):
                prev = self.max_string.get(name, 0)
                self.max_string[name] = max(prev, len(value))

            column = self.terms.setdefault(name, {})
            column[self.index] = value
            attributes_count += 1

        assert len(self.terms) < self.max_terms
        self.index += 1
        return attributes_count


def parse_score(x):
    if x == '.':
        return numpy.nan
    else:
        if '.' in x:
            return float(x)
        else:
            return int(x)


def parse_strand(x):
    if x == '+':
        return 1
    elif x == '-':
        return -1
    else:
        return numpy.nan


def parse_phase(x):
    if x == '.':
        return numpy.nan
    else:
        phase = int(x)
        if phase < 0 or phase > 2:
            raise ValueError('Invalid frame: x')
        return phase


class GFFParser:
    def __init__(self, sep=' ', ignore=None):
        self.attribute_parser = AttributesParser(sep=sep)
        self.gtf = None
        self.ignore = ignore

    def read_gff(self, inname):
        tzero = time.monotonic()
        required_gtf_names = [
            'chromosome', 'source', 'type', 'start', 'stop',
            'score', 'strand', 'frame']
        column_names = required_gtf_names + ['attributes']

        with xopen(inname, "rt") as instream:
            gtf = pandas.read_csv(
                instream,
                sep='\t',
                header=None,
                names=column_names,
                index_col=False,
                comment="#",
                na_values='.',
                converters={
                    'chromosome': str,
                    'score': parse_score,
                    'strand': parse_strand,
                    'frame': parse_phase,
                    'attributes': self.attribute_parser,
                },
            )
        tnow = time.monotonic()
        tprev = tnow
        logger.info("Parsed in {:.3} seconds".format(tnow-tzero))
        print('gtf.shape', gtf.shape)
        for name in ['chromosome', 'type']:
            self.attribute_parser.max_string[name] = max(gtf[name].map(len))

        # drop my synthetic column counting how many records are in the
        # variant column
        gtf.drop('attributes', axis=1, inplace=True)
        # create new table with metadata attributes
        attributes = pandas.DataFrame(self.attribute_parser.terms)
        gtf_with_metadata = gtf.merge(attributes, left_index=True, right_index=True)
        tnow = time.monotonic()
        logger.info("Merged table in {:.3} seconds".format(tnow-tprev))
        tprev = tnow
        self.gtf = gtf_with_metadata

    def write_hdf5(self, outname, table_name):
        if self.gtf is None:
            raise RuntimeError('No gtf table to write. Did you call read_gff?')
        # columns to index
        data_columns = [
            "chromosome",
            "source",
            "type",
            "start",
            "stop",
            "gene_id",
            "transcript_id",
            "gene_name",
            "gene_type",
            "transcript_type"]
        tprev = time.monotonic()
        store = pandas.HDFStore(outname, mode='w', complevel=9, complib='blosc:zstd')
        store.append(
            table_name,
            self.gtf,
            data_columns=data_columns,
            min_itemsize=self.attribute_parser.max_string)
        tnow = time.monotonic()
        logger.info("Wrote table in {:.3} seconds".format(tnow-tprev))
        tprev = tnow
        store.create_table_index(table_name, optlevel=9, kind='full')
        tnow = time.monotonic()
        logger.info("Wrote index in {:.3} seconds".format(tnow-tprev))
        store.close()

    def write_table(self, outname, sep='\t', columns=None):
        if self.gtf is None:
            raise RuntimeError('No gtf table to write. Did you call read_gff?')
        gtf = self.gtf[columns] if columns is not None else self.gtf
        gtf.to_csv(outname, sep=sep)


def format_gtf_record(row, value_sep=' ', field_sep='; '):
    strand = {
        1: '+',
        -1: '-',
    }
    record = []
    record.append(str(row.chromosome))  # 0
    record.append(str(row.source))      # 1
    record.append(str(row.type))        # 2
    record.append(str(int(row.start)))  # 3
    record.append(str(int(row.stop)))   # 4
    record.append(str(int(row.score)) if not pandas.isnull(row.score) else '.')  # 5
    record.append(strand[row.strand])   # 6
    record.append(str(int(row.frame)) if not pandas.isnull(row.frame) else '.')  # 7
    attributes = []
    for name in row.index[8:]:
        if isinstance(row[name], str) or not pandas.isnull(row[name]):
            value = row[name]
            attributes.append(f'{name}{value_sep}"{value}"')
    row = field_sep.join(attributes) + field_sep.rstrip()
    record.append(row)
    return '\t'.join(record)


def make_parser():
    parser = argparse.ArgumentParser(
        description="Convert GTF/GFF file to a binary cache file"
    )
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('-d', '--debug', action='store_true', default=False)
    parser.add_argument('-n', '--name', default='gtf', help='table name in hdf5 file')
    parser.add_argument('--ignore', default=[], action='store_true',
                        help='do not add listed attributes to output h5 file')
    parser.add_argument('-o', '--output',
                        help='specify output name, defaults to name.h5')
    parser.add_argument('--output-type', choices=['hdf5', 'csv'], default='hdf5',
                        help='specify output file format')
    parser.add_argument('filename', nargs='+',
                        help='specify input GFF/GTF filename')
    parser.add_argument('--name-value-sep', default=" ",
                        help='Specify name value separator in the attribute column. '
                             '(GTF and GFF use different delimiters.)')
    return parser


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARN)

    if args.output is None:
        name, ext = os.path.splitext(args.filename[0])
        args.output = name + '.h5'

    for filename in args.filename:
        logger.info("Converting %s to %s", filename, args.output)

        gtf = GFFParser(args.name_value_sep, ignore=args.ignore)
        gtf.read_gff(filename)

    if args.output_type == 'hdf5':
        gtf.write_hdf5(args.output, args.name)
    elif args.output_type == 'csv':
        gtf.write_table(args.output, args.name)
    else:
        raise RuntimeError('Unhandled output-type {}'.format(args.output_type))


if __name__ == '__main__':
    main()
