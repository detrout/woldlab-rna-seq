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

logger = logging.getLogger('gff2table')


class AttributesParser:
    def __init__(self, sep=' ', ignore=None):
        self.index = 0
        self.sep = sep
        self.terms = collections.OrderedDict()
        self.max_string = collections.OrderedDict()
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
                elif character == self.sep:
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
                if character in (self.sep, ';'):
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

    def __call__(self, cell):
        tokens = self.tokenize(cell)
        attributes_count = 0
        for term in tokens:
            name = term

            sep = next(tokens)

            value = next(tokens)

            if name not in self.ignore:
                if name in self.reserved:
                    suffix = self.reserved[name] + 1
                    self.reserved[name] = suffix
                    name = name + str(suffix)

                if value in ('"NULL"', 'NULL', 'nan'):
                    value = None
                elif value[0] == '"':
                    value = value[1:-1]
                    prev = self.max_string.get(name, 0)
                    self.max_string[name] = max(prev, len(value))
                else:
                    try:
                        value = int(value)
                    except ValueError:
                        pass

                column = self.terms.setdefault(name, {})
                column[self.index] = value
                attributes_count += 1

            try:
                end_field = next(tokens)
            except StopIteration:
                break

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
            'score', 'strand', 'frame',]
        column_names = required_gtf_names + ['attributes']

        gtf = pandas.read_csv(
            inname,
            sep='\t',
            header=None,
            names=column_names,
            index_col=False,
            comment="#",
            na_values='.',
            converters={
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
        tprev = time.monotonic()
        store = pandas.HDFStore(outname, mode='w', complevel=9, complib='bzip2')
        store.append(table_name, self.gtf,
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
        if columns is not None:
            gtf = gtf[[columns]]
        gtf.to_csv(outname, sep=sep)


def main(cmdline=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('-d', '--debug', action='store_true', default=False)
    parser.add_argument('-n', '--name', default='gtf', help='table name in hdf5 file')
    parser.add_argument('--ignore', default=[], action='store_true',
                        help='do not add listed attributes to output h5 file')
    parser.add_argument('-o', '--output',
                        help='specify output name, defaults to name.h5')
    parser.add_argument('filename', nargs='+',
                        help='specify input GFF/GTF filename')
    parser.add_argument('--name-value-sep', default=" ",
                        help='name value seperator in attribute column')

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
        gtf.write_hdf5(args.output, args.name)


if __name__ == '__main__':
    main()
