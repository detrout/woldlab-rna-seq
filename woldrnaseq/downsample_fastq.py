#!/usr/bin/python3
"""Downsample fastqs to target number
"""
from argparse import ArgumentParser
import numpy.random
import os
import sys
from subprocess import Popen, PIPE
import shutil
import requests
import gzip


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('filenames', nargs='+')
    parser.add_argument('-s', '--seed', type=int,
                        help='set seed for random number generator')
    parser.add_argument('--reads', type=int,
                        help='set total number of reads')
    parser.add_argument('-t', '--target', type=int,
                        help='target number of reads to include')
    parser.add_argument('-o', '--output')
    args = parser.parse_args(cmdline)

    if args.seed:
        numpy.random.seed(args.seed)

    if args.output is not None:
        output = GzipExt(args.output, 'wt')
    else:
        output = sys.stdout

    if args.reads is None:
        total = count_reads(args.filenames)
        print(f'Phase 1 counted {total} reads', file=sys.stderr)
    else:
        total = args.reads
        print(f'Phase 1 assumed {total} reads', file=sys.stderr)

    fraction = args.target / total
    included = 0
    seen = 0
    limit = int(total * 1.1)
    for read in read_fastqs(args.filenames):
        seen += 1
        if numpy.random.rand() <= fraction:
            included += 1
            assert len(read) == 4
            for line in read:
                output.write(line)
        if seen > limit:
            break

    fraction_included = included / total
    print(f'{included}/{total}={fraction_included:.3}',
          file=sys.stderr)

    if args.output is not None:
        output.close()


class GzipExt:
    def __init__(self, filename, mode):
        assert mode == 'wt'
        gzip_command = shutil.which('gzip')
        self.outstream = open(filename, mode)
        self.process = Popen([gzip_command], stdin=PIPE, stdout=self.outstream)

    def write(self, block):
        self.process.stdin.write(block)

    def close(self):
        self.process.stdin.flush()
        self.process.stdin.close()
        self.outstream.flush()
        self.outstream.close()


def read_fastqs(filenames):
    for filename in filenames:
        if os.path.exists(filename):
            yield from read_fastq_local(filename)
        else:
            yield from read_fastq_remote(filename)


def read_fastq_local(filename):
    gzip_command = shutil.which('gzip')

    if gzip_command is not None:
        proc = Popen([gzip_command, '-d', '-c', filename], stdout=PIPE)
        read = []
        for i, line in enumerate(proc.stdout):
            read.append(line)
            if len(read) == 4:
                yield read
                read = []


def read_fastq_remote(url):
    response = requests.get(url, stream=True)
    stream = gzip.open(response.raw)

    read = []
    for i, line in enumerate(stream):
        read.append(line)
        if len(read) == 4:
            yield read
            read = []


def count_reads(filenames):
    total = 0
    for read in read_fastqs(filenames):
        total += 1
    return total


if __name__ == '__main__':
    main()
