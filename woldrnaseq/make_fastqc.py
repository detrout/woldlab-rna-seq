#!/usr/bin/python3
from __future__ import absolute_import, print_function

from argparse import ArgumentParser
import itertools
import logging
import os
import re

from . import models
from .common import (
    add_debug_arguments,
    configure_logging,
    find_fastqs,
    get_seperator,
)

logger = logging.getLogger('make_fastqc')


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-l', '--library', action='append')
    parser.add_argument('-c', '--contamination',
                        help='contamination filename')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    add_debug_arguments(parser)
    args = parser.parse_args(cmdline)

    configure_logging(args)

    if args.library is None or len(args.library) == 0:
        parser.error('Library metadata file required')

    sep = get_seperator(args.sep)
    libraries = models.load_library_tables(args.library, sep)

    read1 = dict(find_fastqs(libraries, 'read_1'))
    if 'read_2' in libraries.columns:
        read2 = dict(find_fastqs(libraries, 'read_2'))
    else:
        read2 = {}

    fastqs = set()
    for lib_id, lib_fastqs in itertools.chain(read1.items(), read2.items()):
        fastqs.update(lib_fastqs)

    generate_fastqc_condor(fastqs, args.contamination)


def generate_fastqc_condor(fastqs, contamination):
    header = """universe=vanilla

output=fastqc.$(process).out
error=fastqc.$(process).out
log=fastqc.log

executable=/usr/bin/fastqc
transfer_executable=false

request_memory=4G

should_transfer_files=IF_NEEDED
"""
    args = [header]
    for fastq in sorted(fastqs):
        arg = []
        qc_report = make_qc_filename(fastq)
        if not os.path.exists(qc_report):
            if contamination is not None:
                arg.extend(['-c', contamination])
            arg.append(fastq)
            args.append('arguments="' + ' '.join(arg) + '"')
            args.append('queue')
        else:
            logger.info("Skipping %s, fastqc report already exists", fastq)
        
    print(os.linesep.join(args))

def make_qc_filename(fastq_filename):
    path = re.sub('\.fastq(\.gz|\.bz2|\.xz)?$', '', fastq_filename)
    return path + '_fastqc.html'

if __name__ == '__main__':
    main()
