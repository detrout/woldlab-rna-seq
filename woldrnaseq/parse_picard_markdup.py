#!/usr/bin/env python3
"""Parse metrics file produced by picard MarkDuplicates command.

e.g. java -Xmx8G -jar $(PICARD) MarkDuplicates \
        I=sample_genome.bam \
        O=sample_picard_markdup.bam \
        M=sample_picard_markdup.metrics
"""
from argparse import ArgumentParser
import pandas
from pathlib import Path

from woldrnaseq.models import load_library_tables, genome_name_from_library


def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-l', '--library', required=True, action='append',
                        help="library table to load")
    parser.add_argument('-o', '--output', help='filename to write report to')
    args = parser.parse_args(cmdline)

    libraries = load_library_tables(args.library)

    metrics = []
    for library_id, library in libraries.iterrows():
        genome_triple = genome_name_from_library(library)
        filename = library.analysis_name + '-' + genome_triple + '_picard_markdup.metrics'
        pathname = Path(library.analysis_dir) / filename
        if pathname.exists():
            picard_metric = parse_picard_metric(pathname, library_id=library_id)
            metrics.append(picard_metric)
        else:
            print('{} is missing. Skipping'.format(pathname))

    metrics = pandas.DataFrame(metrics)
    metrics.set_index('LIBRARY', inplace=True)

    if args.output:
        metrics.to_csv(args.output, sep='\t')
    else:
        print(metrics)


def parse_picard_metric(filename, library_id=None):
    with open(filename) as instream:
        return parse_picard_metric_stream(instream, library_id)


def parse_picard_metric_stream(stream, library_id=None):
    IGNORE = 0
    METRIC_HEADER = 1
    METRIC_BODY = 2

    state = IGNORE
    for line in stream:
        if line.startswith('## METRICS CLASS'):
            state = METRIC_HEADER
        elif state == METRIC_HEADER:
            state = METRIC_BODY
            header = line.rstrip().split('\t')
        elif state == METRIC_BODY:
            body = line.rstrip().split('\t')
            break

    series = pandas.Series(body, index=header[:len(body)])

    if library_id is not None:
        series.LIBRARY = library_id

    return series


if __name__ == "__main__":
    main()
