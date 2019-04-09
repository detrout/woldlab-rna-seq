#!/usr/bin/python3

import argparse
from itertools import cycle
import logging
import os

import pandas

from bokeh.io import save
from bokeh.layouts import row, widgetbox
from bokeh.models import Legend
from bokeh.plotting import figure, curdoc
from bokeh import resources, palettes

from woldrnaseq.models import (
    load_library_tables,
    load_all_gene_coverage,
)


LOGGER = logging.getLogger("gene_coverage_detail")


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    libraries = []
    if args.libraries:
        libraries = load_library_tables(args.libraries, analysis_root=args.root)
        LOGGER.info("loaded %d libraries", len(libraries))

    if len(libraries) == 0 and len(args.gene_list) == 0:
        parser.error('Please specify a libraries to process')

    with open(args.gtf, 'rt') as stream:
        gene_types = readGeneTypes(stream)
        LOGGER.info("Loaded %s gene types", len(gene_types))

    coverage_by_type = {}
    counts_by_type = {}
    for gene_coverage_table in load_all_gene_coverage(libraries, args.gene_list, args.gene_normalization):
        coverage, counts = sum_gene_coverage_by_type(gene_types, gene_coverage_table)
        coverage_by_type[coverage.name] = coverage
        counts_by_type[coverage.name] = counts

    LOGGER.info('Preparing plot class')
    plot = GeneCoverageDetail(coverage_by_type, counts_by_type, args.gene_normalization)

    if args.save:
        for library_id in plot:
            # avoid names that cause problems for files systems
            assert not library_id.startswith('..')
            assert '/' not in library_id
            assert '\\' not in library_id
            filename = '{}_gene_coverage_detail.html'.format(library_id)
            pathname = os.path.join(args.output_dir, filename)
            LOGGER.info("Saving plot for %s to %s", library_id, pathname)
            save(plot.make_plot(library_id),
                 pathname,
                 resources=resources.CDN,
                 title=library_id,
            )
    return plot


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf', required=True, help='Specify GTF annotation file')
    parser.add_argument('gene_list', nargs='*', help='specify what gene list file to load')
    parser.add_argument('-e', '--experiments', action='append', default=[], help='experiments table')
    parser.add_argument('-l', '--libraries', action='append', default=[], help='library information tables')
    #parser.add_argument('-n', '--use-experiment', help='plot specific experiment name')
    #parser.add_argument('-r', '--remove', nargs='*', action='append',
    #                    help='Libraries to filter out')
    parser.add_argument('--gene-normalization', choices=['None', 'max'], default=None)
    parser.add_argument('--root', default=None,
                        help='analysis_dir will be relative to this path '
                        'instead of library.txt file')
    parser.add_argument('--save', default=False, action='store_true',
                        help='save coverage plots')
    parser.add_argument('-o', '--output-dir', default='.',
                        help='directory to write bokeh plots to <library_id>_gene_coverage_detail.html')
    parser.add_argument('-d', '--debug', action='store_true', default=False,
                        help='Use debug logging level')

    return parser


def readGeneTypes(stream):
    gene_types = {}
    for i, line in enumerate(stream):
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        if fields[2] != 'exon':
            continue
        geneID = getGFFAttributeValueByKey(fields[8], 'gene_id')
        source = fields[1]
        if source == 'spikein':
            gene_type = 'spikein'
        else:
            gene_type = getGFFAttributeValueByKey(fields[8], 'gene_type')
            if gene_type is None:
                gene_type = "n/a"
        gene_types.setdefault(gene_type, []).append(geneID)
    LOGGER.debug('Read %d lines from gtf stream', i)
    return gene_types


def getGFFAttributeValueByKey(field, name):
    start = field.find(name)
    if start == -1:
        return None
    start += len(name)
    start = field.index('"', start)
    start += 1
    end = field.index('"', start)
    return field[start:end]


def sum_gene_coverage_by_type(gene_types, gene_coverage):
    """Compute the coverage by gene type.

    :Parameters:
    - gene_types: (dict) gene_id to gene_type
    - gene_coverage: (DataFrame) gene_id by gene coverage percentile.
    :Returns:
    - sums of gene type bins, number of genes per bin.
    """
    LOGGER.info("Computing coverage for %s", gene_coverage.name)
    sums = {}
    sizes = {}
    for gene_type in gene_types:
        genes = set(gene_types[gene_type])
        current_genes = gene_coverage.reindex(genes).dropna()
        LOGGER.debug('len(%s)=%s shape=%s', gene_type, len(genes), current_genes.shape)
        if len(current_genes) > 0:
            sizes[gene_type] = len(current_genes)
            sums[gene_type] = current_genes.sum()
    # The leading space is so when sort the names Total is always first
    sizes[' Total'] = sum(sizes.values())
    sums[' Total'] = sum(sums.values())
    counts_by_type = pandas.Series(sizes)
    counts_by_type.name = gene_coverage.name
    coverage_by_type = pandas.DataFrame(sums, index=range(100))
    coverage_by_type.name = gene_coverage.name
    LOGGER.debug('coverage for %s done %s', coverage_by_type.name, coverage_by_type.shape)
    return coverage_by_type, counts_by_type


class GeneCoverageDetail:
    """Plot showing the coverage distribution for gene type classes

    :Parameters:
    - gene_types
    """
    def __init__(self, coverage_by_type, counts_by_type, normalization):
        self._coverage_by_type = coverage_by_type
        self._counts_by_type = counts_by_type
        self._library_ids = sorted(self._coverage_by_type.keys())
        self._library_id = None
        self._layout = None
        self._normalization = normalization

    @property
    def library_id(self):
        if self._library_id is None:
            self._library_id = self._library_ids[0]
        return self._library_id

    @library_id.setter
    def library_id(self, value):
        if value in self._library_ids:
            self._library_id = value
        else:
            raise KeyError("Invalid library id")

    def __iter__(self):
        for library_id in self._library_ids:
            yield library_id

    def __getitem__(self, key):
        return self._coverage_by_type[key]

    def __len__(self):
        return len(self._coverage_by_type)

    def make_plot(self, library_id=None):
        """Show read depth coverage over normalized gene regions.
        """
        if library_id is None:
            library_id = self.library_id


        coverage = self._coverage_by_type[library_id]
        gene_counts = self._counts_by_type[library_id]

        if self._normalization == 'max':
            title = "{} normalized to maximum percentile".format(library_id)
        elif self._normalization in ('None', None):
            title = "{} raw coverage".format(library_id)
        else:
            raise RuntimeError("Unrecognized normalization %s", self._normalization)
        f = figure(x_axis_label="position quantile (5' to 3')",
                   y_axis_label="Read depth",
                   toolbar_location='right',
                   title=title,
                   plot_width=800,
                   plot_height=600,
        )
        colors = palettes.Category20[20]
        colorcycler = cycle(colors)
        legend_items = []
        for category in sorted(coverage):
            next_color = next(colorcycler)
            line = f.line(x=coverage.index, y=coverage[category],
                          line_color=next_color,
                          line_width=2,
            )
            legend_items.append(("{} ({})".format(category.strip(), gene_counts[category]), [line]))
        legend_y = int(20 - len(category) * 2)
        legend = Legend(items=legend_items, location=(0, legend_y))
        legend.click_policy = 'hide'
        f.add_layout(legend, 'right')
        return f

    def app_layout(self):
        controls = widgetbox([self.experiments_combo], width=200)
        self._layout = row(self.make_plot(self.library_id), controls)
        self._layout = self.make_plot()
        return self._layout

    def update_library_id(self, attr, old, new):
        if self._layout is not None:
            self._layout.children[0] = self.make_plot()


if __name__ == '__main__':
    plot = main()
elif __name__.startswith('bk_script'):
    plot = main()
    if plot is not None:
        curdoc().add_root(plot.app_layout())
