#!/usr/bin/env python3
from __future__ import print_function, unicode_literals, division

import argparse
import logging
import pandas
import numpy

from bokeh.io import save
from bokeh.layouts import row, widgetbox
from bokeh.models import HoverTool, Legend, LegendItem, Select
from bokeh.plotting import figure, curdoc, ColumnDataSource
from bokeh import palettes

from woldrnaseq.common import (
    add_debug_arguments,
    configure_logging,
)
from woldrnaseq.models import (
    load_experiments,
    load_library_tables,
    load_quantifications,
)
from woldrnaseq.gtfcache import GTFCache, protein_coding_gene_ids

logger = logging.getLogger(__name__)


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    experiments = load_experiments(args.experiments)
    libraries = load_library_tables(args.libraries)
    if args.use_experiment:
        try:
            experiments = experiments.loc[[args.use_experiment]]
        except KeyError:
            logger.error('{} was not found in {}'.format(args.use_experiment, ', '.join(list(experiments.index))))
            return None
    plot = GenesDetectedPlot(experiments, libraries, args.genome_dir, args.quantification)

    if __name__ == '__main__':
        curdoc().add_root(plot.static_layout())
        save(curdoc(), args.output, title=plot.title)

    return plot


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genome-dir', required=True,
                        help='root for genome indexes')
    parser.add_argument('-e', '--experiments', action='append', default=[], help='experiments table')
    parser.add_argument('-l', '--libraries', action='append', default=[], help='library information tables')
    parser.add_argument('-n', '--use-experiment', help='plot specific experiment name')
    parser.add_argument('-q', '--quantification', default='TPM', help='Specify quantification type to use')
    parser.add_argument('-o', '--output', default='genesdetected.html', help='output filename')
    parser.add_argument('filenames', nargs='*',
                        help='Combined quantification file: libraries by genes')
    add_debug_arguments(parser)
    return parser


def load_csv_quantification_file(filename):
    all_quantifications = pandas.read_csv(filename, header=0, index_col=0)
    if 'gene_name' in all_quantifications.columns:
        columns = [c for c in all_quantifications.columns if c != 'gene_name']
        all_quantifications = all_quantifications[columns]
    return all_quantifications


DEFAULT_BINS = [0.1, 1, 2, 5, 10, 50, 500, 5000, 1e9]


class Quantifications:
    DEFAULT_BINS = {
        'FPKM': [0.1, 1, 2, 5, 10, 50, 500, 5000, 1e9],
        'TPM': [0.1, 1, 2, 5, 10, 50, 500, 5000, 1e9],
    }

    def __init__(self, experiments, sep='\t', analysis_root=None):
        self.name = None
        self.experiments = load_experiments(experiments, sep)
        self.quantification_name = None
        self.quantification = None


def bin_library_quantification(quantification, quantification_name, bins=None):
    """Bin a library quantification file

    default bins are [0.1, 1, 2, 5, 10, 50, 500, 5000, 1e9]
    """
    if bins is None:
        bins = DEFAULT_BINS

    histogram = {}
    for col in quantification:
        histogram[col], _ = numpy.histogram(quantification[col], bins)

    histogram = pandas.DataFrame(
        histogram,
        columns=quantification.columns,
        index=['{}_{}'.format(str(x).replace('.', '_'), quantification_name) for x in bins[:-1]])

    libs_vs_threshold = histogram.reindex(histogram.index[::-1]).T
    libs_vs_threshold.index.name = 'library_id'
    return libs_vs_threshold


class GenesDetectedPlot:
    def __init__(self, experiments, libraries, genome_dir, quantification_name='FPKM'):
        self.experiments = experiments
        self.experiment_names = sorted(self.experiments.index)
        self.libraries = libraries
        self.genome_dir = genome_dir
        self._gtf_cache = GTFCache(self.libraries, self.genome_dir)
        self.quantification_name = quantification_name
        self.binned_quantifications = {}
        self.bin_names = {}
        self.experiments_combo = Select(
            title='Experiments',
            value=self.experiment_names[0],
            options=self.experiment_names)
        self.experiments_combo.on_change('value', self.update_plot)
        self._layout = None

        self.load_all_quantifications(self.experiments)

    @property
    def experiment_name(self):
        """Return name of the currently selected experiment
        """
        return self.experiments_combo.value

    @experiment_name.setter
    def experiment_name(self, value):
        self.experiments_combo.value = value

    @property
    def title(self):
        return self.experiment_name + ' genes detected'

    def load_all_quantifications(self, experiments):
        for experiment_name, experiment_row in experiments.iterrows():
            all_quant = load_quantifications(experiment_row, self.quantification_name)
            if all_quant is None:
                continue

            annotation = self._gtf_cache[all_quant.columns[-1]]
            protein_genes = protein_coding_gene_ids(annotation)

            protein_quant = all_quant.loc[protein_genes]

            binned = bin_library_quantification(protein_quant, self.quantification_name)
            self.bin_names[experiment_name] = binned.columns
            binned['total'] = binned.sum(axis=1)
            self.binned_quantifications[experiment_name] = binned

    def make_plot(self, experiment_name=None, title=None):
        self.experiment_name = experiment_name if experiment_name is not None else self.experiment_names[0]

        binned = self.binned_quantifications[self.experiment_name]
        bin_names = self.bin_names[self.experiment_name]
        friendly_names = ['{} {}'.format(k, self.quantification_name) for k in reversed(DEFAULT_BINS[:-1])]
        tooltips = [('library_id', '@library_id'),
                    ('Total', '@total')]
        for name, column_name in zip(friendly_names, bin_names):
            tooltips.append((name, '@' + column_name))
        hover = HoverTool(tooltips=tooltips)

        title = self.title if title is None else title
        source = ColumnDataSource(binned)
        f = figure(title=title,
                   x_range=list(binned.index),
                   plot_width=1200,
        )
        f.add_tools(hover)

        f.vbar_stack(
            bin_names,
            x='library_id',
            width=0.5,
            source=source,
            color=palettes.Oranges8,
        )

        legend_items = []
        for label, color in zip(friendly_names, palettes.Oranges8):
            legend_items.append(
                LegendItem(
                    label=label,
                    renderers=[f.square([1], [1], color=color, line_color='black')]))
        legend = Legend(items=legend_items, location=(30, 0))
        f.add_layout(legend, 'right')
        f.y_range.start = 0
        f.x_range.range_padding = 0.1
        f.xgrid.grid_line_color = None
        f.axis.minor_tick_line_color = None
        f.outline_line_color = None
        #f.legend.location = "top_left"
        #f.legend.orientation = "vertical"
        f.xaxis.major_label_orientation = numpy.pi/4

        return f

    def update_plot(self, attr, old, new):
        if self._layout is not None:
            experiment_name = self.experiments_combo.value
            self._layout.children[1] = self.make_plot(experiment_name)

    def app_layout(self):
        controls = widgetbox([self.experiments_combo])
        f = self.make_plot(self.experiments_combo.value)
        if f is not None:
            self._layout = row(controls, f)
            self._layout.sizing_mode = 'scale_both'
        return self._layout

    def static_layout(self):
        name = self.experiment_name
        f = self.make_plot(name)

        return f


if __name__ == '__main__':
    plot = main()
elif __name__.startswith('bk_script'):
    plot = main()
    if plot is not None:
        curdoc().add_root(plot.app_layout())
