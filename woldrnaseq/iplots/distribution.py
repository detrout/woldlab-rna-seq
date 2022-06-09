#!/usr/bin/env python3
from __future__ import print_function, unicode_literals, division

import argparse
import logging
import numpy

from bokeh.io import save
from bokeh.layouts import row, widgetbox
from bokeh.models import HoverTool, Select
from bokeh.plotting import figure, curdoc, ColumnDataSource
from bokeh import palettes
from bokeh import resources

from ..common import (
    add_debug_arguments,
    configure_logging,
)
from woldrnaseq.models import (
    load_experiments,
    load_library_tables,
    load_all_distribution
)

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
            print('{} was not found in {}'.format(args.use_experiment, ', '.join(list(experiments.index))))
            return None
    plot = DistributionPlot(experiments, libraries)
    return plot


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiments', action='append', default=[], help='experiments table')
    parser.add_argument('-l', '--libraries', action='append', default=[], help='library information tables')
    parser.add_argument('-n', '--use-experiment', help='plot specific experiment name')
    #parser.add_argument('-r', '--remove', nargs='*', action='append',
    #                    help='Libraries to filter out')
    add_debug_arguments(parser)
    return parser


class DistributionPlot:
    def __init__(self, experiments, libraries):
        self.experiments = experiments
        self.libraries = libraries
        self._distribution = load_all_distribution(self.libraries)

        self.experiments_combo = Select(title='Experiments',
                                        value=self.experiment_names[0],
                                        options=self.experiment_names)
        self.experiments_combo.on_change('value', self.update_experiment)
        self._layout = None

    @property
    def experiment_name(self):
        return self.experiments_combo.value

    @experiment_name.setter
    def experiment_name(self, value):
        self.experiments_combo.value = value

    @property
    def experiment_names(self):
        return sorted(self.experiments.index)

    @property
    def library_names(self):
        return self.experiments.loc[self.experiment_name].replicates

    def make_plot(self, experiment_name=None, title=None):
        if experiment_name is not None:
            self.experiment_name = experiment_name
        subset = self._distribution.loc[self.library_names]
        subset.index.name = 'library_id'
        subset = subset.reset_index()
        logger.debug("Using subset of {} records out of {}".format(
            subset.shape[0],
            self._distribution.shape[0]))

        source = ColumnDataSource(subset)
        categories = ['Exonic', 'Intronic', 'Intergenic', ]
        tooltips = [
            ('library_id', '@library_id'),
            ('Exonic', '@Exonic{%0.2f}'),
            ('Intronic', '@Intronic{%0.2f}'),
            ('Intergenic', '@Intergenic{%0.2f}'),
        ]
        hover = HoverTool(tooltips=tooltips)

        plot = figure(
            x_range=list(subset['library_id']),
            title="Distribution for {}".format(self.experiment_name)
        )
        plot.add_tools(hover)
        plot.vbar_stack(
            categories,
            x='library_id',
            width=0.5,
            color=palettes.Set1[3],
            source=source,
            legend_label=categories,
        )
        plot.legend.location = (10, plot.plot_height/2.0)
        plot.add_layout(plot.legend[0], 'right')
        plot.xgrid.grid_line_color = None
        plot.axis.minor_tick_line_color = None
        plot.xaxis.major_label_orientation = numpy.pi/4
        return plot

    def update_experiment(self, attr, old, new):
        if self._layout is not None:
            self._layout.children[1] = self.make_plot()

    def app_layout(self):
        controls = widgetbox([self.experiments_combo])
        f = self.make_plot()
        if f is not None:
            self._layout = row(controls, f)
            self._layout.sizing_mode = 'scale_both'
        return self._layout

    def static_layout(self):
        experiments_names = list(self.experiments.index)
        name = experiments_names[0]
        self._layout = self.make_plot()
        return self._layout


if __name__ == '__main__':
    plot = main()
    if plot is not None:
        curdoc().add_root(plot.static_layout())
        save(curdoc(),
             filename='distribution.html',
             resources=resources.CDN,
             title=plot.experiment_name)
elif __name__.startswith('bk_script'):
    plot = main()
    if plot is not None:
        curdoc().add_root(plot.app_layout())
