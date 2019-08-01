#!/usr/bin/python3

import argparse
from itertools import cycle

from bokeh.io import export_png, save
from bokeh.layouts import row, column, widgetbox
from bokeh.models import HoverTool, Legend, LegendItem, Select
from bokeh.plotting import figure, curdoc, ColumnDataSource
from bokeh import resources, palettes

from woldrnaseq.models import (
    load_experiments,
    load_library_tables,
)

from woldrnaseq.iplots.coverage import Coverage


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    experiments = load_experiments(args.experiments)
    libraries = load_library_tables(args.libraries)

    plot = GeneCoverage(experiments, libraries)
    plot.use_experiment(args.use_experiment)
    return plot

def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiments', action='append', default=[], help='experiments table')
    parser.add_argument('-l', '--libraries', action='append', default=[], help='library information tables')
    parser.add_argument('-n', '--use-experiment', help='plot specific experiment name')
    #parser.add_argument('-r', '--remove', nargs='*', action='append',
    #                    help='Libraries to filter out')
    return parser
    
class GeneCoverage(Coverage):
    """Plot showing the coverage distribution for every library
    """
    def __init__(self, experiments, libraries):
        super().__init__(experiments, libraries)

    def make_plot(self, experiment_name=None):
        """Show read depth coverage over normalized gene regions.
        """
        if experiment_name is not None:
            self.experiment_name = experiment_name

        library_ids = self.library_names
        subset = self.coverage[library_ids]
        source = ColumnDataSource(subset)
        f = figure(x_axis_label="position quantile (5' to 3')",
                   y_axis_label="Read depth",
                   toolbar_location='above',
                   width=900,
                   sizing_mode='scale_width',
        )
        colors = palettes.Category10[10]
        colorcycler = cycle(colors)
        legend_items = []
        for lib in subset:
            next_color = next(colorcycler)
            line = f.line(x=subset[lib].index, y=subset[lib], line_color=next_color)
            legend_items.append((lib, [line]))
        legend = Legend(items=legend_items, location=(0, -60))
        f.add_layout(legend, 'right')
        return f

    def app_layout(self):
        controls = widgetbox([self.experiments_combo], width=200)
        self._layout = row(self.make_plot(self.experiment_name), controls)
        return self._layout

    def update_experiment(self, attr, old, new):
        if self._layout is not None:
            self._layout.children[0] = self.make_plot()
        
if __name__ == '__main__':
    plot = main()
    if plot is not None:
        curdoc().add_root(plot.static_layout())
        #export_png(curdoc(), 'genesdetected.png')
        save(curdoc(), 'coverage_mean.html')
elif __name__.startswith('bk_script'):
    plot = main()
    if plot is not None:
        curdoc().add_root(plot.app_layout())

