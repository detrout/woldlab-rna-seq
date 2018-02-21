#!/usr/bin/python3

import argparse

from bokeh.io import export_png, save
from bokeh.layouts import row, column, widgetbox
from bokeh.models import HoverTool, Legend, LegendItem, Select
from bokeh.plotting import figure, curdoc, ColumnDataSource
from bokeh import resources, palettes

from woldrnaseq.models import (
    load_experiments,
    load_library_tables,
)

# This can't be a relative import because if bokeh serve loads it, its
# not in a proper module
from woldrnaseq.bokeh.coverage import Coverage

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    experiments = load_experiments(args.experiments)
    libraries = load_library_tables(args.libraries)

    plot = MeanGeneCoverage(experiments, libraries)
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
    
class MeanGeneCoverage(Coverage):
    """Plot showing the coverage distribution, and highlighting 
    """
    def __init__(self, experiments, libraries):
        super().__init__(experiments, libraries)

        self.library_combo = Select(
            title='Libraries',
            value='None',
            options=['None'] + self.library_names)
        self.library_combo.on_change('value', self.update_library)

    @property
    def library_id(self):
        """return currently selected library id"""
        return self.library_combo.value
        
    def make_plot(self, experiment_name, library_id=None):
        """Make a coverage summary plot

        :Parameters:
           experiment_name - name of experiment to plot
           library_id - optional library id to highlight.
        """
        library_ids = self.experiments.loc[experiment_name].replicates
        subset = self.coverage[library_ids]

        mean = subset.mean(axis=1)
        dev = subset.std(axis=1)
        ymax = mean+dev
        ymin = mean-dev
        f = figure(toolbar_location='above')
        f.line(x=subset.index, y=mean, legend='mean', line_color=palettes.Blues4[0])
        f.quad(left=subset.index,
               right=subset.index+1,
               top=ymax,
               bottom=ymin,
               alpha=0.5,
               color=palettes.Blues4[2],
               legend='+/- stdev')
        if library_id is not None and library_id != 'None' and library_id in self.coverage.columns:
            f.line(x=subset.index,
                   y=subset[library_id],
                   legend=library_id,
                   line_color=palettes.Reds4[0])

        f.legend.location='bottom_center'
        return f

    def app_layout(self):
        controls = widgetbox([self.experiments_combo, self.library_combo], width=200)
        self._layout = row(self.make_plot(self.experiment_name), controls)
        return self._layout

    def update_experiment(self, attr, old, new):
        old_value = self.library_combo.value
        self.library_combo.value = 'None'

        self.library_combo.options = ['None'] + self.library_names
        self.update_library('value', old_value, self.library_combo.value)

    def update_library(self, attr, old, new):
        self._layout.children[0] = self.make_plot(
            self.experiment_name,
            self.library_id)
        
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

