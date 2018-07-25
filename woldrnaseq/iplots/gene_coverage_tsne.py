#!/usr/bin/python3

import argparse

import pandas

from sklearn.manifold import TSNE

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

    plot = TSNEGeneCoverage(experiments, libraries)
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
    
class TSNEGeneCoverage(Coverage):
    """Plot showing the coverage distribution, and highlighting 
    """
    def __init__(self, experiments, libraries):
        super().__init__(experiments, libraries)
        self._tsne_cache = {}


    def get_tsne(self, experiment_name):
        if experiment_name in self._tsne_cache:
            return self._tsne_cache[experiment_name]
        else:
            library_ids = self.library_names
            subset = self.coverage[library_ids]
            tsne_subset = pandas.DataFrame(TSNE(init='pca').fit_transform(subset))
            tsne_subset.columns=['x', 'y']
            tsne_subset['library'] = subset.index
            self._tsne_cache[experiment_name] = tsne_subset
            return tsne_subset

        
    def make_plot(self, experiment_name=None):
        if experiment_name is not None:
            self.experiment_name = experiment_name
        tsne_subset = self.get_tsne(self.experiment_name)

        source = ColumnDataSource(tsne_subset)
        hover = HoverTool(tooltips=[('index', '$index'),
                                    ('(x,y)', '($x, $y)'),
                                    ('library', '@library'),
        ])

        f = figure(toolbar_location='above')
        f.add_tools(hover)
        points = f.scatter(x='x', y='y', source=source, size=10)
        return f

    def app_layout(self):
        controls = widgetbox([self.experiments_combo], width=200)
        self._layout = row(self.make_plot(self.experiment_name), controls)
        return self._layout

    def update_experiment(self, attr, old, new):
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

