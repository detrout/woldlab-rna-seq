from __future__ import print_function, unicode_literals, division

import argparse
from collections import OrderedDict
import pandas
import math
import numpy
import itertools

from bokeh.io import save
from bokeh.layouts import row, column, widgetbox
from bokeh.models import (
    ColorBar,
    LinearColorMapper,
    Select
)
from bokeh.plotting import figure, curdoc
from bokeh import palettes

from woldrnaseq.models import (
    load_experiments,
    load_correlations
)


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    experiments = load_experiments(args.experiments)
    #libraries = load_library_tables(args.libraries)
    return ScoreCorrelationPlot(experiments)


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-e', '--experiments', action='append', default=[], help='experiments table')
    parser.add_argument('-n', '--use-experiment', help='plot specific experiment name')
    #parser.add_argument('-r', '--remove', nargs='*', action='append',
    #                    help='Libraries to filter out')
    return parser


class ScoreCorrelationPlot:
    def __init__(self, experiments):
        """Try to intellgently format our heatmap.
        """
        #, experiment_name, vmin=None, vmax=None, cmap="coolwarm"):
        self.experiments = experiments
        self.experiment_names = sorted(self.experiments.index)
        self._scores = {}

        self.experiments_combo = Select(
            title='Experiments',
            value=self.experiment_names[0],
            options=self.experiment_names)
        self.experiments_combo.on_change('value', self.update_experiment)

        self._layout = None

    @property
    def experiment_name(self):
        """Return name of the currently selected experiment
        """
        return self.experiments_combo.value

    @experiment_name.setter
    def experiment_name(self, value):
        self.experiments_combo.value = value

    def make_plot(self, experiment_name=None):
        if experiment_name is not None:
            self.experiment_name = experiment_name

        score = self._scores.get(
            self.experiment_name,
            load_correlations(self.experiments.loc[self.experiment_name]))
        spearman = score['rafa_spearman']
        
        factors = list(spearman.index)
        x = list(itertools.chain(*[itertools.repeat(name, len(factors)) for name in factors]))
        y = list(itertools.chain(*itertools.repeat(factors, len(factors))))
        colors = []
        #spearman_min = numpy.min(spearman.unstack())
        #spearman_max = numpy.max(spearman.unstack())
        spearman_min = 0
        spearman_max = 1
        spearman_range = spearman_max-spearman_min
        for c in spearman.unstack():
            c = (c - spearman_min)/spearman_range * 255
            if pandas.isnull(c):
                colors.append('#ffffff')
            else:
                colors.append(palettes.Plasma256[int(c)])

        plot = figure(
            title='Correlations for {}'.format(self.experiment_name),
            x_range=factors,
            y_range=factors,
            plot_width=600,
            x_axis_location='below',
            toolbar_location='above',
        )
        plot.xaxis.major_label_orientation = math.pi/2
        plot.rect(x, y, color=colors, width=1, height=1)

        mapper = LinearColorMapper(palette=palettes.Plasma256, low=spearman_min, high=spearman_max)
        cb = ColorBar(color_mapper=mapper, location=(0,0))
        plot.add_layout(cb, 'right')

        hist = self.make_score_histogram(spearman)
        layout = column(row(plot, hist))
        return layout

    def make_score_histogram(self, scores, bins=10):
        triangle = score_upper_triangular(scores)
        hist, edges = numpy.histogram(triangle, bins=bins)
        #spearman_min = numpy.min(scores.unstack())
        #spearman_max = numpy.max(scores.unstack())
        spearman_min = 0
        spearman_max = 1
        spearman_range = spearman_max-spearman_min
        vmax = max(hist)
        
        plot = figure(
            #toolbar_location=None,
            #plot_width=200,
            #plot_height=main_plot.plot_height,
            #x_range=(min(edges), max(edges)),
            x_range=(0, 1),
            #y_range=(0, vmax),
            #min_border=10,
            #y_axis_location="right"
            x_axis_label="Spearman Correlation",
            y_axis_label="Count of Scores",
        )
        plot.quad(bottom=0, left=edges[:-1], right=edges[1:], top=hist)

        return plot

    def app_layout(self):
        controls = widgetbox([self.experiments_combo], width=200)
        plot = self.make_plot(self.experiment_name)
        self._layout = row(plot, controls)
        return self._layout

    def update_experiment(self, attr, old, new):
        if self._layout is not None:
            self._layout.children[0] = self.make_plot()

    def make_plot_matplotlib(self):
        score = scores[score_name]
        #score.to_csv(experiment_name + '_' + score_name + '.csv')
        figure, ax = pyplot.subplots(1,1)
        filename = score_name + '-' + experiment_name + '.png'
        ax.set_title('Pairwise {} for {}'.format(score_name, experiment_name))
        ticks = range(len(score.columns))
        cax = ax.imshow(score, cmap='coolwarm', interpolation='none', origin='lower')
        cax.axes.set_xticks(ticks)
        cax.axes.set_yticks(ticks)
        cax.axes.set_xticklabels(score.columns, rotation=90)
        cax.axes.set_yticklabels(score.columns)

        divider = make_axes_locatable(ax)
        div_ax = divider.append_axes("right", size='5%', pad=0.05)
        pyplot.colorbar(cax, cax=div_ax)
        figure.savefig(filename)

        return filename


def make_correlation_histogram(scores, score_name, experiment_name):
    scores = score_upper_triangular(scores[score_name])
    plot = figure(title='{} scores for {}'.format(
        score_name, experiment_name))
    plot.vbar(scores, bins=min(10, len(scores)))
    return plot


def score_upper_triangular(df):
    """Return the cells from the upper triangular indicies of a data frame.
    """
    scores = []
        scores.append(df.ix[i,j])
    for i, j in zip(*numpy.triu_indices(len(df), k=1)):
    return scores


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
