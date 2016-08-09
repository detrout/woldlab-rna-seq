#!/usr/bin/python3
import argparse
import itertools
import logging

from jinja2 import Environment, PackageLoader
import matplotlib
matplotlib.use('Agg')
import numpy
from matplotlib import pyplot
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas
from bokeh import mpl
from bokeh.models import HoverTool
from bokeh.plotting import figure, ColumnDataSource
from bokeh.charts import Line, Bar, Dot, HeatMap, BoxPlot
from bokeh.charts.attributes import CatAttr
from bokeh.embed import components

import models
import madqc

logger = logging.getLogger('QC Report')

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)
    env = Environment(loader=PackageLoader('woldrnaseq', 'templates'))
    sep = models.get_seperator(args.sep)

    experiments = models.load_experiments([args.experiments], sep)
    libraries = models.load_library_tables([args.libraries], sep)

    samstats = models.load_all_samstats(libraries)
    distribution = models.load_all_distribution(libraries)
    coverage = models.load_all_coverage(libraries)


    seen_libraries = set()
    experiment_report = {}
    transcript_library_plots = []
    plots = {}
    plot_handle = itertools.count()
    for experiment in experiments:
        scores = models.load_correlations(experiment)
        quantifications = models.load_quantifications(experiment, args.quantification)

        coverage_handle = str(next(plot_handle))
        plots[coverage_handle] = make_coverage_plot(coverage, experiments, experiment)
        distribution_handle = str(next(plot_handle))
        plots[distribution_handle] = make_distribution_plot(distribution, experiments, experiment)
        spearman_filename = make_correlation_heatmap(scores, 'rafa_spearman', experiment)
        spike_variance_handle = str(next(plot_handle))
        plots[spike_variance_handle] = make_spikein_variance_plot(quantifications, experiments, experiment)

        library_ids = experiments[experiment]
        seen_libraries.update(set(library_ids))
        for library_id in library_ids:
            transcript_handle = str(next(plot_handle))
            plots[transcript_handle] = make_spikein_per_transcript_plot(
                quantifications, library_id, args.quantification)
            transcript_library_plots.append(transcript_handle)

        experiment_report[experiment] = {
            'samstats': samstats.select(lambda x: x in library_ids).to_html(),
            'spearman': scores.rafa_spearman.to_html(),
            'coverage': coverage_handle,
            'distribution': distribution_handle,
            'spearman_plot': spearman_filename,
            'spike_variance': spike_variance_handle,
        }

    spare_libraries = set(libraries.index).difference(seen_libraries)
    print(libraries.index)
    print(seen_libraries)
    print(spare_libraries)

    script, plot_divs = components(plots)

    template = env.get_template('rnaseq.html')
    page = template.render(
        experiments=experiments,
        experiment_report=experiment_report,
        transcript_library_plots=transcript_library_plots,
        plot_divs=plot_divs,
        bokeh_script=script,
        )
    with open(args.output, 'wt') as outstream:
        outstream.write(page)


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--quantification', choices=['FPKM', 'TPM'],
                        default='FPKM',
                        help='which quantification value to use')
    parser.add_argument('-l', '--libraries', help='library information table')
    parser.add_argument('-e', '--experiments', help='experiment information table')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    parser.add_argument('-o', '--output', default='report.html',
                        help='output html filename')

    return parser


def make_coverage_plot(coverage, experiments, experiment):
    library_ids = experiments[experiment]
    subset = coverage[library_ids]
    plot = Line(subset,
                title="Coverage for {}".format(experiment),
                xlabel="position quantile (5' to 3')",
                ylabel="Read depth",
                legend=True)
    return plot


def make_distribution_plot(distribution, experiments, experiment):
    library_ids = experiments[experiment]
    subset = distribution.select(lambda x: x in library_ids).unstack()
    subset.index.names = ['class', 'library_id']
    subset.name = 'fraction'
    subset = subset.reset_index()
    return Bar(subset,
               label='library_id',
               values='fraction',
               stack='class',
               title="Distribution for {}".format(experiment),
               legend=True)


def make_spikein_per_transcript_plot(quantifications, library_id, quantification='FPKM'):
    spikein_cpc = pandas.DataFrame(models.get_single_spike_cpc(), columns=['copies'])
    library = quantifications[library_id]
    spikes = library[library.index.isin(spikein_cpc.index)]
    copies_vs_quant = pandas.concat([spikein_cpc, spikes], axis=1)

    source = ColumnDataSource(
        data=dict(
            x = copies_vs_quant['copies'],
            y = copies_vs_quant[library_id],
            desc = copies_vs_quant.index,
        )
    )
    hover = HoverTool(
        tooltips=[
            ('index', '$index'),
            ('transcripts', '$x'),
            ('fpkms', '$y'),
            ('desc', '@desc'),
        ]
    )
    plot = figure(
        title="Spike-in for library {}".format(library_id),
        tools=["pan","wheel_zoom","reset","resize",hover]
    )
    plot.xaxis.axis_label = "transcripts spiked in"
    plot.yaxis.axis_label = "RNA-Seq RSEM ({})".format(quantification)
    plot.circle('x', 'y', source=source)
    return plot

def make_spikein_variance_plot(quantifications, experiments, experiment, quantification='FPKM'):
    spikein_cpc = models.get_single_spike_cpc().sort_values(inplace=False)
    library_ids = experiments[experiment]
    libraries = quantifications[library_ids]
    spikes = libraries[libraries.index.isin(spikein_cpc.index)]
    spikes = spikes.stack()
    spikes.name = quantification
    spikes.index.names = ['spike-ins', 'library']
    spikes = spikes.reset_index(level=1)
    spikes = spikes.join(spikein_cpc)
    spikes_sorted = spikes.sort_values('cpc')
    spikes_sorted.index.name = 'spike-ins'
    spikes_sorted = spikes_sorted.reset_index()

    # Thanks to
    # https://github.com/bokeh/bokeh/issues/2924#issuecomment-166969014
    # for how to get around the auto sorting of the category with bokeh 0.10+
    plot = BoxPlot(
        spikes_sorted,
        values=quantification,
        label=CatAttr(columns=['spike-ins'], sort=False),
        color = 'cpc',
        title="Spike-in variance for experiment {}".format(experiment),
        xlabel="spike-in",
        ylabel="RNA-Seq RSEM ({})".format(quantification),
        width=900,
    )
    return plot


def make_correlation_heatmap(scores, score_name, experiment_name, vmin=None, vmax=None, cmap="coolwarm"):
    """Try to intellgently format our heatmap.
    """
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

if __name__ == "__main__":
    main()
