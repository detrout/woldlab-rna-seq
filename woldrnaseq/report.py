#!/usr/bin/python3
from __future__ import absolute_import, print_function

import argparse
import itertools
import logging
import os

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
from bokeh.charts import Line, Bar, Dot, HeatMap, BoxPlot, Histogram
from bokeh.charts.attributes import CatAttr
from bokeh.embed import components

from woldrnaseq import models
from .models import (load_experiments,
                     load_library_tables,
                     load_correlations,
                     load_all_coverage,
                     load_all_distribution,
                     load_all_samstats,
                     load_quantifications,
                     get_single_spike_cpc,
                     genome_name_from_library,
)
from woldrnaseq import madqc
from .common import (add_default_path_arguments,
                     add_debug_arguments,
                     configure_logging,
                     get_seperator,
)
from .plot_genes_detected import (load_gtf_cache,
                                  bin_library_quantification,
                                  protein_coding_gene_ids,
                                  plot_gene_detection_histogram,
)

logger = logging.getLogger('QC Report')

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    sep = get_seperator(args.sep)

    if not args.experiments:
        parser.error("Please provide experiment table filename")

    if not args.libraries:
        parser.error("Please provide library table filename")

    report = QCReport(args.experiments, args.libraries,
                      quantification=args.quantification,
                      genome_dir=args.genome_dir,
                      sep=sep)

    with open(args.output, 'wt') as outstream:
        outstream.write(report.render())


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--quantification', choices=['FPKM', 'TPM'],
                        default='FPKM',
                        help='which quantification value to use')
    parser.add_argument('-l', '--libraries', action='append',
                        help='library information table')
    parser.add_argument('-e', '--experiments', action='append',
                        help='experiment information table')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    parser.add_argument('-o', '--output', default='report.html',
                        help='output html filename')
    add_default_path_arguments(parser)
    add_debug_arguments(parser)

    return parser


class QCReport:
    def __init__(self, experiments, libraries, quantification, genome_dir, sep="\t"):
        # user set parameters
        self.experiments = load_experiments(experiments, sep)
        self.libraries = load_library_tables(libraries, sep)
        self.quantification_name = quantification
        self.genome_dir = genome_dir

        # cached values
        self._samstats = load_all_samstats(self.libraries)
        self._distribution = load_all_distribution(self.libraries)
        self._coverage = load_all_coverage(self.libraries)
        self._gtf_cache = {}

        # working space for generating report
        self._plot_handle = itertools.count()
        self._plots = {}
        self._experiment_report = {}
        self._transcript_library_plots = []

    def _load_gtf_cache(self, genome_name):
        """Look through list of libraries and attempt to load GTF caches
        """
        if genome_name not in self._gtf_cache:
            cache_pathname = os.path.join(self.genome_dir, genome_name, genome_name + '.h5')
            if os.path.exists(cache_pathname):
                load_gtf_cache(cache_pathname)
            else:
                logging.error('Unable to load gene cache %s', cache_pathname)

    @property
    def next_plot_handle(self):
        return str(next(self._plot_handle))

    def render(self):
        env = Environment(loader=PackageLoader('woldrnaseq', 'templates'))

        self.generate_report()
        script, plot_divs = components(self._plots)
        template = env.get_template('rnaseq.html')
        page = template.render(
            experiments=self.experiments,
            experiment_report=self._experiment_report,
            transcript_library_plots=self._transcript_library_plots,
            plot_divs=plot_divs,
            bokeh_script=script,
            )
        return page

    def generate_report(self):
        seen_libraries = set()
        for experiment in sorted(self.experiments):
            quantifications = load_quantifications(
                experiment,
                self.quantification_name)

            library_ids = self.experiments[experiment]
            seen_libraries.update(set(library_ids))
            seen_genomes = set()
            for library_id in library_ids:
                transcript_handle = self.make_spikein_per_transcript_plot(
                    quantifications,
                    library_id)
                self._transcript_library_plots.append(transcript_handle)
                genome_name = genome_name_from_library(self.libraries.loc[library_id])
                self._load_gtf_cache(genome_name)
                seen_genomes.add(genome_name)

            if len(seen_genomes) > 1:
                logger.warning('%s has mixed genome types %s',
                               experiment, ','.join(seen_genomes))

            cur_experiment = {
                'samstats': self.make_samstats_html(library_ids),
                'coverage': self.make_coverage_plot(experiment),
                'distribution': self.make_distribution_plot(experiment),
            }

            if genome_name in self._gtf_cache:
                cur_experiment['protein_genes'] = self.plot_genes_detected(
                    quantifications,
                    genome_name,
                    experiment,
                )

            handle = self.make_spikein_variance_plot(quantifications, experiment)
            if handle:
                cur_experiment['spike_variance'] = handle
            else:
                print("Didn't generate spike report for:", experiment)

            cur_experiment.update(self.make_correlation_plots(experiment))
            self._experiment_report[experiment] = cur_experiment

        spare_libraries = set(self.libraries.index).difference(seen_libraries)

    def make_correlation_plots(self, experiment):
        report = {}
        scores = load_correlations(experiment)
        if scores.shape[0] > 0:
            report['spearman_plot'] = make_correlation_heatmap(
                scores, 'rafa_spearman', experiment)
            report['spearman'] = scores.rafa_spearman.to_html()
            if scores.shape[1] > 2:
                handle = self.next_plot_handle
                self._plots[handle] = make_correlation_histogram(
                    scores, 'rafa_spearman', experiment)
                report['correlation_histogram'] = handle
        else:
            print('scores failed')

        return report

    def make_coverage_plot(self, experiment):
        """Show read depth coverage over normalized gene regions.
        """
        library_ids = self.experiments[experiment]
        subset = self._coverage[library_ids]
        plot = Line(subset,
                    title="Coverage for {}".format(experiment),
                    xlabel="position quantile (5' to 3')",
                    ylabel="Read depth",
                    legend=True)
        handle = self.next_plot_handle
        logger.debug('coverage plot handle for %s: %s', experiment, handle)
        self._plots[handle] = plot
        return handle

    def make_distribution_plot(self, experiment):
        """Show fraction of reads landing in exon, intron, and intergenic regions.
        """
        library_ids = self.experiments[experiment]
        subset = self._distribution.select(lambda x: x in library_ids).unstack()
        subset.index.names = ['class', 'library_id']
        subset.name = 'fraction'
        subset = subset.reset_index()
        plot = Bar(subset,
                   label='library_id',
                   values='fraction',
                   stack='class',
                   title="Distribution for {}".format(experiment),
                   legend=True,
        )
        handle = self.next_plot_handle
        logger.debug('distribution plot handle for %s: %s', experiment, handle)
        self._plots[handle] = plot
        return handle

    def make_samstats_html(self, library_ids):
        return self._samstats.select(lambda x: x in library_ids).to_html()

    def make_spikein_per_transcript_plot(self, quantifications, library_id):
        spikein_cpc = pandas.DataFrame(get_single_spike_cpc(),
                                       columns=['cpc'])
        library = quantifications[library_id]
        spikes = library[library.index.isin(spikein_cpc.index)]
        if len(spikes) == 0:
            logger.warning("No spike-ins detected for %s", library_id)
            return None

        copies_vs_quant = pandas.concat([spikein_cpc, spikes], axis=1)

        source = ColumnDataSource(
            data=dict(
                copies = copies_vs_quant['cpc'],
                fpkms = copies_vs_quant[library_id],
                desc = copies_vs_quant.index,
            )
        )
        hover = HoverTool(
            tooltips=[
                ('index', '$index'),
                ('transcripts', '@copies'),
                ('fpkms', '@fpkms'),
                ('desc', '@desc'),
            ]
        )
        plot = figure(
            title="Spike-in for library {}".format(library_id),
            tools=["pan","wheel_zoom","reset","resize",hover]
        )
        plot.xaxis.axis_label = "transcripts spiked in"
        plot.yaxis.axis_label = "RNA-Seq RSEM ({})".format(self.quantification_name)
        #plot.circle('copies', 'fpkms', source=source)
        plot.circle(x=copies_vs_quant['cpc'], y=copies_vs_quant[library_id])

        handle = self.next_plot_handle
        logger.debug('spikein per transcript plot handle for %s: %s',
                     library_id,
                     handle)
        self._plots[handle] = plot
        return handle

    def make_spikein_variance_plot(self, quantifications, experiment):
        spikein_cpc = models.get_single_spike_cpc().sort_values(inplace=False)
        library_ids = self.experiments[experiment]
        libraries = quantifications[library_ids]
        spikes = libraries[libraries.index.isin(spikein_cpc.index)]
        if len(spikes) == 0:
            logger.warning("No spikes detected for %s", str(experiment))
            return None

        spikes = spikes.stack()
        spikes.name = self.quantification_name
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
            values=self.quantification_name,
            label=CatAttr(columns=['spike-ins'], sort=False),
            color = 'cpc',
            title="Spike-in variance for experiment {}".format(experiment),
            xlabel="spike-in",
            ylabel="RNA-Seq RSEM ({})".format(self.quantification_name),
            width=900,
        )

        handle = self.next_plot_handle
        self._plots[handle] = plot
        return handle

    def plot_genes_detected(self, quantifications, genome_name, experiment_name):
        annotation = self._gtf_cache[genome_name]
        protein_coding = protein_coding_gene_ids(annotation)
        png_name = experiment_name + '-genes-detected.png'
        csv_name = experiment_name + '-genes-detected.csv'

        protein_quantifications = quantifications[protein_coding]
        binned_quantifications = bin_library_quantification(
            protein_quantifications,
            self.quantification_name)
        binned_quantifications.to_csv(csv_name)
        f = plot_gene_detection_histogram(binned_quantifications, experiment_name)
        f.savefig(png_name, bbox_inches='tight')

        return png_name

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


def make_correlation_histogram(scores, score_name, experiment_name):
    scores = score_upper_triangular(scores[score_name])
    plot = Histogram(scores, bins=min(10, len(scores)),
                     title='{} scores for {}'.format(
                     score_name, experiment_name))
    return plot


def score_upper_triangular(df):
    """Return the cells from the upper triangular indicies of a data frame.
    """
    scores = []
    for i,j in zip(*numpy.triu_indices(len(df), k=1)):
        scores.append(df.ix[i,j])
    return scores


if __name__ == "__main__":
    main()
