#!/usr/bin/python3
from __future__ import absolute_import, print_function

import argparse
import datetime
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
import bokeh
from bokeh.models import HoverTool
from bokeh.plotting import figure, ColumnDataSource
from bokeh.embed import components

from .models import (
    load_experiments,
    load_library_tables,
    load_correlations,
    load_all_coverage,
    load_all_distribution,
    load_all_samstats,
    load_quantifications,
    get_single_spike_cpc,
    genome_name_from_library,
)
from . import madqc
from .version import get_git_version

from .common import (
    add_default_path_arguments,
    add_debug_arguments,
    add_version_argument,
    configure_logging,
    get_seperator,
)
from .iplots.correlation import ScoreCorrelationPlot
from .iplots.genes_detected import GenesDetectedPlot
from .iplots.gene_coverage import GeneCoverage
from .iplots.distribution import DistributionPlot

logger = logging.getLogger('QC Report')

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    if args.version:
        parser.exit(0, 'Version: %s' % (get_git_version(),))

    configure_logging(args)

    sep = get_seperator(args.sep)

    if not args.experiments:
        parser.error("Please provide experiment table filename")

    if not args.libraries:
        parser.error("Please provide library table filename")

    report = QCReport(args.experiments, args.libraries,
                      quantification=args.quantification,
                      genome_dir=args.genome_dir,
                      sep=sep,
                      root=args.root,
    )

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
    parser.add_argument('--root', default=None,
                        help='analysis_dir will be relative to this path '\
                        'instead of library.txt file')
    add_default_path_arguments(parser)
    add_version_argument(parser)
    add_debug_arguments(parser)

    return parser


class QCReport:
    def __init__(self, experiments, libraries, quantification, genome_dir,
                 sep="\t",
                 root=None
    ):
        # user set parameters
        self.experiments = load_experiments(experiments, sep, analysis_root=root)
        self.libraries = load_library_tables(libraries, sep, analysis_root=root)
        self.quantification_name = quantification
        self.genome_dir = genome_dir
        self.analysis_root = root

        # cached values
        self._samstats = load_all_samstats(self.libraries)

        # working space for generating report
        self._plot_handle = itertools.count()
        self._plots = {}
        self._experiment_report = {}
        self._transcript_library_plots = []

    @property
    def next_plot_handle(self):
        return str(next(self._plot_handle))

    def render(self):
        env = Environment(loader=PackageLoader('woldrnaseq', 'templates'))

        self.generate_report()
        script, plot_divs = components(self._plots)
        template = env.get_template('rnaseq.html')
        page = template.render(
            experiments=self.experiments.index,
            experiment_report=self._experiment_report,
            transcript_library_plots=self._transcript_library_plots,
            plot_divs=plot_divs,
            bokeh_script=script,
            username=os.getlogin(),
            timestamp=datetime.datetime.now().isoformat(),
            woldrnaseq_version=get_git_version(),
            bokeh_version=bokeh.__version__,
            )
        return page

    def generate_report(self):
        seen_libraries = set()

        coverage_plots = GeneCoverage(self.experiments, self.libraries)
        distribution_plots = DistributionPlot(self.experiments, self.libraries)
        genes_detected_plots = GenesDetectedPlot(self.experiments, self.libraries, self.genome_dir)
        score_correlation_plots = ScoreCorrelationPlot(self.experiments)

        for experiment_name in sorted(self.experiments.index):
            experiment = self.experiments.loc[experiment_name]
            quantifications = load_quantifications(
                experiment,
                self.quantification_name)
            scores = load_correlations(experiment)

            if quantifications is None:
                raise FileNotFoundError(
                    "Unable to load quantification {} for {}".format(
                        self.quantification_name,
                        experiment_name))

            library_ids = experiment['replicates']
            seen_libraries.update(set(library_ids))
            seen_genomes = set()
            for library_id in library_ids:
                #transcript_handle = self.make_spikein_per_transcript_plot(
                #    quantifications,
                #    library_id)
                #self._transcript_library_plots.append(transcript_handle)
                genome_name = genome_name_from_library(self.libraries.loc[library_id])
                seen_genomes.add(genome_name)

            if len(seen_genomes) > 1:
                logger.warning('%s has mixed genome types %s',
                               experiment_name, ','.join(seen_genomes))

            cur_experiment = {
                'spearman': scores.rafa_spearman.to_html(na_rep=''),
                'samstats': self.make_samstats_html(library_ids),
                'coverage': self.make_plot_handle(coverage_plots, experiment_name),
                'distribution': self.make_plot_handle(distribution_plots, experiment_name),
                'correlation': self.make_plot_handle(score_correlation_plots, experiment_name),
                'protein_genes': self.make_plot_handle(genes_detected_plots, experiment_name),
            }

            self._experiment_report[experiment_name] = cur_experiment
            #handle = self.make_spikein_variance_plot(quantifications, experiment)
            #if handle:
            #    cur_experiment['spike_variance'] = handle
            #else:
            #    print("Didn't generate spike report for:", experiment)

        spare_libraries = set(self.libraries.index).difference(seen_libraries)

    def make_plot_handle(self, plot, experiment):
        handle = self.next_plot_handle
        logger.debug('%s plot handle %s: %s', type(plot), experiment, handle)
        self._plots[handle] = plot.make_plot(experiment)
        return handle

    def make_samstats_html(self, library_ids):
        def reads_format(x):
            return "{:,}".format(int(x))

        def fraction_mapped(x):
            return "{:.3}".format(x)

        def read_length(x):
            return "{:,}".format(int(x))

        return self._samstats.loc[library_ids].to_html(
            formatters={
                'Unique': reads_format,
                'Unique Splices': reads_format,
                'Multi': reads_format,
                'Multi Splices': reads_format,
                'Total Aligned Reads': reads_format,
                'Fraction Mapped': fraction_mapped,
                'Read Length, Minimum': read_length,
                'Read Length, Maximum': read_length,
            })


    def make_spikein_per_transcript_plot(self, quantifications, library_id):
        """WARNING: Outdated plot

        Make scatter plot of spike-in FPKM vs expected FPKM
        """
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
            tools=["pan","wheel_zoom","reset",hover]
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
        """WARNING: Outdated plot

        Make box plot showing typical performance of the spike-ins
        """
        spikein_cpc = get_single_spike_cpc().sort_values(inplace=False)
        library_ids = self.experiments.loc[experiment]['replicates']
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
        plot = figure(
            x_range=spikes_sorted.index,
            title="Spike-in variance for experiment {}".format(experiment),
        )
        plot = BoxPlot(
            spikes_sorted,
            values=self.quantification_name,
            label=CatAttr(columns=['spike-ins'], sort=False),
            color = 'cpc',
            xlabel="spike-in",
            ylabel="RNA-Seq RSEM ({})".format(self.quantification_name),
            width=900,
        )

        handle = self.next_plot_handle
        self._plots[handle] = plot
        return handle


if __name__ == "__main__":
    main()
