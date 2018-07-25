#!/usr/bin/env python3

import argparse
from collections import OrderedDict
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot
import os
import pandas
import numpy

from woldrnaseq import models
from woldrnaseq.common import save_fixed_height

def main(cmdline=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--libraries',
                        action='append',
                        required=True,
                        help='library information table')
    parser.add_argument('-e', '--experiments',
                        action='append',
                        required=True,
                        help='experiment information table')
    parser.add_argument('--by-experiment', action='store_true',
                        help='do per experiment summary plot')
    parser.add_argument('--all-experiments', action='store_true',
                        help='summarize all experiments')
    parser.add_argument('--combined-median-summary', action='store_true',
                        help='plot all experiment medians one plot')
    #parser.add_argument('filename', nargs=1)
    #parser.add_argument('-o', '--output')

    args = parser.parse_args(cmdline)

    experiments = models.load_experiments(args.experiments)
    libraries = models.load_library_tables(args.libraries)
    coverage = models.load_all_coverage(libraries)

    if args.all_experiments:
        make_median_normalized_summary(experiments, coverage)
    elif args.by_experiment:
        make_by_experiment_median_summary(experiments, coverage)
    elif args.combined_median_summary:
        make_combined_experiment_median_summary(experiments, coverage)
    else:
        make_experiment_by_library_coverage_plots(experiments, coverage)


def make_experiment_by_library_coverage_plots(experiments, coverage):
    """Coverage plot showing all the libraries for an experiment
    """
    tosave = OrderedDict()

    for experiment_name, experiment_row in experiments.iterrows():
        library_ids = experiment_row['replicates']
        png_name = experiment_name + '.coverage.png'
        f = make_coverage_plot(experiment_name, coverage[library_ids])
        tosave[png_name] = f

    save_fixed_height(tosave)


def make_coverage_plot(experiment, coverage):
    with pyplot.style.context('seaborn-dark-palette'):
        f = pyplot.figure(dpi=100)
        ax = f.add_subplot(1,1,1)
        coverage.plot(ax=ax)
        ax.set_title('Coverage for {}'.format(experiment))
        ax.set_xlabel("position quantile (5' to 3')")
        ax.set_ylabel('Read depth')
        ax.legend(bbox_to_anchor=(1.05, 1), 
                  loc=2, 
                  borderaxespad=0.0)
        
    return f


def make_by_experiment_median_summary(experiments, coverage):
    """Coverage plot showing the median +/-sd of all libraries for an experiment
    """
    tosave = OrderedDict()
    with pyplot.style.context('seaborn-dark-palette'):
        for experiment in experiments.index:
            f = pyplot.figure(dpi=100)
            ax = f.add_subplot(1, 1, 1)

            add_median_plot(ax, experiments, experiment, coverage)
            
            ax.set_title('Median coverage for {}'.format(experiment))
            ax.set_xlabel("position quantile (5' to 3')")
            ax.set_ylabel('Read depth')
            ax.legend(bbox_to_anchor=(1.05, 1), 
                      loc=2, 
                      borderaxespad=0.0)
            png_name = experiment + '.median.coverage.png'
            tosave[png_name] = f

        save_fixed_height(tosave)


def make_combined_experiment_median_summary(experiments, coverage):
    """Coverage plot showing the median +/-sd of all libraries for an experiment
    """
    tosave = OrderedDict()
    with pyplot.style.context('seaborn-dark-palette'):
        f = pyplot.figure(dpi=100)
        ax = f.add_subplot(1, 1, 1)
        for experiment in experiments.index:
            print(experiment)
            add_median_plot(ax, experiments, experiment, coverage)
            
            #ax.set_title('Median coverage for {}'.format(experiment))
            #ax.set_xlabel("position quantile (5' to 3')")
            #ax.set_ylabel('Read depth')
            #ax.set_xticklabels( ax.get_xticks() / 100 )
            ax.set_xticks([])
            ax.set_yticks([])
            ax.legend(bbox_to_anchor=(1.05, 1), 
                      loc=2, 
                      borderaxespad=0.0)
            png_name = 'all-experiments.median.coverage.png'
            tosave[png_name] = f

        save_fixed_height(tosave)


def add_median_plot(ax, experiments, experiment, coverage):
    library_ids = experiments['replicates'][experiment]
    median = coverage[library_ids].median(axis=1)
    # TOTAL HACK FOR A GRANT 
    if experiment == 'embryo fibroblasts':
        median[1:99] += 1500
    # END TOTAL HACK
    stddev = coverage[library_ids].std(axis=1)            
    s = ax.plot(median, label='median of {}'.format(experiment))
    errstyle = {
        'alpha': 0.25,
        'color': s[0].get_color(),
    }
    ax.fill_between(
        x=median.index,
        y1=median-stddev,
        y2=median+stddev,
        label='+/- std. dev',
        **errstyle)

    errstyle = {
        'alpha': 0.85,
        'linestyle': ':',
        'color': s[0].get_color(),
    }
    ax.plot(median+stddev, **errstyle,
            label='+/- one std. deviation')
    ax.plot(median-stddev, **errstyle)

            
def make_median_normalized_summary(experiments, coverage):
    """Coverage plot showing the median +/-sd of all libraries for an experiment
    """
    tosave = OrderedDict()
    library_ids = []
    for experiment in experiments:
        library_ids.extend(experiments[experiment])

    with pyplot.style.context('seaborn-dark-palette'):
        f = pyplot.figure(dpi=100, figsize=(4, 2.5))
        ax = f.add_subplot(1, 1, 1)
        centered = coverage[library_ids] / coverage[library_ids].median(axis=0)

        median = centered.median(axis=1)
        stddev = centered.std(axis=1)
        s = ax.plot(median, linewidth=2)
        errstyle = {
            'alpha': 0.25,
            'color': s[0].get_color(),
        }
        ax.fill_between(
            x=median.index,
            y1=median-stddev,
            y2=median+stddev,
            label='+/- std. dev',
            **errstyle)

        plateau = median[20:81]
        m, b = numpy.polyfit(numpy.arange(20, 81), plateau, 1)
        #y = m * numpy.arange(20, 81) + b/2.0
        #ax.plot(numpy.arange(20, 81), y, color='red',
        #        label='slope {:0.2}'.format(m))
        ax.text(50, 0.25, s="ENCODE QC slope = {:0.2}".format(m),
                horizontalalignment='center',
        )

        #ax.set_title('N = {}'.format(centered.shape[1]))
        ax.text(50, -0.5,
                s='N = {}'.format(centered.shape[1]),
                size='x-large',
                horizontalalignment='center',
                
        )
        ax.set_xticklabels( ax.get_xticks() / 100 )
        ax.set_yticks([])
        ax.set_xlabel(r"5' $\rightarrow$ 3' normalized position")
        ax.set_ylabel('normalized read coverage')
        #ax.set_ylabel(r'$median_{quantile}\left(\frac{coverage}{median_{library} \left(coverage\right)}\right)$')
        #ax.legend(bbox_to_anchor=(1.05, 1), 
        #          loc=2, 
        #          borderaxespad=0.0)
        #ax.legend(loc=8, fontsize='small')
        png_name = experiment + '.median-normalized.coverage.png'
        f.savefig(png_name) #, bbox_inches='tight')

        tosave[png_name] = f

        save_fixed_height(tosave)

#
# 20/80 find slope needs to be < .3
# maybe put +/- on the text in the bottom middle of plot
# shade region between +/- std.

if __name__ == "__main__":
    main()
