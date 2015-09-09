#!/usr/bin/python3
import argparse
from jinja2 import Environment, PackageLoader
from bokeh import mpl
from bokeh.charts import Line
from bokeh.resources import CDN
from bokeh.embed import components

import models

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)
    env = Environment(loader=PackageLoader('woldrnaseq', 'templates'))
    sep = models.get_seperator(args.sep)

    experiments = models.load_experiments(['experiments.tsv'], sep)
    scores = models.load_all_correlations(experiments)

    libraries = models.load_library_tables(['libraries.tsv'], sep)
    samstats = models.load_all_samstats(libraries)
    distribution = models.load_all_distribution(libraries)
    coverage = models.load_all_coverage(libraries)

    spare_libraries = set(libraries.index).difference(set(experiments))

    experiment_report = []
    scripts = []
    for experiment, library_ids in experiments.items():
        print(scores[experiment].keys())
        script, coverage_div = components(make_coverage_plot(coverage, library_ids))
        experiment_report.append((
            experiment,
            samstats.select(lambda x: x in library_ids).to_html(),
            scores[experiment].rafa_spearman.to_html(),
            coverage_div
        ))
        scripts.append(script)

    template = env.get_template('rnaseq.html')
    print(template.render(
        experiment_report=experiment_report,
        scripts=scripts,
        ))


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--library', help='library information table')
    parser.add_argument('-e', '--experiments', help='experiment information table')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')

    return parser


def make_coverage_plot(coverage, library_ids):
    subset = coverage #coverage.select(lambda x: x in library_ids)
    subset.to_csv('/tmp/subset.csv')
    print(subset)
    return Line(subset,
                title="Coverage",
                xlabel="position percentile",
                ylabel="Read depth",
                legend=True)
if __name__ == "__main__":
    main()
