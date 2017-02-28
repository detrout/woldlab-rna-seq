#!/usr/bin/env python3
from __future__ import print_function, unicode_literals, division

import argparse
from collections import OrderedDict
import pandas
import numpy
import os
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot

from woldrnaseq.common import save_fixed_height
from woldrnaseq.models import load_gtf_cache

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    annotation = load_gtf_cache(args.gtf_cache)
    protein_coding = protein_coding_gene_ids(annotation)

    toplot = OrderedDict()
    for filename in args.filenames:
        all_quantifications = pandas.read_csv(filename, header=0, index_col=0)
        if 'gene_name' in all_quantifications.columns:
            columns = [ c for c in all_quantifications.columns if c != 'gene_name']
            all_quantifications = all_quantifications[columns]
        protein_quantifications = all_quantifications.loc[protein_coding]
        
        _, filename = os.path.split(filename)
        basename, _ = os.path.splitext(filename)
        png_name = basename + '.png'
        csv_name = 'genes-detected_' + basename + '.csv'

        binned_quantifications = bin_library_quantification(protein_quantifications, 'FPKM')
        binned_quantifications.to_csv(csv_name)

        f = plot_gene_detection_histogram(binned_quantifications,
                                          basename,
                                          show_genes_detected=not args.hide_detected_sum)
        toplot[png_name] = f

    save_fixed_height(toplot)


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gtf-cache', required=True, help='name of HDF5 GTF file')
    parser.add_argument('--hide-detected-sum', default=False,
                        help='hide the total genes detected')
    parser.add_argument('filenames', nargs='+',
                        help='Combined quantification file: libraries by genes')
    return parser


def protein_coding_gene_ids(annotation):
    """Filter GTF just protein coding genes
    """
    entry_type = (annotation['type'] == 'gene')
    gene_type = (annotation['gene_type'] == 'protein_coding')
    return annotation[entry_type & gene_type]['gene_id']


def bin_library_quantification(quantification, quantification_name, bins=None):
    """Bin a library quantification file

    default bins are [0.1, 1, 2, 5, 10, 50, 500, 5000, 1e9]
    """
    if bins is None:
        bins = [0.1, 1, 2, 5, 10, 50, 500, 5000, 1e9]
    
    histogram = {}
    for col in quantification:
        histogram[col], _ = numpy.histogram(quantification[col], bins)

    histogram = pandas.DataFrame(
        histogram, 
        columns=quantification.columns,
        index=['{} {}'.format(x, quantification_name) for x in bins[:-1]])

    return histogram.reindex(histogram.index[::-1]).T


def plot_gene_detection_histogram(binned_quantifications, basename,
                                  show_genes_detected=True):
    """Apply formatting to gene detction histogram
    """
    with pyplot.style.context('seaborn-talk'):
        width = max(len(binned_quantifications.index) * 0.5, 6)
        f = pyplot.figure(figsize=(width, 6), dpi=100)
        ax = f.add_subplot(1,1,1)

        matplotlib.rcParams['patch.force_edgecolor'] = True
        matplotlib.rcParams['patch.facecolor'] = 'b'
        gene_detection_histogram(ax, binned_quantifications,
                                 show_genes_detected=show_genes_detected)
        ax.set_title(basename)

    return f


def gene_detection_histogram(ax, binned,
                             cm=pyplot.cm.OrRd_r,
                             show_genes_detected=True):
    ax.set_ylabel('Number of genes')
    binned.plot.bar(
        stacked=True, 
        cmap=cm,
        ax=ax)

    if show_genes_detected:
        for rect, total in zip(ax.patches, list(binned.sum(axis=1))):
            x = rect.get_x() + rect.get_width()/2
            y = total + 5
            label = "{:.2}k".format(total/1000)
            ax.text(x, y, label, ha='center', va='bottom')

    ax.legend(bbox_to_anchor=(1.05, 1), 
              loc=2, 
              borderaxespad=0.0)
    return ax

if __name__ == '__main__':
    main()
