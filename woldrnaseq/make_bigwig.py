#!/usr/bin/python3
from argparse import ArgumentParser
import logging
import os
import sys
import subprocess
from dask import delayed
from dask.distributed import Client, LocalCluster
from pysam import AlignmentFile

from woldrnaseq.common import (
    add_debug_arguments,
    add_default_path_arguments,
    configure_logging,
)

logger = logging.getLogger(__name__)


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    logger.debug("current directory {}".format(os.getcwd()))
    logger.debug("env: {}".format(os.environ))


    if args.prefix:
        prefix = args.prefix
    else:
        bam_extension = '_genome.bam'
        if args.bam.endswith(bam_extension):
            base = make_basename(args.bam)
            prefix = base[:-len('_genome')]
        else:
            parser.error('Target prefix must be provided for non-standard bam names')

    if args.stranded:
        targets = {
            'Signal.UniqueMultiple.str1.out.bg': prefix + '_minusAll.bw',
            'Signal.Unique.str1.out.bg': prefix + '_minusUniq.bw',
            'Signal.UniqueMultiple.str2.out.bg': prefix + '_plusAll.bw',
            'Signal.Unique.str2.out.bg': prefix + '_plusUniq.bw',
        }
    else:
        targets = {
            'Signal.UniqueMultiple.str1.out.bg': prefix + '_all.bw',
            'Signal.Unique.str1.out.bg': prefix + '_uniq.bw',
        }

    local_cluster = LocalCluster(n_workers=1, threads_per_worker=len(targets))
    client = Client(local_cluster)
    star_p = client.submit(
        run_star_to_bedgraph,
        args.bam,
        args.stranded,
        args.reference_prefix,
        args.star_dir,
    )
    star_p.result()

    chrom_info = delayed(make_chrom_info)(args.bam)
    bigwigs = [star_p]
    for target in targets:
        bedsort = delayed(run_bedsort)(target)
        bigwig = delayed(run_bedgraph2bigwig)(bedsort, chrom_info, targets[target], args.ucsc_tools_dir)
        bigwigs.append(bigwig)

    results = client.compute(bigwigs)
    status, generated = client.gather(results)
    print('Generated', ','.join([str(x) for x in generated]))
    client.close()
    return 0


def make_parser():
    parser = ArgumentParser()
    parser.add_argument('--bam', required=True, help='Source bam file')
    parser.add_argument('--prefix', help='set output name prefix')
    parser.add_argument('--stranded', action='store_true', default=False,
                        help='Process as stranded data')
    parser.add_argument('--reference-prefix', default='chr')
    add_default_path_arguments(parser)
    add_debug_arguments(parser)
    return parser


def run_star_to_bedgraph(bam, stranded, reference_prefix, star_dir=None):
    if stranded:
        stranded_option = 'Stranded'
    else:
        stranded_option = 'Unstranded'

    star = 'STAR' if star_dir is None else os.path.join(star_dir, 'STAR')
    cmd = [star,
           '--runMode', 'inputAlignmentsFromBAM',
           '--inputBAMfile', bam,
           '--outWigType', 'bedGraph',
           '--outWigReferencesPrefix', reference_prefix,
           '--outWigStrand', stranded_option]

    print('running', ' '.join(cmd))
    subprocess.run(cmd, check=True)
    return True


def run_bedsort(filename):
    cmd = ['bedSort', filename, filename]
    print('running', ' '.join(cmd))
    subprocess.run(cmd, check=True)
    return filename


def run_bedgraph2bigwig(source, chrom_info, destination, ucsc_tools_dir=None):
    bedgraph_cmd = 'bedGraphToBigWig' if ucsc_tools_dir is None else os.path.join(ucsc_tools_dir, 'bedGraphToBigWig')
    cmd = [bedgraph_cmd, source, chrom_info, destination]
    print('running', ' '.join(cmd))
    subprocess.run(cmd, check=True)
    os.unlink(source)
    return destination


def make_chrom_info(bam):
    base = make_basename(bam)
    chrom_info_filename = base + '.chrom_info'
    with AlignmentFile(bam, 'rb') as alignment, open(chrom_info_filename, 'wt') as chrom_info:
        for row in alignment.header['SQ']:
            name = row['SN']
            length = row['LN']
            chrom_info.write(name)
            chrom_info.write('\t')
            chrom_info.write(str(length))
            chrom_info.write(os.linesep)

    return chrom_info_filename


def make_basename(bam):
    path, name = os.path.split(bam)
    base, ext = os.path.splitext(name)
    return base


if __name__ == '__main__':
    sys.exit(main())
