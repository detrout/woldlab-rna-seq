#!/usr/bin/python3
from __future__ import print_function, absolute_import
import argparse
import os
import logging
from urllib.parse import urljoin

import trackhub

from woldrnaseq.common import (
    add_trailing_slash,
    add_debug_arguments,
    add_version_argument,
    configure_logging,
    get_seperator,
)
from woldrnaseq.models import (
    genome_name_from_library,
    load_experiments,
    load_library_tables,
)

logger = logging.getLogger('make_tracks')


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    if args.email is None:
        parser.error("Please either set email with --email or EMAIL environment variable")

    if args.long_name is None:
        args.long_name = args.short_name

    sep = get_seperator(args.sep)

    args.web_root = add_trailing_slash(args.web_root)

    library_filenames = args.libraries
    if len(library_filenames) == 0:
        parser.error('Need library information table')
    libraries = load_library_tables(library_filenames, sep)

    experiment_filenames = args.experiments
    if len(experiment_filenames) == 0:
        parser.error('Need experiment information table')
    experiments = load_experiments(experiment_filenames, sep)

    genomes = set(libraries['genome'])
    if len(genomes) > 1:
        raise ValueError('We can only generate tracks for one genome')
    genome = ucsc_genome_conversion(genomes.pop())

    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name=args.hub,
        short_label=args.short_name,
        long_label=args.long_name,
        genome=genome,
        email=args.email)

    tracks_added = False
    if args.bigwig:
        make_bigwig_trackhub(experiments, libraries, trackdb, args.web_root)
        tracks_added = True
    if args.bam:
        make_bam_trackhub(experiments, libraries, trackdb, args.web_root)

    if not tracks_added:
        print("Did you want to add tracks? use --bigwig and/or --bam")

    if not args.dry_run:
        trackhub.upload.upload_hub(
            hub=hub,
            host='localhost',
            remote_dir=args.output)
        print('trackhub: ' + args.web_root + hub.hub + '.hub.txt')
    else:
        print(trackdb)

    #print(trackdb)
    #print(hub)
    #print(genomes_file)
    #print(genome)
    #print('trackhub: ' +  args.web_root + hub.hub + '.hub.txt')


def make_parser():
    parser = argparse.ArgumentParser(
        description='Generate a track hub of bigwigs or bam files for '
                    'visualization with the UCSC Genome Browser'
    )
    email = os.environ.get('EMAIL', None)
    parser.add_argument('--hub', help='hub name', required=True)
    parser.add_argument('-n', '--short-name', required=True)
    parser.add_argument('--long-name')
    parser.add_argument('--email', default=email)
    parser.add_argument('-w', '--web-root', required=True,
                        help='base URL for the track hub. ')
    parser.add_argument('-o', '--output', default='./',
                        help='base directory to write the track hub to')
    parser.add_argument('--bigwig', action='store_true', default=False,
                        help='generate track blocks for bigwigs')
    parser.add_argument('--bam', action='store_true', default=False,
                        help='generate track blocks for bam files')

    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    parser.add_argument('-l', '--libraries', action='append', default=[])
    parser.add_argument('-e', '--experiments', action='append', default=[])
    parser.add_argument('--dry-run', action='store_true', default=False)
    add_version_argument(parser)
    add_debug_arguments(parser)

    return parser


def make_bigwig_trackhub(experiments, libraries, trackdb, baseurl):
    experiment_mapping = {}
    for key in experiments.index:
        experiment_mapping[key] = key

    experiment_group = trackhub.SubGroupDefinition(
            name='experiment',
            label='Experiment',
            mapping=experiment_mapping)

    subgroups = [
        experiment_group,
        trackhub.SubGroupDefinition(
            name='multi',
            label='multi',
            mapping={
                'uniq': 'uniq',
                'all': 'all',
            }),
    ]

    composite = trackhub.CompositeTrack(
        name='composite',
        short_label='signal',
        dimensions='dimX=experiment dimY=multi',
        tracktype='bigWig',
        visibility='dense',
    )
    composite.add_subgroups(subgroups)
    trackdb.add_tracks(composite)

    signal_view = trackhub.ViewTrack(
        name='signalviewtrack',
        view='signal',
        visibility='dense',
        tracktype='bigWig',
        short_label='Signal')
    composite.add_view(signal_view)

    priority = 0
    for experiment_name, experiment in experiments.iterrows():
        for library_id in experiment['replicates']:
            row = libraries.loc[library_id]
            for track_type in ['uniq', 'all']:
                track = trackhub.Track(
                    url=make_bigwig_url(baseurl, row, track_type),
                    name="{:02d}".format(priority) + '_' + row.analysis_name + '_' + track_type,
                    visibility='dense',
                    tracktype='bigWig',
                    subgroups={
                        'experiment': experiment_name,
                        'multi': track_type},
                    priority=priority,
                )
                signal_view.add_tracks(track)
                priority += 1


def make_bigwig_url(baseurl, row, track_type):
    analysis_name = add_trailing_slash(row.analysis_name)
    url = urljoin(baseurl, analysis_name)
    url = urljoin(url,
                  row.analysis_name + '-' +
                  genome_name_from_library(row) + '_' +
                  track_type + '.bw')
    return url


def make_bam_custom_track(library, web_root, analysis_root):
    """Generate a bigwig custom track record

    this is singular 'track' because it returns a single row

    :param Series library: row from a library table DataFrame
    :param str web_root: base url to be prepended to paths
    :param str analysis_root: base directory to look for track files
    """
    track_template = 'track type=bam name={library_id} description={description} visibility=dense db={genome} bigDataUrl={url}'

    pathname = make_bam_track_name(library, analysis_root)
    url = web_root + pathname.replace(analysis_root, '')
    track = track_template.format(library_id=library.name,
                                  description=library.analysis_name,
                                  url=url,
                                  genome=library.genome,
    )

    return track


def make_bam_track_name(library, analysis_root=None):
    """Generate the base path where the bam track is.

    :param Series library: row from a library table DataFrame
    :param str analysis_root: root directory to be searching for track files
    :returns: path of bam file relative to analysis_root
    """
    genome_triplet = genome_name_from_library(library)
    track_name = library.analysis_name + '-' + genome_triplet + '_genome.bam'
    old_name = 'Aligned.sortedByCoord.out.bam'
    to_check = [
        os.path.join(library.analysis_dir, track_name),
        os.path.join(analysis_root, track_name),
        os.path.join(library.analysis_dir, old_name),
    ]
    for pathname in to_check:
        if os.path.exists(pathname):
            bai = pathname + '.bai'
            if not os.path.exists(bai):
                logger.warning('Missing index file for {}'.format(pathname))
            return return_subpath(pathname, analysis_root)

    logger.warning("Couldn't find track file %s", track_name)


def make_bigwig_custom_tracks(library, web_root, analysis_root):
    """Generate a bigwig custom track record

    this is plural 'tracks' because it returns uniq and all bigwig tracks

    :param Series library: row from a library table DataFrame
    :param str web_root: base url to be prepended to paths
    :param str analysis_root: base directory to look for track files
    """
    track_template = 'track type=bigWig name={library_id} description={description} visibility=full color=255,0,0 db={genome} bigDataUrl={url}'

    tracks = []
    for signal_type in ['uniq', 'all']:
        pathname = make_bigwig_track_name(library, signal_type, analysis_root)
        url = web_root + pathname.replace(analysis_root, '')
        track = track_template.format(library_id=library.name,
                                      description=library.analysis_name,
                                      url=url,
                                      genome=library.genome,
        )

        tracks.append(track)
    return tracks


def make_bigwig_track_name(library, signal_type, analysis_root):
    """Generate the base path where the bigwig track is

    :param Series library: row from a library table DataFrame
    :param str signal_type: either uniq or all to specify bigwig type.
    :param str analysis_root: root directory to be searching for track files
    :returns: list of paths of bigWig files relative to analysis_root
    """
    assert signal_type in ('uniq', 'all')

    genome_triplet = genome_name_from_library(library)
    track_name = library.analysis_name + '-' + genome_triplet + '_' + signal_type + '.bw'

    for pathname in [os.path.join(library.analysis_dir, track_name),
                     os.path.join(analysis_root, track_name)]:
        if os.path.exists(pathname):
            return return_subpath(pathname, analysis_root)

    logger.warning("Couldn't find track file %s", track_name)


def return_subpath(pathname, analysis_root):
    """Strip off analysis_root from path to file

    :param str pathname: absolute path to file of interest
    :param str analysis_root: root directory to be searching for track files
    :returns: relative path rooted at analysis_root
    """
    if pathname.startswith(analysis_root):
        common = os.path.commonpath([pathname, analysis_root])
        assert common[-1] != '/'
        return pathname.replace(common + '/', '')
    else:
        raise ValueError("Path {} does not start with {}".format(pathname, analysis_root))


def ucsc_genome_conversion(genome_name):
    """Convert genome names to UCSC convention

    :Parameters:
        genome_name: str our genome names
    :Returns:
        str UCSC genome name
    """
    conversions = {
        'GRCh38': 'hg38',
    }
    if genome_name in conversions:
        return conversions[genome_name]
    else:
        return genome_name


if __name__ == '__main__':
    main()
