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
    sanitize_name,
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

    if args.short_name is None:
        args.short_name = args.hub
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

    hub = build_trackhub(experiments, libraries, args)

    if args.dry_run:
        print(trackhub)
    else:
        ucsc_name_map = build_ucsc_name_map_for_libraries(libraries)
        publish_trackhub(hub, ucsc_name_map, args)


def build_trackhub(experiments, libraries, args):
    ucsc_name_map = build_ucsc_name_map_for_libraries(libraries)

    hub = trackhub.Hub(
        hub=args.hub,
        short_label=args.short_name,
        long_label=args.long_name,
        email=args.email)

    genome_column = detect_genome_column_name(libraries)
    genomes = {}
    genome_files = {}
    trackdbs = {}
    tracks_added = False

    for genome_name in ucsc_name_map:
        ucsc_name = ucsc_name_map[genome_name]
        genomes[ucsc_name] = trackhub.Genome(ucsc_name)
        genome_files[genome_name] = trackhub.GenomesFile()
        trackdbs[genome_name] = trackhub.TrackDb()

        hub.add_genomes_file(genome_files[genome_name])
        genome_files[genome_name].add_genome(genomes[ucsc_name])
        genomes[ucsc_name].add_trackdb(trackdbs[genome_name])

        current_libraries = libraries[libraries[genome_column] == genome_name]
        current_experiment_filter = []
        for replicates in experiments.replicates:
            current = True
            for replicate in replicates:
                if replicate not in current_libraries.index:
                    current = False
            current_experiment_filter.append(current)
        current_experiments = experiments[current_experiment_filter]
        if args.bigwig:
            make_bigwig_trackhub(
                current_experiments,
                current_libraries,
                trackdbs[genome_name],
                args.web_root,
                args.stranded)
            tracks_added = True

        if args.bam:
            make_bam_trackhub(
                current_experiments,
                current_libraries,
                trackdbs[genome_name],
                args.web_root)

    if not tracks_added:
        print("Did you want to add tracks? use --bigwig and/or --bam")

    return hub


def publish_trackhub(hub, ucsc_name_map, args):
    trackhub.upload.upload_hub(
        hub=hub,
        host='localhost',
        remote_dir=args.output)
    hub_url = args.web_root + hub.hub + '.hub.txt'

    print('trackhub: {}'.format(hub_url))
    clickable = "clickable: "\
        "http://genome.ucsc.edu/cgi-bin/hgTracks?db={}&hubUrl={}"
    for genome_name in ucsc_name_map:
        ucsc_name = ucsc_name_map[genome_name]
        print(clickable.format(ucsc_name, hub_url))


def make_parser():
    parser = argparse.ArgumentParser(
        description='Generate a track hub of bigwigs or bam files for '
                    'visualization with the UCSC Genome Browser'
    )
    email = os.environ.get('EMAIL', None)
    parser.add_argument('--hub', required=True,
                        help="Name of hub, used as prefix for filenames")
    parser.add_argument('-n', '--short-name',
                        help="Short name for hub, defaults to --hub")
    parser.add_argument('--long-name', help="long name for hub, defaults to short-name")
    parser.add_argument('--email', default=email,
                        help="Specify responsible person for hub, "
                             "defaults to EMAIL environment variable")
    parser.add_argument('-w', '--web-root', required=True,
                        help='base URL for the track hub. ')
    parser.add_argument('-o', '--output', default='./',
                        help='base directory to write the track hub to')
    parser.add_argument('--bigwig', action='store_true', default=False,
                        help='generate track blocks for bigwigs')
    parser.add_argument('--bam', action='store_true', default=False,
                        help='generate track blocks for bam files')
    parser.add_argument('--stranded', action='store_true', default='False',
                        help='hack to force generation of stranded tracks')
    parser.add_argument('-s', '--sep', choices=['TAB', ','], default='TAB')
    parser.add_argument('-l', '--libraries', action='append', default=[],
                        help="list of library metadata tables to generate tracks for")
    parser.add_argument('-e', '--experiments', action='append', default=[],
                        help="list of experiment metadata tables to generate tracks for")
    parser.add_argument('--dry-run', action='store_true', default=False,
                        help="avoid uploading trackhub files and instead just print them")
    add_version_argument(parser)
    add_debug_arguments(parser)

    return parser


def make_bigwig_trackhub(experiments, libraries, trackdb, baseurl, stranded=False):
    experiment_mapping = {}
    for key in experiments.index:
        experiment_mapping[sanitize_name(key)] = sanitize_name(key)

    experiment_group = trackhub.SubGroupDefinition(
            name='experiment',
            label='Experiment',
            mapping=experiment_mapping)

    if stranded:
        mapping = {
            'minusUniq': 'minusUniq',
            'minusAll': 'minusAll',
            'plusUniq': 'plusUniq',
            'plusAll': 'plusAll',
        }
        track_types = ['minusUniq', 'plusUniq', 'minusAll', 'plusAll']
    else:
        mapping = {
            'uniq': 'uniq',
            'all': 'all',
        }
        track_types = ['uniq', 'all']

    subgroups = [
        experiment_group,
        trackhub.SubGroupDefinition(
            name='multi',
            label='multi',
            mapping=mapping,
        )
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
        visibility='full',
        tracktype='bigWig',
        short_label='Signal')
    composite.add_view(signal_view)

    priority = 0
    for experiment_name, experiment in experiments.iterrows():
        extra = {}
        if 'color' in experiment.keys():
            extra['color'] = experiment['color']
        for library_id in experiment['replicates']:
            row = libraries.loc[library_id]
            if 'color' in row:
                extra['color'] = row.color
            for track_type in track_types:
                track = trackhub.Track(
                    url=make_bigwig_url(baseurl, row, track_type),
                    name="{:03d}".format(priority) + '_' + row.analysis_name + '_' + track_type,
                    visibility='full',
                    tracktype='bigWig',
                    subgroups={
                        'experiment': sanitize_name(experiment_name),
                        'multi': track_type},
                    priority=priority,
                    **extra
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


def detect_genome_column_name(libraries):
    for name in ["genome_name", "genome"]:
        if name in libraries.columns:
            return name

    raise ValueError("Unrecognized library file, need a genome_name or genome column")


def build_ucsc_name_map_for_libraries(libraries):
    genome_column = detect_genome_column_name(libraries)
    unique_names = libraries[genome_column].unique()
    ucsc_name_map = {x: ucsc_genome_conversion(x) for x in unique_names}

    return ucsc_name_map


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
        'GRCh38-V29-male': 'hg38',
        'mm10-M21-male': 'mm10',
    }
    if genome_name in conversions:
        return conversions[genome_name]
    else:
        return genome_name


if __name__ == '__main__':
    main()
