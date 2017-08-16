#!/usr/bin/python3
from __future__ import print_function, absolute_import
import argparse
import os
import logging
from urllib.parse import urljoin

from woldrnaseq.common import (
    add_default_path_arguments,
    add_debug_arguments,
    add_version_argument,
    configure_logging,
    find_fastqs,
    get_seperator,
    validate_args
)
from woldrnaseq.models import load_library_tables
from woldrnaseq.version import get_git_version

logger = logging.getLogger('make_tracks')

def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    if args.version:
        parser.exit(0, 'version: %s\n' % (get_git_version(),))

    sep = get_seperator(args.sep)

    library_filenames = args.libraries
    if len(library_filenames) == 0:
        parser.error('Need library information table')

    libraries = load_library_tables(library_filenames, sep)
    
    custom_tracks = []
    for library_id, library in libraries.iterrows():
        if args.bigwig:
            custom_tracks.extend(
                make_bigwig_custom_tracks(library, args.web_root, args.root))

    print(os.linesep.join(custom_tracks))


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-w', '--web-root', required=True,
                        help='base url for tracks')
    parser.add_argument('--root', default=os.getcwd(),
                        help='root directory to search for tracks')
    parser.add_argument('--bigwig', action='store_true', default=False,
                        help='generate track blocks for bigwigs')
    
    parser.add_argument('-s', '--sep', choices=['TAB',','], default='TAB')
    parser.add_argument('-l', '--libraries', action='append', default=[])

    add_version_argument(parser)
    add_debug_arguments(parser)

    return parser


def make_bigwig_custom_tracks(library, web_root, analysis_root):
    """Generate a bigwig custom track record

    :param Series library: row from a library table DataFrame
    :param str web_root: base url to be prepended to paths
    :param str analysis_root: base directory to look for track files
    """
    track_template = 'track type=bigWig name={library_id} description={description} visibility=full color=255,0,0 bigDataUrl={url}'

    tracks = []
    for signal_type in ['uniq', 'all']:
        pathname = make_bigwig_track_name(library, signal_type, analysis_root)
        analysis_name = make_analysis_name(library)
        url = web_root +  pathname.replace(analysis_root, '')
        track = track_template.format(library_id=library.name,
                                      description=analysis_name,
                                      url=url)

        tracks.append(track)
    return tracks


def make_bigwig_track_name(library, signal_type, analysis_root):
    """Generate the base path where the bigwig track is
    """
    assert signal_type in ('uniq', 'all')
    
    analysis_name = make_analysis_name(library)
    genome_triplet = make_genome_triplet(library)
    track_name = analysis_name + '-' + genome_triplet + '_' + signal_type + '.bw'

    for pathname in [ os.path.join(library.analysis_dir, track_name), track_name]:
        abspathname = os.path.join(analysis_root, pathname)
        if os.path.exists(abspathname):
            return pathname

    logger.warning("Couldn't find track file %s", track_name)


def make_analysis_name(library):
    """generate analysis name from library table

    FIXME: this is duplicate functionaliy in make_star_rsem_dag
    """
    if 'analysis_name' in library.index:
        return library.analysis_name
    else:
        return library.analysis_dir


def make_genome_triplet(library):
    """Make genome triplet string from library table row

    FIXME: this is duplicate functionality in make_star_rsem_dag
    """
    return '-'.join([library.genome, library.annotation, library.sex])
    
if __name__ == '__main__':
    main()
    
