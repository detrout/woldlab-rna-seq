#!/usr/bin/python3
"""Download fastqs from illumina runfolders.
"""
from argparse import ArgumentParser
from datetime import datetime
from collections import namedtuple
import hashlib
import logging
import os
from pathlib import Path
import pandas
import sys
import requests
from urllib.parse import urljoin
from lxml.html import parse
import gzip
from xopen import xopen

from .models import (
    load_library_tables,
)
from .common import (
    add_metadata_arguments,
    add_debug_arguments,
    add_separator_argument,
    configure_logging,
    get_seperator,
)
from htsworkflow.util.api import (
    add_auth_options,
    make_auth_from_opts,
    HtswApi,
)

logger = logging.getLogger('downloader')


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)
    configure_logging(args)

    apidata = make_auth_from_opts(args)
    api = HtswApi(args.host, apidata)

    sep = get_seperator(args.sep)
    libraries = load_library_tables(args.libraries, sep=sep)

    #download_fastqs(args.url, args.flowcell, libraries)
    fragments = []
    flowcell_urls = {}
    library_details = {}
    for library_id, row in libraries.iterrows():
        library = api.get_library(library_id)
        for lane in library.get('lane_set', []):
            if lane['status'] in ('Good', 'Unknown'):
                flowcell_id = lane['flowcell']

                if args.flowcell is not None and flowcell_id not in args.flowcell:
                    # if there's a flowcell filter skip other flowcells
                    continue

                if flowcell_id not in flowcell_urls:
                    flowcell_entry = find_flowcell(args.root_url, flowcell_id)
                    projects = list(find_projects(flowcell_entry.url))
                    flowcell_urls[flowcell_id] = projects

                for project_entry in flowcell_urls[flowcell_id]:
                    _, project_id, project_index = project_entry.name.split('_')
                    if project_id == library_id:
                        fastq_entries = {}
                        for sample_entry in find_sample_url(project_entry.url):
                            for fastq in find_fastqs(sample_entry.url):
                                shortened_name = make_short_fastq_name(fastq.name)
                                fastq_entries.setdefault(shortened_name, []).append(fastq.url)

                        for shortened_name in fastq_entries:
                            target = Path(row.read_1).parent / shortened_name
                            fragments.extend(fast_download_merged_fastq(target, fastq_entries[shortened_name]))

    pandas.DataFrame(fragments, columns=['name', 'url', 'length', 'md5']).to_csv('fragments.csv')


def make_parser():
    parser = ArgumentParser()
    #add_metadata_arguments(parser)
    parser.add_argument('-l', '--libraries', default=[], action='append',
                        help='library metadata table to load')
    parser.add_argument('-f', '--flowcell', default=[], action='append',
                        help='limit to listed flowcells, otherwise try everything')
    parser.add_argument('-u', '--root-url', help='root url to search')
    add_separator_argument(parser)
    add_auth_options(parser)
    add_debug_arguments(parser)
    return parser


def download_fastqs(root_url, flowcell, libraries):
    for runfolder in find_flowcell(root_url, flowcell):
        logger.debug('Scanning runfolder {}'.format(runfolder.name))
        for project in find_projects(runfolder.url):
            name = project.name.replace('Project_', '')
            for library_id, library in libraries.iterrows():
                if name.startswith(library_id):
                    fastq_entries = {}
                    for sample in find_sample_url(project.url):
                        for fastq in find_fastqs(sample.url):
                            shortened_name = make_short_fastq_name(fastq.name)
                            fastq_entries.setdefault(shortened_name, []).append(fastq)

                    if len(fastq_entries) > 0:
                        for shortened_name in fastq_entries:
                            read_fastqs = sorted(fastq_entries[shortened_name], key=lambda f: f.name)
                            if not os.path.exists(library.analysis_dir):
                                os.mkdir(library.analysis_dir)
                            target = os.path.join(library.analysis_dir, shortened_name)
                            download_merged_fastq(target, [f.url for f in read_fastqs])


def fast_download_merged_fastq(fastq_name, fastq_urls):
    fragments = []
    BLOCK_SIZE = 65535
    if os.path.exists(fastq_name):
        logger.error('{} already exists'.format(fastq_name))
    else:
        with open(fastq_name, mode='wb') as outstream:
            block = bytearray(BLOCK_SIZE)
            full_digest = hashlib.md5()
            total_count = 0
            for url in fastq_urls:
                logger.debug('Downloading: {}'.format(url))
                fragment_hash = hashlib.md5()
                count = 0
                response = requests.get(url, stream=True)
                for block in response.iter_content(BLOCK_SIZE):
                    outstream.write(block)
                    count += len(block)
                    fragment_hash.update(block)
                    full_digest.update(block)
                logger.debug('Read {} bytes. md5={}'.format(count, fragment_hash.hexdigest()))
                fragments.append((fastq_name, url, count, fragment_hash.hexdigest()))
                total_count += count
        logger.debug('Final size {} bytes. md5={}'.format(total_count, full_digest.hexdigest()))
    print(fragments)
    return fragments

def download_merged_fastq(fastq_name, fastq_urls):
    BLOCK_SIZE = 65535
    if os.path.exists(fastq_name):
        logger.error('{} already exists'.format(fastq_name))
    else:
        with xopen(fastq_name, mode='wb', compresslevel=9, threads=6) as outstream:
            block = bytearray(BLOCK_SIZE)
            for url in fastq_urls:
                logger.debug('Downloading: {}'.format(url))
                fragment_hash = hashlib.md5()
                count = 0
                response = requests.get(url, stream=True)
                with gzip.open(response.raw, 'rb') as instream:
                    #read = instream.readinto(block)
                    #while read != 0:
                    #    outstream.write(block)
                    #    count += read
                    #    read = instream.readinto(block)
                    #    fragment_hash.update(block)
                    for line in instream:
                        count += len(line)
                        outstream.write(line)
                        fragment_hash.update(line)
                logger.debug('Read {} bytes. md5={}'.format(count, fragment_hash.hexdigest()))


def find_flowcell(root_url, flowcell):
    for entry in parse_apache_dirindex(root_url):
        if entry.link_type == 'subdirectory' and flowcell in entry.name:
            return entry


def find_projects(url, maxdepth=1):
    if maxdepth < 0:
        return

    for entry in parse_apache_dirindex(url):
        if entry.link_type == 'subdirectory' and entry.name.startswith('Project_'):
            yield entry
        elif entry.link_type == 'subdirectory':
            yield from find_projects(entry.url, maxdepth-1)


def find_sample_url(url):
    """Return the Sample subdirectory of a illumina project runfolder"""
    for sample in parse_apache_dirindex(url):
        if sample.link_type == 'subdirectory' and sample.name.startswith('Sample_'):
            yield sample


def find_fastqs(url):
    """Search apache dirindex for fastq files"""
    for entry in parse_apache_dirindex(url):
        if entry.link_type == 'compressed' and entry.name.endswith('.fastq.gz'):
            yield entry


direntry = namedtuple('direntry', ['link_type', 'name', 'url', 'last_modified'])


def parse_apache_dirindex(url):
    """Read an apache dirindex and return the directories
    """
    icon_src_type_map = {
        '/icons/folder.gif': 'subdirectory',
        '/icons/back.gif': 'parent',
        '/icons/unknown.gif': 'file',
        '/icons/text.gif': 'file',
        '/icons/compressed.gif': 'compressed',
    }

    root = parse(url)
    for row in root.xpath('/html/body/table/tr'):
        td = row.getchildren()
        icon = td[0].getchildren()
        if len(icon) != 1:
            continue
        icon = icon[0]
        if icon.tag != 'img':
            continue
        icon_src = icon.attrib.get('src')
        link_type = icon_src_type_map.get(icon_src)
        a = td[1].getchildren()
        if len(a) != 1:
            continue
        a = a[0]
        name = a.text
        href = a.attrib.get('href')
        if href is not None:
            href = urljoin(url, href)
        last_modified = td[2].text
        if last_modified is not None and len(last_modified) > 0:
            last_modified = last_modified.strip()
            if len(last_modified) > 0:
                last_modified = datetime.strptime(last_modified, '%d-%b-%Y %H:%M')
            else:
                last_modified = None
        yield direntry(link_type, name, href, last_modified)


def make_short_fastq_name(fastq_name):
    """Generate shortened name for combined fastq files
    """
    ext = '.fastq.gz'

    if not fastq_name.endswith(ext):
        raise ValueError('Unexpected fastq extension {}'.format(fastq_name))

    fastq_name = fastq_name[:-len(ext)]
    parts = fastq_name.split('_')
    if not parts[2].startswith('L'):
        raise ValueError('Unexpected fastq lane {}'.format(fastq_name))

    if not parts[3].startswith('R'):
        raise ValueError('Unexpected fastq read {}'.format(fastq_name))

    return '_'.join((parts[0], parts[1], parts[2], parts[3])) + ext


if __name__ == '__main__':
    sys.exit(main())
