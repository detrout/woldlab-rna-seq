#!/usr/bin/python3
"""Download fastqs from illumina runfolders.
"""
from argparse import ArgumentParser
from datetime import datetime
from collections import namedtuple
import hashlib
import logging
import os
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
    configure_logging,
    get_seperator,
)


logger = logging.getLogger(__name__)


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    libraries = load_library_tables(args.libraries)

    download_fastqs(args.url, args.flowcell, libraries)


def make_parser():
    parser = ArgumentParser()
    add_metadata_arguments(parser)
    parser.add_argument('-f', '--flowcell',
                        help='Name of flowcell to look for')
    parser.add_argument('-u', '--url', help='root url to search')
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
            yield entry


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

    return '_'.join((parts[0], parts[1], parts[3])) + ext


if __name__ == '__main__':
    sys.exit(main())
