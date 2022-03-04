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
import sys
import re
import requests
from urllib.parse import urljoin, urlsplit
from lxml.html import parse
import gzip
from xopen import xopen

from .models import (
    load_library_tables,
)
from .common import (
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

    if len(args.libraries) == 0:
        parser.error('No library metadata files specified, so nothing to do')

    sep = get_seperator(args.sep)
    libraries = load_library_tables(args.libraries, sep=sep)

    if len(libraries) == 0:
        parser.error('No libraries loaded from metadata files, nothing to do.')

    apidata = make_auth_from_opts(args)
    api = HtswApi(args.host, apidata)
    target_flowcells = get_library_flowcells(api, libraries, args.flowcell)

    runfolders = search_for_runfolders(args.root_url, target_flowcells)

    target_fastqs = {}
    for flowcell_id in runfolders:
        url = runfolders[flowcell_id]
        runfolder_fastqs = search_runfolder_for_fastqs(
            url, libraries, merge_lanes=args.merge_lanes)
        _update_fastq_entries(target_fastqs, runfolder_fastqs)

    download_fastqs_by_entry(target_fastqs, libraries)

def make_parser():
    parser = ArgumentParser()
    parser.add_argument('-l', '--libraries', default=[], action='append',
                        help='library metadata table to load')
    parser.add_argument('-f', '--flowcell', default=[], action='append',
                        help='limit to listed flowcells, otherwise try everything')
    parser.add_argument('-u', '--root-url', default=[], action="append",
                        help='root url to search')
    parser.add_argument('--merge-lanes', default=False, action='store_true',
                        help='merge multiple lanes into one file')
    add_separator_argument(parser)
    add_auth_options(parser)
    add_debug_arguments(parser)
    return parser


def get_library_flowcells(htsw, libraries, limit=[]):
    """Retrieve flowcells used for a list of libraries was run on

    Parameters:
      htsw: initialized HtswApi object
      libraries: pandas.DataFrame with library ids as index.
      limit: list of flowcell ids to include or everything if None.
    Returns:
      {flowcell1: {library_id1, library_id2}, flowcell2: ...}
      if the library id wasn't found or was filtered out it may not appear in
      in the return dictionary
    """
    flowcells = set()
    for library_id, row in libraries.iterrows():
        library = htsw.get_library(library_id)
        for lane in library.get('lane_set', []):
            if lane['status'] in ('Good', 'Unknown'):
                flowcell_id = lane['flowcell']

                if len(limit) > 0 and flowcell_id not in limit:
                    # if there's a flowcell filter skip other flowcells
                    continue
                else:
                    flowcells.add(flowcell_id)
            else:
                logger.info("Skipping {} {} because status {}".format(
                    library_id, lane["flowcell"], lane["status"]))

    logger.info("Found {}".format(",".join(flowcells)))
    return flowcells


def search_for_runfolders(runfolder_roots, flowcells):
    """Search a list of apache indexes for a list of flowcell folders

    Parameters:
      - urls: (list) apache directory indexes urls containing runfolders
      - flowcells: (list) flowcell ids to look for.

    Returns:
      dictionary of flowcell ids and the urls to find the runfolder
      {flowcell: url, flowcell: url}
    """
    flowcell_runfolder = {}
    for url in runfolder_roots:
        if not url.endswith("/"):
            url += "/"
        for entry in parse_apache_dirindex(url):
            if entry.link_type == "subdirectory":
                for flowcell in flowcells:
                    if entry.name.endswith("{}/".format(flowcell)):
                        flowcell_runfolder[flowcell] = entry.url

    return flowcell_runfolder


def search_runfolder_for_fastqs(runfolder_url, libraries, merge_lanes=False):
    logger.debug('Scanning runfolder {}'.format(runfolder_url))
    fastq_entries = {}

    contents = list(parse_apache_dirindex(runfolder_url))
    runfolder_type = guess_runfolder_type(contents)
    logger.debug("Guessing runfolder type: {}".format(runfolder_type))

    if runfolder_type == "HiSeq":
        for entry in contents:
            if entry.link_type == "subdirectory" and entry.name.startswith("Unaligned"):
                _update_fastq_entries(fastq_entries, find_hiseq_fastqs(entry.url, libraries, merge_lanes))
    elif runfolder_type == "nextseq":
        for entry in contents:
            if entry.link_type == "subdirectory" and entry.name.startswith("Analysis"):
                _update_fastq_entries(fastq_entries, find_nextseq_fastqs(entry.url, libraries, merge_lanes))
    else:
        raise ValueError("Unrecognized folder type {}".format(runfodler_url))

    return fastq_entries


def _update_fastq_entries(merged_fastqs, new_fastqs):
    """Add new fastq urls to the dictionary indexed by destination file.

    This is intended to group fastqs by what our desired output file is,
    grouping by things like lane, chunk, or even flowcell if desired.
    """
    for name in new_fastqs:
        if name in merged_fastqs:
            logger.warning(
                "Multiple locations are contributing to {}".format(name))
        merged_fastqs.setdefault(name, []).extend(new_fastqs[name])


def find_nextseq_fastqs(url, libraries, merge_lanes):
    """Search nextseq runfolder for fastqs
    """
    logger.debug("Scanning for fastqs in {}".format(url))

    #Nextseq runfolders appear to be much flatter than the hi-seq runfolders
    fastq_entries = {}
    #21154_index34_S1_R1_001.fastq.gz
    fastq_re = re.compile("(?P<library_id>[^_]+)_.*\.fastq\.gz")
    for fastq_dir in find_nextseq_fastq_dir(url):
        for entry in parse_apache_dirindex(fastq_dir.url):
            match = fastq_re.match(entry.name)
            if match:
                logger.debug("Found fastq {}".format(entry.url))
                library_id = match.group("library_id")
                if library_id in libraries.index:
                    fastq_entry = FastqFragment(entry.url, merge_lanes=merge_lanes)
                    fastq_entries.setdefault(fastq_entry.key, []).append(fastq_entry)
    return fastq_entries

def find_nextseq_fastq_dir(url):
    logger.debug("Scanning {}".format(url))
    contents = list(parse_apache_dirindex(url))
    done = False
    for entry in contents:
        if entry.link_type == "subdirectory" and entry.name == "fastq/":
            logger.debug("Found fastq directory at {}".format(url))
            done = True
            yield entry

    if not done:
        for entry in contents:
            if entry.link_type == "subdirectory":
                yield from find_nextseq_fastq_dir(entry.url)



def find_hiseq_fastqs(url, libraries, merge_lanes):
    """Recurse though a runfolder looking for fastqs associated with libraries
    """
    fastq_entries = {}
    project_re = re.compile("Project_(?P<library_id>[^_]+)")
    for project in find_projects(url):
        match = project_re.match(project.name)
        if match:
            library_id = match.group("library_id")
            if library_id in libraries.index:
                for sample in find_sample_url(project.url):
                    for fastq in find_fastqs(sample.url):
                        fastq_entry = FastqFragment(fastq.url, merge_lanes=merge_lanes)
                        fastq_entries.setdefault(fastq_entry.key, []).append(fastq_entry)
    return fastq_entries


def find_projects(url, maxdepth=1):
    """Look for Project directories
    """
    if maxdepth < 0:
        return

    print("find project scanning: {}".format(url))
    for entry in parse_apache_dirindex(url):
        # for hiseq 2500
        if entry.link_type == 'subdirectory' and entry.name.startswith('Project_'):
            yield entry


def find_sample_url(url):
    """Return the Sample subdirectory of a illumina project runfolder

    For the HiSeq runfolders its found within a Project_ directory
    """
    for sample in parse_apache_dirindex(url):
        if sample.link_type == 'subdirectory' and sample.name.startswith("Data"):
            yield sample
        # for hiseq 2500
        elif sample.link_type == 'subdirectory' and sample.name.startswith('Sample_'):
            yield sample


def find_fastqs(url):
    """Search apache dirindex for fastq files"""
    for entry in parse_apache_dirindex(url):
        if entry.link_type == 'compressed' and entry.name.endswith('.fastq.gz'):
            yield entry


def download_fastqs_by_entry(fastq_entries, libraries):
    if len(fastq_entries) > 0:
        for key in sorted(fastq_entries):
            read_fastqs = sorted(fastq_entries[key], key=lambda f: f.name)
            library = libraries.loc[read_fastqs[0].library_id]

            if not os.path.exists(library.analysis_dir):
                os.mkdir(library.analysis_dir)
            target = os.path.join(library.analysis_dir, read_fastqs[0].short_name)
            download_merged_fastq(target, [f.url for f in read_fastqs])


def download_merged_fastq(fastq_name, fastq_urls):
    if os.path.exists(fastq_name):
        logger.error('{} already exists'.format(fastq_name))
    else:
        with xopen(fastq_name, mode='wb', compresslevel=9, threads=6) as outstream:
            logger.info("Downloading: {}".format(fastq_name))
            for url in fastq_urls:
                logger.debug('Downloading part: {}'.format(url))
                fragment_hash = hashlib.md5()
                count = 0
                response = requests.get(url, stream=True)
                with gzip.open(response.raw, 'rb') as instream:
                    for line in instream:
                        count += len(line)
                        outstream.write(line)
                        fragment_hash.update(line)
                logger.debug('Read {} bytes. md5={}'.format(count, fragment_hash.hexdigest()))


direntry = namedtuple('direntry', ['link_type', 'name', 'url', 'last_modified'])


def parse_apache_dirindex(url):
    """Read an apache dirindex and return the directories
    """
    icon_src_type_map = {
        "/icons/blank.gif": "header",
        '/icons/folder.gif': 'subdirectory',
        '/icons/back.gif': 'parent',
        '/icons/unknown.gif': 'file',
        '/icons/text.gif': 'file',
        '/icons/compressed.gif': 'compressed',
    }

    try:
        root = parse(url)
    except TypeError as e:
        print("Unable to parse {}".format(url))
        print(e)
        return

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


def guess_runfolder_type(direntries):
    for entry in direntries:
        # the nextseq version of this file is named RunParameters
        if entry.name == "runParameters.xml":
            return "HiSeq"
        elif entry.name == "RTA3.cfg":
            return "nextseq"
    raise RuntimeError("Unknown folder type")


class FastqFragment:
    hiseq_re = re.compile("(?P<library_id>[^_]+)_(?P<index>[AGCT]+)_(?P<lane>L[0-9]+)_(?P<read>R[0-9])_(?P<chunk>[0-9]+)\.fastq\.gz")
    nextseq_re = re.compile("(?P<library_id>[^_]+)_((?P<index>.*)_)?(?P<sample>S[0-9]+)_(?P<read>R[0-9])_(?P<chunk>[0-9]+)\.fastq\.gz")
    patterns = [hiseq_re, nextseq_re]

    def __init__(self, url, merge_lanes=False):
        self.url = url
        self.merge_lanes = merge_lanes
        self.extension = ".fastq.gz"
        parts = urlsplit(url)
        self.name = Path(parts.path).name

        if not self.name.endswith(self.extension):
            raise ValueError('Unexpected fastq extension {}'.format(self.name))

        if self.name.startswith("Undetermined_"):
            raise ValueError(
                "We do not want to download fastqs that were not demultiplexed")

        for fastq_re in FastqFragment.patterns:
            match = fastq_re.match(self.name)
            if match:
                break
        if match is None:
            raise ValueError("filename does not match known fastq patterns")

        groups = match.groupdict()
        self.library_id = groups["library_id"]
        self.index = groups.get("index", None)
        if self.index is not None:
            self.index = self.index.replace("_", ".")
        self.lane = groups.get("lane", None)
        self.read = groups["read"]
        self.chunk = groups["chunk"]

    @property
    def key(self):
        if self.lane is None or self.merge_lanes is True:
            return (self.library_id, self.index, self.read)
        else:
            return (self.library_id, self.index, self.lane, self.read)

    def __repr__(self):
        return self.url

    @property
    def short_name(self):
        parts = []
        parts.append(self.library_id)
        if self.index is not None:
            parts.append(self.index)
        if not (self.lane is None or self.merge_lanes):
            parts.append(self.lane)
        parts.append(self.read)

        return '_'.join(parts) + ".fastq.gz"


if __name__ == '__main__':
    sys.exit(main())
