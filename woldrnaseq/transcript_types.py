#!/usr/bin/python3
"""Read annotation bam and count reads per gene_type
"""
from argparse import ArgumentParser
from collections import Counter
import logging
from pathlib import Path
import pandas
import pysam
import time

from .common import (
    add_debug_arguments,
    add_metadata_arguments,
    add_separator_argument,
    add_version_argument,
    configure_logging,
    get_seperator,
)

from .models import (
    load_library_tables,
    load_experiments,
    find_library_bam_file,
)

logger = logging.getLogger(__name__)


def main(cmdline=None):
    parser = make_parser()
    args = parser.parse_args(cmdline)

    configure_logging(args)

    transcript_type_map = make_transcript_type_map(args.gtf_cache)

    if len(args.filenames) > 0:
        if args.output is None:
            parser.perror("Output filename is required when listing bam files directly")

        scores = make_transcript_type_scores(args.filenames, transcript_type_map)
        scores.to_csv(args.output, sep="\t")

    if not (args.libraries is None or args.experiments is None):
        sep = get_seperator(args.sep)
        libraries = load_library_tables(args.libraries, sep=sep)
        experiments = load_experiments(args.experiments, sep=sep)

        for i, experiment in experiments.iterrows():
            logging.info("Processing: %s", experiment.name)
            scores = make_experiment_transcript_type_scores(
                experiment, libraries, transcript_type_map
            )
            name = "{}_gene_type.tsv".format(experiment.name)
            scores.to_csv(name, sep="\t")
    elif args.libraries is None and args.experiments is None:
        # Neither provided
        pass
    else:
        # only one provided
        parser.perror(
            "You need to provide both a libraries and experiment table to use this mode")


def make_parser():
    parser = ArgumentParser()
    parser.add_argument(
        "filenames", nargs="*", help="Names of transcriptome bam files to score"
    )
    parser.add_argument(
        "-o", "--output", help="Name of output score table for directly read bam files"
    )
    add_metadata_arguments(parser)
    parser.add_argument(
        "--gtf-cache", required=True, help="name of gtf-cache file to read"
    )
    add_separator_argument(parser)
    add_version_argument(parser)
    add_debug_arguments(parser)
    return parser


def make_transcript_type_map(cache_filename):
    tstart = time.monotonic()
    logger.info("Loading GTF Cache")
    type_name = "gene_type"
    store = pandas.HDFStore(cache_filename)
    trna = store.select(
        "gtf", columns=["transcript_id", type_name], where=["type==tRNA"]
    )
    transcripts = store.select(
        "gtf", columns=["transcript_id", type_name], where=["type==transcript"]
    )
    spikes = store.select("gtf", columns=["transcript_id"], where=["source==spikein"])
    store.close()

    transcript_type_map = {k: "spikein" for k in spikes["transcript_id"]}
    transcript_series = transcripts.set_index("transcript_id")[type_name]
    transcript_type_map.update(transcript_series.to_dict())
    trna_series = trna.set_index("transcript_id")[type_name]
    transcript_type_map.update(trna_series.to_dict())
    assert len(transcript_type_map) == (transcripts.shape[0] + spikes.shape[0] + trna.shape[0])
    logging.debug("Loading finished {:.3} sec".format(time.monotonic() - tstart))
    return transcript_type_map


def make_experiment_transcript_type_scores(experiment, libraries, transcript_type_map):
    scores = {}
    for library_id in experiment.replicates:
        library = libraries.loc[library_id]
        anno = find_library_bam_file(library, "transcriptome")
        scores[library_id] = score_bam_transcript_type(anno, transcript_type_map)

    return pandas.DataFrame(scores)


def make_transcript_type_scores(urls, transcript_type_map):
    scores = {}
    for url in urls:
        basename = Path(url).name
        scores[basename] = score_bam_transcript_type(url, transcript_type_map)

    return pandas.DataFrame(scores)


def score_bam_transcript_type(alignment_filename, transcript_type_map):
    tstart = time.monotonic()
    logger.info("Counting {}".format(alignment_filename))
    counts = Counter()
    with pysam.AlignmentFile(alignment_filename, "r") as aligned:
        for read in aligned.fetch(until_eof=True):
            if not (read.is_secondary or read.is_unmapped or read.is_qcfail or read.is_duplicate):
                transcript_type = transcript_type_map.get(read.reference_name, "no id")
                counts[transcript_type] += 1

    logging.debug(
        "Finished counting {} {:.3f}".format(
            alignment_filename, time.monotonic() - tstart
        )
    )
    return counts


if __name__ == "__main__":
    main()
