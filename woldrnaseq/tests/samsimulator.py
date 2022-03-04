"""Utilities to help generate some test data
"""
from collections import abc
from contextlib import contextmanager
import itertools
import os
import tempfile
import pysam


# Op BAM Description Consumes Consumes
#                    query    reference
# M 0 alignment match (can be a sequence match  yes yes
#     or mismatch)
# I 1 insertion to the reference                yes no
# D 2 deletion from the reference               no yes
# N 3 skipped region from the reference         no yes
# S 4 soft clipping (clipped sequences present in SEQ) yes no
# H 5 hard clipping (clipped sequences NOT present in SEQ) no no
# P 6 padding (silent deletion from padded reference) no no
# = 7 sequence match yes yes
# X 8 sequence mismatch yes yes


def tokenize_cigar(cigar):
    tokens = []
    number = []
    for c in cigar:
        if c.isdigit():
            number.append(c)
        else:
            tokens.append((int("".join(number)), c))
            number = []

    return tokens


# References from https://samtools.github.io/hts-specs/SAMv1.pdf
# SAM fields
# Col Field Type   Regexp/Range                Brief description
# 1   QNAME String [!-?A-~]{1,254}             Query template NAME
# 2   FLAG  Int    [0, 216 − 1]                bitwise FLAG
# 3   RNAME String \*|[:rname:∧*=][:rname:]*   Reference sequence NAME11
# 4   POS   Int    [0, 231 − 1]                1-based leftmost mapping POSition
# 5   MAPQ  Int    [0, 28 − 1]                 MAPping Quality
# 6   CIGAR String \*|([0-9]+[MIDNSHPX=])+     CIGAR string
# 7   RNEXT String \*|=|[:rname:∧*=][:rname:]* Reference name of the mate/next read
# 8   PNEXT Int    [0, 231 − 1]                Position of the mate/next read
# 9   TLEN  Int    [−231 + 1, 231 − 1]         observed Template LENgth
# 10  SEQ   String \*|[A-Za-z=.]+              segment SEQuence
# 11  QUAL  String [!-~]+                      ASCII of Phred-scaled base QUALity+33
# 12  TAGS?

# Defined SAM flags
# 1    0x1   template having multiple segments in sequencing
# 2    0x2   each segment properly aligned according to the aligner
# 4    0x4   segment unmapped
# 8    0x8   next segment in the template unmapped
# 16   0x10  SEQ being reverse complemented
# 32   0x20  SEQ of the next segment in the template being reverse complemented
# 64   0x40  the first segment in the template
# 128  0x80  the last segment in the template
# 256  0x100 secondary alignment
# 512  0x200 not passing filters, such as platform/vendor quality controls
# 1024 0x400 PCR or optical duplicate
# 2048 0x800 supplementary alignment


def make_sam_flags(
    multi=False,
    aligned=True,
    unmapped=False,
    next_unmapped=False,
    is_reverse=False,
    next_is_reverse=False,
    first_segment=False,
    last_segment=False,
    secondary=False,
    qc_fail=False,
    duplicate=False,
    supplementary=False,
    **kwargs,
):
    return (
        0x1 * int(multi)
        | 0x2 * int(aligned)
        | 0x4 * int(unmapped)
        | 0x8 * int(next_unmapped)
        | 0x10 * int(is_reverse)
        | 0x20 * int(next_is_reverse)
        | 0x40 * int(first_segment)
        | 0x80 * int(last_segment)
        | 0x100 * int(secondary)
        | 0x200 * int(qc_fail)
        | 0x400 * int(duplicate)
        | 0x800 * int(supplementary)
    )


class SAMGenerator:
    def __init__(self):
        self._reads = []
        self._sq = {}

    def add_single_read(
        self,
        name,
        chromosome,
        pos,
        flag=None,
        mapq=None,
        cigar=None,
        rnext=None,
        pnext=None,
        tlen=None,
        seq=None,
        qual=None,
        tags=None,
        **kwargs,
    ):
        if seq is None:
            seq = "A" * 10
        if qual is None:
            qual = "@" * len(seq)
        assert len(seq) == len(qual)
        if mapq is None:
            mapq = 0
        if cigar is None:
            cigar = "{}M".format(len(seq))
        if rnext is None:
            rnext = "*"
        if pnext is None:
            pnext = "0"
        if tlen is None:
            tlen = len(seq)
        if tags is None:
            tags = {"NH": 1}

        if flag is None:
            if chromosome == "*":
                aligned = False
                tags["NH"] = 0
                pos = 0
                aligned = False
            else:
                aligned = True
                multi = True if tags["NH"] > 1 else False

            flag = make_sam_flags(multi=multi, aligned=aligned, **kwargs)

        items = [
            name,
            "{:g}".format(flag),
            chromosome,
            "{:g}".format(pos),
            "{:g}".format(mapq),
            cigar,
            rnext,
            pnext,
            "{:g}".format(tlen),
            seq,
            qual,
        ]
        for tag in tags:
            if tag == "NH":
                items.append("NH:i:{}".format(tags[tag]))

        cigar_len = sum((x[0] for x in tokenize_cigar(cigar)))
        current_sq = self._sq.get(chromosome, 0)
        self._sq[chromosome] = max(current_sq, pos + cigar_len, pos + len(seq))
        self._reads.append("\t".join(items))

    def add_paired_read(
        self,
        name,
        chromosome,
        pos,
        flag=None,
        mapq=None,
        cigar=None,
        rnext=None,
        pnext=None,
        tlen=None,
        seq=None,
        qual=None,
        tags=None,
        is_reverse=None,
        **kwargs,
    ):
        if not isinstance(pos, abc.Sequence):
            raise ValueError("Adding paired ends requires multiple locations")

        if seq is None:
            seq = []
            for base, count in zip(
                itertools.cycle("AGCT"),
                [
                    10,
                ]
                * len(pos),
            ):
                seq.append(base * count)
        if isinstance(seq, str):
            raise ValueError("Paired end reads need a list of sequences")
        if len(pos) != len(seq):
            raise ValueError("Length of sequence list and positions must be equal")

        if qual is None:
            qual = ["@" * len(x) for x in seq]

        if len(seq) != len(qual):
            raise ValueError("Lengths of sequence and qualities must be equal")
        for s, q in zip(seq, qual):
            if len(s) != len(q):
                raise ValueError(
                    "Lengths of sequence and qualities pairs must equal {} {}".format(
                        s, q
                    )
                )

        if is_reverse is None:
            is_reverse = [False] * len(pos)
        elif len(is_reverse) != len(pos):
            raise ValueError("Lengths of is_reverse and position must be the same")

        if mapq is None:
            mapq = [0] * len(pos)
        elif len(mapq) != len(pos):
            raise ValueError("Lengths of mapq must match number of positions")

        if cigar is None:
            cigar = ["{}M".format(len(x)) for x in seq]
        elif len(cigar) != len(pos):
            raise ValueError("Lengths of cigar and positions must be the same")

        if rnext is None:
            rnext = ["="] * len(seq)
        elif len(rnext) != len(pos):
            raise ValueError("Lengths of next reads and positions must be the same")

        if pnext is None:
            pnext = [str(p) for p in itertools.chain(pos[1:], [pos[0]])]
        elif len(pnext) != len(pos):
            raise ValueError(
                "Lengths of next next positions and positions must be the same"
            )

        if tlen is None:
            strand = [-1 if reverse else 1 for reverse in is_reverse]
            tlen = [r * len(s) for r, s in zip(strand, seq)]
        elif len(tlen) != len(pos):
            raise ValueError("Lenghts of tlen need to match number of positions")

        if tags is None:
            tags = {"NH": 1}

        if flag is None:
            if chromosome == "*":
                aligned = False
                tags["NH"] = 0
                pos = 0
                aligned = False
            else:
                aligned = True
                multi = True if tags["NH"] > 1 else False

            flag = []
            for i in range(len(pos)):
                is_secondary = True if i > 0 else False
                kwargs["is_reverse"] = is_reverse[i]
                flag.append(
                    make_sam_flags(
                        multi=multi,
                        aligned=aligned,
                        is_secondary=is_secondary,
                        **kwargs,
                    )
                )

        for i in range(len(pos)):
            kwargs["is_reverse"] = is_reverse[i]
            self.add_single_read(
                name,
                chromosome,
                pos[i],
                flag=flag[i],
                mapq=mapq[i],
                cigar=cigar[i],
                rnext=rnext[i],
                pnext=pnext[i],
                tlen=tlen[i],
                seq=seq[i],
                qual=qual[i],
                tags=tags,
                **kwargs,
            )

    def make_header(self):
        yield "\t".join(("@HD", "VN:1.4", "SO:coordinate"))
        yield "\t".join(("@PG", "ID:DS", "PN:datasimulator"))
        for key in sorted(self._sq):
            if key != "*":
                yield "\t".join(
                    ("@SQ", "SN:{}".format(key), "LN:{}".format(self._sq[key]))
                )

    @contextmanager
    def to_alignedfile(self):
        with tempfile.NamedTemporaryFile(mode="w+t", suffix=".sam") as sam:
            sam_name = sam.name
            for line in itertools.chain(self.make_header(), self._reads):
                sam.write(line)
                sam.write(os.linesep)
            sam.seek(0)

            with pysam.AlignmentFile(sam.name, "r") as instream:
                yield instream
