"""
Utils for chain-related processing

Nae-Chyun Chen
Johns Hopkins University

Nancy Fisher Hansen
NIH/NHGRI

2021-2022
"""
import sys

import intervaltree
import pysam


def read_fasta(ref_fn: str) -> dict:
    """Reads a FASTA file as a dict.

    If given an empty path, return an empty dict.
    """
    ref = {}
    if ref_fn != "":
        f = pysam.FastaFile(ref_fn)
        for r in f.references:
            ref[r] = f[r].upper()
    return ref


def fasta_reader(ref_fn: str) -> pysam.FastaFile:
    if ref_fn:
        return pysam.FastaFile(ref_fn)
    else:
        return None


def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    nuc_trans = str.maketrans("AaCcGgTtN", "TTGGCCAAN")
    return seq.translate(nuc_trans)[::-1]


def compute_hamming_dist(
    forward: bool,
    ref1: dict,
    contig1: str,
    start1: int,
    end1: int,
    ref2: dict,
    contig2: str,
    start2: int,
    end2: int,
) -> float:
    if contig1 in ref1:
        s1 = ref1[contig1][start1:end1]
    else:
        print(f"Warning : {contig1} not in ref1", file=sys.stderr)
        return 0
    if contig2 in ref2:
        s2 = ref2[contig2][start2:end2]
        if not forward:
            s2 = reverse_complement(s2)
    else:
        print(f"Warning : {contig2} not in ref2", file=sys.stderr)
        return 0

    try:
        assert len(s1) == len(s2)
    except:
        print("Error: lengths do not match", file=sys.stderr)
        print(len(s1), len(s2), file=sys.stderr)
        print(f"{contig1}:{start1}-{end1}", file=sys.stderr)
        print(f"{contig2}:{start2}-{end2}", file=sys.stderr)
        exit(1)

    if len(s1) == 0:
        return 0
    idy = 0
    for i_s, s in enumerate(s1):
        if s == s2[i_s]:
            idy += 1
    idy /= len(s1)
    return idy


def get_target_entries(fn: str) -> dict:
    """Reads a chain file and gets target entries.

    Args:
        - fn: chain filename

    Returns:
        - a dict of:
            key: target contig
            value: target contig length
    """
    CC = ChainConst()
    dict_contig_length = {}
    with open(fn, "r") as f:
        for line in f:
            line = line.split()
            if len(line) == 13:
                dict_contig_length[line[CC.CHDR_TARGET]] = line[CC.CHDR_TLEN]
    return dict_contig_length


class ChainConst:
    CHDR_SCORE = 1
    CHDR_TARGET = 2
    CHDR_TLEN = 3
    CHDR_TSTRAND = 4
    CHDR_TSTART = 5
    CHDR_TEND = 6
    CHDR_QUERY = 7
    CHDR_QLEN = 8
    CHDR_QSTRAND = 9
    CHDR_QSTART = 10
    CHDR_QEND = 11
    CHDR_ID = 12
    C_SIZE = 0
    C_DT = 1
    C_DQ = 2


class Chain(ChainConst):
    def __init__(self, fields=None):
        self.ttree = intervaltree.IntervalTree()
        self.qtree = intervaltree.IntervalTree()
        if not fields:
            return
        self.score = int(fields[self.CHDR_SCORE])
        self.target = fields[self.CHDR_TARGET]
        self.tlen = int(fields[self.CHDR_TLEN])
        self.tstart = int(fields[self.CHDR_TSTART])
        self.toffset = self.tstart
        self.tend = int(fields[self.CHDR_TEND])
        self.tstrand = fields[self.CHDR_TSTRAND]
        self.query = fields[self.CHDR_QUERY]
        self.qlen = int(fields[self.CHDR_QLEN])
        self.qstrand = fields[self.CHDR_QSTRAND]
        self.qstart = int(fields[self.CHDR_QSTART])
        if self.qstrand == "+":
            self.qoffset = self.qstart
        else:
            self.qoffset = self.qlen - self.qstart
        self.qend = int(fields[self.CHDR_QEND])

        self.strand = (
            "+" if (self.tstrand == "+" and self.qstrand == "+") else "-"
        )
        self.id = fields[self.CHDR_ID]
        self.dt = 0
        self.dq = 0
        self.seglen = 0
        # Variables for converting uses
        self.cigarstring = ""
        self.nm_tag_val = 0
        self.total_num_match = 0
        self.total_block_t = 0

    def add_record(self, fields) -> None:
        segment_size = int(fields[self.C_SIZE])
        self.seglen += segment_size
        if segment_size > 0:
            self.ttree[self.toffset : self.toffset + segment_size] = (
                self.qoffset - self.toffset,
                self.dt,
                self.dq,
            )
            if self.strand == "+":
                self.qtree[self.qoffset : self.qoffset + segment_size] = 1
            else:
                self.qtree[self.qoffset - segment_size : self.qoffset] = 1

        if len(fields) == 3:
            self.dt = int(fields[self.C_DT])
            self.dq = int(fields[self.C_DQ])
        else:
            self.dt = 0
            self.dq = 0
        self.toffset += segment_size + self.dt
        if self.strand == "+":
            self.qoffset += segment_size + self.dq
        else:
            self.qoffset -= segment_size + self.dq

    def record_to_bed(self, fields: list, coord: str) -> str:
        segment_size = int(fields[0])
        if segment_size > 0:
            if coord == "target":
                return f"{self.target}\t{self.toffset}\t{self.toffset + segment_size}"
            elif coord == "query":
                if self.strand == "+":
                    return f"{self.query}\t{self.qoffset}\t{self.qoffset + segment_size}"
                else:
                    return f"{self.query}\t{self.qoffset - segment_size}\t{self.qoffset}"
        return ""

    def print_hdr(self) -> str:
        return (
            f"chain {self.score} "
            f"{self.target} {self.tlen} {self.tstrand} {self.tstart} {self.tend} "
            f"{self.query} {self.qlen} {self.qstrand} {self.qstart} {self.qend} "
            f"{self.id}\n"
        )

    def print_chain(self) -> str:
        msg = (
            f"chain {self.score} "
            f"{self.target} {self.tlen} {self.tstrand} {self.tstart} {self.tend} "
            f"{self.query} {self.qlen} {self.qstrand} {self.qstart} {self.qend} "
            f"{self.id}\n"
        )
        intervals = sorted(self.ttree.all_intervals)
        for i, intvl in enumerate(intervals):
            if i == 0:
                # If the size of the first interval is not zero
                if intvl.data[1:3] != (0, 0):
                    msg += (
                        f"{intvl.begin - self.tstart - intvl.data[1]}\t"
                        f"{intvl.data[1]}\t{intvl.data[2]}\n"
                    )
            else:
                pintvl = intervals[i - 1]
                msg += f"{pintvl.end - pintvl.begin}\t{intvl.data[1]}\t{intvl.data[2]}\n"

        # If the size of the last segment is greater than 0
        if intvl.end == self.tend:
            msg += f"{intvl.end - intvl.begin}\n"
        else:
            msg += (
                f"{intvl.end - intvl.begin}\t"
                f"{self.tend - intvl.end}\t"
                f"{self.qend - (intvl.end+intvl.data[0])}\n0"
            )
        return msg

    # Reset a Chain object at a given breakpoint
    def reset_at_break(self, dt: int, dq: int, tend: int, qend: int) -> None:
        self.ttree = intervaltree.IntervalTree()
        self.qtree = intervaltree.IntervalTree()
        self.tstart = self.tend + dt
        self.toffset = self.tstart
        self.tend = tend
        if self.strand == "+":
            self.qstart = self.qend + dq
            self.qoffset = self.qstart
            self.qend = qend
        else:
            self.qoffset -= dq
            # self.qoffset = self.qstart + dq
            # self.qoffset = self.qend - dq
            self.qstart = self.qlen - self.qoffset
            self.qend = qend
            # self.qend = self.qlen - qend
        # Update ID
        if len(self.id.split(".")) > 1:
            l = self.id.split(".")
            self.id = l[0] + "." + str(int(l[1]) + 1)
        else:
            self.id = self.id + ".1"

    # TODO: naechyun
    def try_merge(self, c):
        # print(len(list(self.stree)))
        # print(self.stree.items())
        # print(c.stree.items())

        self_sorted_intervals = sorted(self.ttree.all_intervals)
        c_sorted_intervals = sorted(c.stree.all_intervals)

        clone_tr = self.ttree.copy()
        for c_intvl in c_sorted_intervals:
            for i_s, s_intvl in enumerate(self_sorted_intervals):
                if c_intvl.begin > s_intvl.end or c_intvl.end < s_intvl.begin:
                    continue
                else:
                    # If offsets are the same, good to merge
                    if c_intvl.data[0] == s_intvl.data[0]:
                        if (
                            c_intvl.begin < s_intvl.begin
                            and c_intvl.end < s_intvl.end
                        ):
                            # Extend from beginning
                            # print(c_intvl)
                            # print(s_intvl)
                            start = min(s_intvl.begin, c_intvl.begin)
                            diff = s_intvl.begin - c_intvl.begin
                            data = (
                                s_intvl.data[0],
                                s_intvl.data[1] - diff,
                                s_intvl.data[2] - diff,
                            )
                            updated_s_intvl = intervaltree.Interval(
                                start, s_intvl.end, data
                            )
                            clone_tr.remove(s_intvl)
                            c.stree.remove(c_intvl)
                            clone_tr.add(updated_s_intvl)
                        elif (
                            c_intvl.begin > s_intvl.begin
                            and c_intvl.end > s_intvl.end
                        ):
                            # Extend from end
                            # This operation is trickier b/c the chain diffs are stored in the
                            # next interval. We need to pop the next interval and update it.
                            # print(c_intvl)
                            # print(s_intvl)
                            clone_tr.remove(s_intvl)
                            c.stree.remove(c_intvl)
                            end = max(s_intvl.end, c_intvl.end)
                            diff = c_intvl.end - s_intvl.end
                            clone_tr.add(
                                intervaltree.Interval(
                                    s_intvl.begin, end, s_intvl[2]
                                )
                            )
                            next_s_intvl = self_sorted_intervals[i_s + 1]
                            next_data = (
                                next_s_intvl.data[0],
                                next_s_intvl.data[1] - diff,
                                next_s_intvl.data[2] - diff,
                            )
                            updated_next_s_intvl = intervaltree.Interval(
                                next_s_intvl.begin, next_s_intvl.end, next_data
                            )
                            clone_tr.remove(next_s_intvl)
                            clone_tr.add(updated_next_s_intvl)
                        elif (
                            c_intvl.begin < s_intvl.begina
                            and c_intvl.end > s_intvl.end
                        ):
                            print(
                                "During merging, the new interval consumes the existing one",
                                file=sys.stderr,
                            )
                            print(
                                "This operation is not yet supported",
                                file=sys.stderr,
                            )
                            return False
                    else:
                        return False
        self.ttree = clone_tr
        # print('\nmerged:')
        # print(self.stree.items())
        sorted_intervals = sorted(self.ttree.all_intervals)

        # Update chain info
        self.score += c.score
        if sorted_intervals[0][0] < self.tstart:
            self.tstart = sorted_intervals[0][0]
            self.qstart = self.tstart + sorted_intervals[0][2]
        if sorted_intervals[-1][1] > self.tend:
            self.tend = sorted_intervals[-1][1]
            self.qend = self.tend + sorted_intervals[-1][2]

        return True

    def update_cigar_indel(
        self, deltaq: int, deltat: int, queryendpos: int = 0
    ) -> int:
        """Updates CIGAR string for a gap

        Args:
            - deltaq: the dq field in a chain line
            - deltat: the dt field in a chain line
            - queryendpos: pos_end of the processed segments

        Returns:
            updated queryendpos
        """
        if deltaq < 0:
            # TODO use `f"{deltaq=}"` when we don't support py3.7
            raise ValueError(
                f"deltaq should not be negative. Got deltaq: {deltaq}"
            )
        if deltat < 0:
            # TODO use `f"{deltat=}"` when we don't support py3.7
            raise ValueError(
                f"deltat should not be negative. Got deltat: {deltat}"
            )
        if deltaq > 0:
            self.cigarstring += f"{deltaq}I"
            self.nm_tag_val += deltaq
        if deltat > 0:
            self.cigarstring += f"{deltat}D"
            self.nm_tag_val += deltat
        if self.strand == "+":
            queryendpos += deltaq
        else:
            queryendpos -= deltaq
        return queryendpos

    def update_cigar_match(
        self,
        intvl: intervaltree.Interval,
        queryref: pysam.FastaFile,
        qname: str,
        targetref: pysam.FastaFile,
        tname: str,
        queryendpos: int = 0,
    ) -> int:
        """Updates CIGAR string for a matched chain segment.

        An interval represents zero-based, half-open start/end of a matched
        segment.

        Args:
            - intvl: the chain segment in an interval format
            - queryref: query FASTA object
            - qname: query contig name
            - targetref: target FASTA object
            - tname: target contig name
            - queryendpos

        Returns:
            updated queryendpos
        """
        segmentlength = intvl.end - intvl.begin
        # If reference sequences are provided
        if queryref and targetref:
            if self.strand == "+":
                alignseq_qstart = intvl.begin + intvl.data[0]
                alignseq_qend = intvl.end + intvl.data[0]
                alignseq_q = queryref.fetch(
                    reference=qname, start=alignseq_qstart, end=alignseq_qend
                ).upper()
                queryendpos += segmentlength
            else:
                alignseq_qend = intvl.begin + intvl.data[0]
                alignseq_qstart = alignseq_qend - (intvl.end - intvl.begin)
                alignseq_q = queryref.fetch(
                    reference=qname, start=alignseq_qstart, end=alignseq_qend
                ).upper()
                alignseq_q = reverse_complement(alignseq_q)
                queryendpos -= segmentlength

            alignseq_tstart = intvl.begin
            alignseq_tend = intvl.end
            alignseq_t = targetref.fetch(
                reference=tname, start=alignseq_tstart, end=alignseq_tend
            ).upper()

            nummatches = 0
            nummismatches = 0
            for i in range(segmentlength):
                if alignseq_q[i] == alignseq_t[i]:
                    if nummismatches > 0:
                        self.cigarstring += f"{nummismatches}X"
                        self.nm_tag_val += nummismatches
                        nummismatches = 0
                    nummatches += 1
                else:
                    if nummatches > 0:
                        self.cigarstring += f"{nummatches}="
                        nummatches = 0
                    nummismatches += 1
            if nummatches > 0:
                self.cigarstring += f"{nummatches}="
                self.total_num_match += nummatches
                self.total_block_t += nummatches
            elif nummismatches > 0:
                self.cigarstring += f"{nummismatches}X"
                self.nm_tag_val += nummismatches
                self.total_block_t += nummismatches
        else:
            # If no reference sequences are provided
            segmentlength = intvl.end - intvl.begin
            self.cigarstring += f"{segmentlength}M"
            self.total_num_match += segmentlength
            self.total_block_t += segmentlength
            if self.strand == "+":
                queryendpos += segmentlength
            else:
                queryendpos -= segmentlength
        return queryendpos

    def to_paf(
        self,
        targetref: pysam.FastaFile,
        queryref: pysam.FastaFile,
        preserve_breakpoint: bool = False,
    ) -> str:
        for i, intvl in enumerate(sorted(self.ttree.all_intervals)):
            # deal with indels first, because they are to the left of the
            # matching segment:
            if i > 0:
                deltat = intvl.data[1]
                self.total_block_t += deltat
                deltaq = intvl.data[2]

                # Reduces breakpoints
                if deltaq > 0 and deltat > 0 and not preserve_breakpoint:
                    min_delta = min([deltat, deltaq])
                    self.cigarstring += f"{min_delta}X"
                    if deltaq > min_delta:
                        self.cigarstring += f"{deltaq - min_delta}I"
                    elif deltat > min_delta:
                        self.cigarstring += f"{deltat - min_delta}D"
                    self.nm_tag_val += deltaq
                else:
                    self.update_cigar_indel(deltaq=deltaq, deltat=deltat)

            # update matched segment
            self.update_cigar_match(
                intvl=intvl,
                queryref=queryref,
                qname=self.query,
                targetref=targetref,
                tname=self.target,
            )

        if self.qstrand == "+":
            msg = (
                f"{self.query}\t{self.qlen}\t{self.qstart}\t{self.qend}\t{self.qstrand}\t"
                f"{self.target}\t{self.tlen}\t{self.tstart}\t{self.tend}\t{self.total_num_match}\t"
                f"{self.total_block_t}\t255\tcg:Z:{self.cigarstring}\tNM:i:{self.nm_tag_val}"
            )
        else:
            msg = (
                f"{self.query}\t{self.qlen}\t{self.qlen - self.qend}\t{self.qlen - self.qstart}\t{self.qstrand}\t"
                f"{self.target}\t{self.tlen}\t{self.tstart}\t{self.tend}\t{self.total_num_match}\t"
                f"{self.total_block_t}\t255\tcg:Z:{self.cigarstring}\tNM:i:{self.nm_tag_val}"
            )
        return msg

    def to_vcf(self, targetref, queryref) -> None:
        msg = ""
        rname = self.target
        qname = self.query

        intervals = sorted(self.ttree.all_intervals)
        for i, intvl in enumerate(intervals):
            # deal with indels first, because they are to the left of the matching segment:
            if i > 0:
                # calculate zero-based, half-open ref and alt start/end (subtracting one to add base prior to event)
                if self.strand == "+":
                    alignseq_qend = intvl.begin + intvl.data[0]
                    alignseq_qstart = (
                        intvl.begin + intvl.data[0] - intvl.data[2] - 1
                    )
                    alignseq_q = queryref.fetch(
                        reference=qname,
                        start=alignseq_qstart,
                        end=alignseq_qend,
                    ).upper()
                    qpos = alignseq_qstart + 1
                else:
                    alignseq_qstart = intvl.begin + intvl.data[0]
                    alignseq_qend = (
                        intvl.begin + intvl.data[0] + intvl.data[2] + 1
                    )
                    alignseq_q = queryref.fetch(
                        reference=qname,
                        start=alignseq_qstart,
                        end=alignseq_qend,
                    ).upper()
                    alignseq_q = reverse_complement(alignseq_q)
                    qpos = alignseq_qend

                # End position of indel allele is position to left of
                # matching segment's start, start has a base appended to the left:
                alignseq_tend = intvl.begin
                alignseq_tstart = intvl.begin - intvl.data[1] - 1
                alignseq_t = targetref.fetch(
                    reference=rname, start=alignseq_tstart, end=alignseq_tend
                ).upper()
                rpos = alignseq_tstart + 1

                msg += (
                    f"{rname}\t{rpos}\t.\t{alignseq_t}\t{alignseq_q}\t."
                    f"\tAUTO\tALN_SCORE={self.score};"
                    f"ALN_QUERY={qname};ALN_QPOS={qpos};"
                    f"ALN_STRAND={self.strand};ALN_DQ={intvl.data[2]};"
                    f"ALN_DT={intvl.data[1]}\n"
                )

            # now print SNPs within matched segment:

            # zero-based, half-open start/end of matched segment
            if self.strand == "+":
                alignseq_qstart = intvl.begin + intvl.data[0]
                alignseq_qend = intvl.end + intvl.data[0]
                alignseq_q = queryref.fetch(
                    reference=qname, start=alignseq_qstart, end=alignseq_qend
                ).upper()
            else:
                alignseq_qend = intvl.begin + intvl.data[0]
                alignseq_qstart = alignseq_qend - (intvl.end - intvl.begin)
                alignseq_q = queryref.fetch(
                    reference=qname, start=alignseq_qstart, end=alignseq_qend
                ).upper()
                alignseq_q = reverse_complement(alignseq_q)

            alignseq_tstart = intvl.begin
            alignseq_tend = intvl.end
            alignseq_t = targetref.fetch(
                reference=rname, start=alignseq_tstart, end=alignseq_tend
            ).upper()

            segmentlength = intvl.end - intvl.begin

            # print(f'Tstart {alignseq_tstart} Tend {alignseq_tend} Sstart {alignseq_sstart} Send {alignseq_send} Intvl {intvl.begin}:{intvl.end} Toffset {intvl.data[1]} Seglength {segmentlength}')

            for i in range(segmentlength):
                if alignseq_q[i] != alignseq_t[i]:
                    # convert to one-based start position:
                    rpos = alignseq_tstart + i + 1
                    if self.strand == "+":
                        qpos = alignseq_qstart + i + 1
                    else:
                        qpos = alignseq_qend - i
                    msg += (
                        f"{rname}\t{rpos}\t.\t{alignseq_t[i]}\t{alignseq_q[i]}"
                        f"\t.\tAUTO\tALN_SCORE={self.score};ALN_QUERY={qname};"
                        f"ALN_QPOS={qpos};ALN_STRAND={self.strand}\n"
                    )
        return msg.rstrip()

    def to_sam(self, targetref, queryref) -> None:
        msg = ""

        # these are all correct for positive or negative strand (positions 1-based):
        qname = self.query
        rname = self.target
        pos = self.tstart + 1
        qlen = self.qlen
        qstart = self.qstart
        chainid = self.id
        segmentid = 1

        # hard clipping at start of alignment and 1-based target starts/ends:
        if self.strand == "+":
            lefthardclip = qstart
            querystartpos = qstart + 1
            queryendpos = querystartpos
            flag = 0
        else:
            lefthardclip = qstart
            querystartpos = qlen - qstart
            queryendpos = querystartpos
            flag = 16

        intervals = sorted(self.ttree.all_intervals)
        for i, intvl in enumerate(intervals):
            # deal with indels first, because they are to the left of the matching segment:
            if i > 0:
                deltat = intvl.data[1]
                deltaq = intvl.data[2]

                # split the alignment if there are unaligned bases in both target and query
                if deltaq > 0 and deltat > 0:
                    # write alignment line, increment segment number, reset start positions
                    # of alignment and left hard clipping
                    if self.strand == "+":
                        righthardclip = qlen - queryendpos + 1
                        seq = queryref.fetch(
                            reference=qname,
                            start=querystartpos - 1,
                            end=queryendpos - 1,
                        ).upper()
                    else:
                        righthardclip = queryendpos
                        seq = queryref.fetch(
                            reference=qname,
                            start=queryendpos,
                            end=querystartpos,
                        ).upper()
                        seq = reverse_complement(seq)
                    fullcigar = (
                        f"{lefthardclip}H"
                        + self.cigarstring
                        + f"{righthardclip}H"
                    )
                    msg += (
                        f"{qname}.{chainid}.{segmentid}\t{flag}\t{rname}\t{pos}\t"
                        f"0\t{fullcigar}\t*\t0\t0\t{seq}\t*\tNM:i:{self.nm_tag_val}\n"
                    )
                    pos = intvl.begin + 1
                    self.cigarstring = ""
                    self.nm_tag_val = 0
                    if self.strand == "+":
                        lefthardclip = intvl.begin + intvl.data[0]
                        querystartpos = lefthardclip + 1
                    else:
                        lefthardclip = qlen - intvl.begin - intvl.data[0]
                        querystartpos = intvl.begin + intvl.data[0]
                    queryendpos = querystartpos
                    segmentid += 1
                else:
                    queryendpos = self.update_cigar_indel(
                        deltaq=deltaq, deltat=deltat, queryendpos=queryendpos
                    )

            # update matched segment
            queryendpos = self.update_cigar_match(
                intvl=intvl,
                queryendpos=queryendpos,
                queryref=queryref,
                qname=qname,
                targetref=targetref,
                tname=rname,
            )

        # write last alignment line
        if self.strand == "+":
            righthardclip = qlen - queryendpos + 1
            seq = queryref.fetch(
                reference=qname, start=querystartpos - 1, end=queryendpos - 1
            ).upper()
        else:
            righthardclip = queryendpos
            seq = queryref.fetch(
                reference=qname, start=queryendpos, end=querystartpos
            ).upper()
            seq = reverse_complement(seq)
        fullcigar = f"{lefthardclip}H" + self.cigarstring + f"{righthardclip}H"
        msg += (
            f"{qname}.{chainid}.{segmentid}\t{flag}\t{rname}\t{pos}"
            f"\t0\t{fullcigar}\t*\t0\t0\t{seq}\t*\tNM:i:{self.nm_tag_val}"
        )

        return msg


# def break_chain(chain: Chain, break_unaligned: bool=True, break_indel: int=0) -> list:
#     list_chains = []
#     return list_chains


def vcf_header(dict_contig_length) -> str:
    header = ""
    header += f"##fileformat=VCFv4.3\n"
    header += f'##FILTER=<ID=AUTO,Description="Generated automatically.">\n'

    for contig in sorted(
        dict_contig_length,
        key=lambda i: int(dict_contig_length[i]),
        reverse=True,
    ):
        header += (
            f"##contig=<ID={contig}, " f"length={dict_contig_length[contig]}>\n"
        )
    header += f'##INFO=<ID=ALN_SCORE,Number=A,Type=Integer,Description="Score of chain alignment.">\n'
    header += f'##INFO=<ID=ALN_QUERY,Number=A,Type=String,Description="Variant query contig.">\n'
    header += f'##INFO=<ID=ALN_QPOS,Number=A,Type=Integer,Description="Variant position on query contig.">\n'
    header += f'##INFO=<ID=ALN_STRAND,Number=A,Type=String,Description="Strand of chain alignment.">\n'
    header += f'##INFO=<ID=ALN_DT,Number=A,Type=Integer,Description="Length of gap on target sequence.">\n'
    header += f'##INFO=<ID=ALN_DQ,Number=A,Type=Integer,Description="Length of gap on query sequence.">\n'
    header += (
        f"#CHROM       POS     ID      REF     ALT     QUAL    FILTER  INFO\n"
    )

    return header.rstrip()


def sam_header(dict_contig_length) -> str:
    header = ""
    header += f"@HD\tVN:1.6\tSO:unsorted\n"

    for contig in sorted(
        dict_contig_length,
        key=lambda i: int(dict_contig_length[i]),
        reverse=True,
    ):
        header += f"@SQ\tSN:{contig}\tLN:{dict_contig_length[contig]}\n"

    return header.rstrip()
