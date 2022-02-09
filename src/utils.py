'''
Utils for chain-related processing

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import argparse
import copy
import intervaltree
import pysam
import sys

'''
Read a FASTA file as a dict if a file name is given. If not, return an empty dict.
'''
def read_fasta(ref_fn: str) -> dict:
    ref = {}
    if ref_fn != '':
        f = pysam.FastaFile(ref_fn)
        for r in f.references:
            ref[r] = f[r].upper()
    return ref


def reverse_complement(seq: str) -> str:
    d = {'A': 'T', 'a': 'T', 'C': 'G', 'c': 'G',
         'G': 'C', 'g': 'G', 'T': 'A', 't': 'A',
         'N': 'N'}
    rc = ''
    for s in seq:
        if s in d:
            rc += d[s]
        else:
            print(f'Base "{s}" is not a known nucleotide and is converted to N',
                  file=sys.stderr)
            rc += 'N'
    return rc[::-1]


def compute_hamming_dist(
    forward: bool,
    ref1: dict, contig1: str, start1: int, end1: int,
    ref2: dict, contig2: str, start2: int, end2: int
) -> float:
    if contig1 in ref1:
        s1 = ref1[contig1][start1: end1]
    else:
        print(f'Warning : {contig1} not in ref1', file=sys.stderr)
        return 0
    if contig2 in ref2:
        s2 = ref2[contig2][start2: end2]
        if not forward:
            s2 = reverse_complement(s2)
    else:
        print(f'Warning : {contig2} not in ref2', file=sys.stderr)
        return 0

    try:
        assert (len(s1) == len(s2))
    except:
        print('Error: lengths do not match', file=sys.stderr)
        print(len(s1), len(s2), file=sys.stderr)
        print(f'{contig1}:{start1}-{end1}', file=sys.stderr)
        print(f'{contig2}:{start2}-{end2}', file=sys.stderr)
        exit(1)

    if len(s1) == 0:
        return 0
    idy = 0
    for i_s, s in enumerate(s1):
        if s == s2[i_s]:
            idy += 1
    idy /= len(s1)
    return idy


class ChainConst():
    CHDR_SCORE = 1
    CHDR_SOURCE = 2
    CHDR_SLEN = 3
    CHDR_SSTRAND = 4
    CHDR_SSTART = 5
    CHDR_SEND = 6
    CHDR_TARGET = 7
    CHDR_TLEN = 8
    CHDR_TSTRAND = 9
    CHDR_TSTART = 10
    CHDR_TEND = 11
    CHDR_ID = 12
    C_SIZE = 0
    C_DS = 1
    C_DT = 2


class Chain(ChainConst):
    def __init__(self, fields=None):
        self.tr = intervaltree.IntervalTree()
        if not fields:
            return
        self.score = int(fields[self.CHDR_SCORE])
        self.source = fields[self.CHDR_SOURCE]
        self.slen = int(fields[self.CHDR_SLEN])
        self.sstart = int(fields[self.CHDR_SSTART])
        self.soffset = self.sstart
        self.send = int(fields[self.CHDR_SEND])
        self.sstrand = fields[self.CHDR_SSTRAND]
        # assert fields[self.CHDR_SSTRAND] == '+'
        # if fields[self.CHDR_SSTRAND] != '+':
        #     print(fields)

        self.target = fields[self.CHDR_TARGET]
        self.tlen = int(fields[self.CHDR_TLEN])
        self.tstrand = fields[self.CHDR_TSTRAND]
        # assert fields[self.CHDR_TSTRAND] == '+'
        # if fields[self.CHDR_TSTRAND] != '+':
        #     print(fields)
        self.tstart = int(fields[self.CHDR_TSTART])
        if self.tstrand == '+':
            self.toffset = self.tstart
        else:
            self.toffset = self.tlen - self.tstart
        self.tend = int(fields[self.CHDR_TEND])

        self.strand = '+' if (self.sstrand == '+' and self.tstrand == '+') \
                          else '-'
        self.id = fields[self.CHDR_ID]
        self.ds = 0
        self.dt = 0

    def add_record_three(self, fields):
        segment_size = int(fields[self.C_SIZE])
        if segment_size > 0:
            self.tr[self.soffset : self.soffset + segment_size] = \
                (self.toffset - self.soffset, self.ds, self.dt, False)
        self.ds = int(fields[self.C_DS])
        self.dt = int(fields[self.C_DT])
        self.soffset += (segment_size + self.ds)
        if self.strand == '+':
            self.toffset += (segment_size + self.dt)
        else:
            self.toffset -= (segment_size + self.dt)

    def add_record_one(self, fields):
        segment_size = int(fields[self.C_SIZE])
        if segment_size > 0:
            self.tr[self.soffset : self.soffset + segment_size] = \
                (self.toffset-self.soffset, self.ds, self.dt)

    def print_chain(self) -> str:
        msg = (f'chain {self.score} '
               f'{self.source} {self.slen} {self.sstrand} {self.sstart} {self.send} '
               f'{self.target} {self.tlen} {self.tstrand} {self.tstart} {self.tend} '
               f'{self.id}\n')
        intervals = sorted(self.tr.all_intervals)
        for i, intvl in enumerate(intervals):
            if i == 0:
                if intvl[2][1:3] != (0, 0):
                    msg += (f'{intervals[0].begin - self.sstart - intervals[0].data[1]}\t'
                            f'{intervals[0].data[1]}\t{intervals[0].data[2]}\n')
            else:
                pintvl = intervals[i-1]
                msg += (f'{pintvl.end - pintvl.begin}\t{intvl[2][1]}\t{intvl[2][2]}\n')

        if intvl.end == self.send:
        # if intvl.end == self.send and self.send + intvl.data[0] == self.tend:
            msg += (f'{intvl.end - intvl.begin}\n\n')
        else:
            msg += (f'{intvl.end - intvl.begin}\t'
                    f'{self.send - intvl.end}\t'
                    f'{self.tend - (intvl.end+intvl.data[0])}\n0\n')
        return msg

    def try_merge(self, c):
        # print(len(list(self.tr)))
        # print(self.tr.items())
        # print(c.tr.items())

        self_sorted_intervals = sorted(self.tr.all_intervals)
        c_sorted_intervals = sorted(c.tr.all_intervals)
        
        clone_tr = self.tr.copy()
        for c_intvl in c_sorted_intervals:
            for i_s, s_intvl in enumerate(self_sorted_intervals):
                if c_intvl.begin > s_intvl.end or c_intvl.end < s_intvl.begin:
                    continue
                else:
                    # If offsets are the same, good to merge
                    if c_intvl.data[0] == s_intvl.data[0]:
                        if c_intvl.begin < s_intvl.begin and c_intvl.end < s_intvl.end:
                            # Extend from beginning
                            # print(c_intvl)
                            # print(s_intvl)
                            start = min(s_intvl.begin, c_intvl.begin)
                            diff = s_intvl.begin - c_intvl.begin
                            data = (s_intvl.data[0], s_intvl.data[1]-diff, s_intvl.data[2]-diff)
                            updated_s_intvl = intervaltree.Interval(start, s_intvl.end, data)
                            clone_tr.remove(s_intvl)
                            c.tr.remove(c_intvl)
                            clone_tr.add(updated_s_intvl)
                        elif c_intvl.begin > s_intvl.begin and c_intvl.end > s_intvl.end:
                            # Extend from end
                            # This operation is trickier b/c the chain diffs are stored in the
                            # next interval. We need to pop the next interval and update it.
                            # print(c_intvl)
                            # print(s_intvl)
                            clone_tr.remove(s_intvl)
                            c.tr.remove(c_intvl)
                            end = max(s_intvl.end, c_intvl.end)
                            diff = c_intvl.end - s_intvl.end
                            clone_tr.add(intervaltree.Interval(s_intvl.begin, end, s_intvl[2]))
                            next_s_intvl = self_sorted_intervals[i_s+1]
                            next_data = (
                                next_s_intvl.data[0], next_s_intvl.data[1]-diff,
                                next_s_intvl.data[2]-diff)
                            updated_next_s_intvl = intervaltree.Interval(
                                next_s_intvl.begin, next_s_intvl.end, next_data)
                            clone_tr.remove(next_s_intvl)
                            clone_tr.add(updated_next_s_intvl)
                        elif c_intvl.begin < s_intvl.begina and c_intvl.end > s_intvl.end:
                            print('During merging, the new interval consumes the existing one',
                                  file=sys.stderr)
                            print('This operation is not yet supported', file=sys.stderr)
                            return False
                    else:
                        return False
        self.tr = clone_tr
        # print('\nmerged:')
        # print(self.tr.items())
        sorted_intervals = sorted(self.tr.all_intervals)

        # Update chain info
        self.score += c.score
        if sorted_intervals[0][0] < self.sstart:
            self.sstart = sorted_intervals[0][0]
            self.tstart = self.sstart + sorted_intervals[0][2]
        if sorted_intervals[-1][1] > self.send:
            self.send = sorted_intervals[-1][1]
            self.tend = self.send + sorted_intervals[-1][2]

        return True

    def to_paf(self, fn_paf) -> None:
        pass

