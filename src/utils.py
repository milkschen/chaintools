'''
Utils for chain-related processing

Nae-Chyun Chen
Johns Hopkins University
2021-2022
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

def fasta_reader(ref_fn: str):
    f = pysam.FastaFile(ref_fn)
    return f


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

def get_source_entries(fn: str)->dict:
    CC = ChainConst()
    dict_contig_length = {}
    with open(fn, 'r') as f:
        for line in f:
            line = line.split()
            if len(line) == 13:
                dict_contig_length[line[CC.CHDR_SOURCE]] = line[CC.CHDR_SLEN]
    return dict_contig_length


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
        self.stree = intervaltree.IntervalTree()
        self.ttree = intervaltree.IntervalTree()
        if not fields:
            return
        self.score = int(fields[self.CHDR_SCORE])
        self.source = fields[self.CHDR_SOURCE]
        self.slen = int(fields[self.CHDR_SLEN])
        self.sstart = int(fields[self.CHDR_SSTART])
        self.soffset = self.sstart
        self.send = int(fields[self.CHDR_SEND])
        self.sstrand = fields[self.CHDR_SSTRAND]
        self.target = fields[self.CHDR_TARGET]
        self.tlen = int(fields[self.CHDR_TLEN])
        self.tstrand = fields[self.CHDR_TSTRAND]
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

        self.seglen = 0

    def add_record_three(self, fields):
        segment_size = int(fields[self.C_SIZE])
        if segment_size > 0:
            self.stree[self.soffset : self.soffset + segment_size] = \
                (self.toffset - self.soffset, self.ds, self.dt, False)
            if self.strand == '+':
                self.ttree[self.toffset: self.toffset + segment_size] = 1
            else:
                self.ttree[self.toffset - segment_size: self.toffset] = 1

        self.ds = int(fields[self.C_DS])
        self.dt = int(fields[self.C_DT])
        self.soffset += (segment_size + self.ds)
        self.seglen += segment_size
        if self.strand == '+':
            self.toffset += (segment_size + self.dt)
        else:
            self.toffset -= (segment_size + self.dt)

    def add_record_one(self, fields):
        segment_size = int(fields[self.C_SIZE])
        self.seglen += segment_size
        if segment_size > 0:
            self.stree[self.soffset : self.soffset + segment_size] = \
                (self.toffset-self.soffset, self.ds, self.dt)
            if self.strand == '+':
                self.ttree[self.toffset: self.toffset + segment_size] = 1
            else:
                self.ttree[self.toffset - segment_size: self.toffset] = 1

    def print_chain(self) -> str:
        msg = (f'chain {self.score} '
               f'{self.source} {self.slen} {self.sstrand} {self.sstart} {self.send} '
               f'{self.target} {self.tlen} {self.tstrand} {self.tstart} {self.tend} '
               f'{self.id}\n')
        intervals = sorted(self.stree.all_intervals)
        for i, intvl in enumerate(intervals):
            if i == 0:
                # If the size of the first interval is not zero
                if intvl.data[1:3] != (0, 0):
                    msg += (f'{intvl.begin - self.sstart - intvl.data[1]}\t'
                            f'{intvl.data[1]}\t{intvl.data[2]}\n')
            else:
                pintvl = intervals[i-1]
                msg += (f'{pintvl.end - pintvl.begin}\t{intvl.data[1]}\t{intvl.data[2]}\n')

        if intvl.end == self.send:
        # if intvl.end == self.send and self.send + intvl.data[0] == self.tend:
            msg += (f'{intvl.end - intvl.begin}\n')
        else:
            msg += (f'{intvl.end - intvl.begin}\t'
                    f'{self.send - intvl.end}\t'
                    f'{self.tend - (intvl.end+intvl.data[0])}\n0')
        return msg

    def try_merge(self, c):
        # print(len(list(self.stree)))
        # print(self.stree.items())
        # print(c.stree.items())

        self_sorted_intervals = sorted(self.stree.all_intervals)
        c_sorted_intervals = sorted(c.stree.all_intervals)
        
        clone_tr = self.stree.copy()
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
                            c.stree.remove(c_intvl)
                            clone_tr.add(updated_s_intvl)
                        elif c_intvl.begin > s_intvl.begin and c_intvl.end > s_intvl.end:
                            # Extend from end
                            # This operation is trickier b/c the chain diffs are stored in the
                            # next interval. We need to pop the next interval and update it.
                            # print(c_intvl)
                            # print(s_intvl)
                            clone_tr.remove(s_intvl)
                            c.stree.remove(c_intvl)
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
        self.stree = clone_tr
        # print('\nmerged:')
        # print(self.stree.items())
        sorted_intervals = sorted(self.stree.all_intervals)

        # Update chain info
        self.score += c.score
        if sorted_intervals[0][0] < self.sstart:
            self.sstart = sorted_intervals[0][0]
            self.tstart = self.sstart + sorted_intervals[0][2]
        if sorted_intervals[-1][1] > self.send:
            self.send = sorted_intervals[-1][1]
            self.tend = self.send + sorted_intervals[-1][2]

        return True


    def to_paf(self) -> None:
        def update_cigar(msg, num_m, ds, dt) -> str:
            if ds != 0 and dt != 0:
                shared = min(ds, dt)
                num_m += shared
                ds -= shared
                dt -= shared
            if num_m > 0:
                msg += f'{num_m}M'
            if max(ds, dt) > 0:
                if ds > dt:
                    msg += f'{ds - dt}I'
                elif ds < dt:
                    msg += f'{dt - ds}D'
            return msg

        if self.tstrand == '+':
            msg = (f'{self.source}\t{self.slen}\t{self.sstart}\t{self.send}\t{self.tstrand}\t'
                   f'{self.target}\t{self.tlen}\t{self.tstart}\t{self.tend}\t{self.score}\t'
                   f'{self.score}\t60\tcg:Z:')
        else:
            msg = (f'{self.source}\t{self.slen}\t{self.sstart}\t{self.send}\t{self.tstrand}\t'
                   f'{self.target}\t{self.tlen}\t{self.tlen - self.tend}\t{self.tlen - self.tstart}\t{self.score}\t'
                   f'{self.score}\t60\tcg:Z:')
        intervals = sorted(self.stree.all_intervals)
        for i, intvl in enumerate(intervals):
            if i == 0:
                num_m = intvl.begin - self.sstart - intvl.data[1]
            else:
                pintvl = intervals[i-1]
                num_m = pintvl.end - pintvl.begin
            ds = intvl.data[1]
            dt = intvl.data[2]
            msg = update_cigar(msg, num_m, ds, dt)

        if intvl.end == self.send:
            msg += f'{intvl.end - intvl.begin}M'
        else:
            num_m = intvl.end - intvl.begin
            ds = self.send - intvl.end
            dt = self.tend - (intvl.end+intvl.data[0])
            msg = update_cigar(msg, num_m, ds, dt)

        return msg
    
    def to_vcf(self, sourceref, targetref) -> None:
        msg = ''
        intervals = sorted(self.stree.all_intervals)
        for i, intvl in enumerate(intervals):
            # deal with indels first, because they are to the left of the matching segment:
            if i > 0:
                # calculate zero-based, half-open ref and alt start/end (subtracting one to add base prior to event)
                if self.strand == "+":
                    alignseq_tend = intvl.begin + intvl.data[0]
                    alignseq_tstart = intvl.begin + intvl.data[0] - intvl.data[2] - 1
                    alignseq_t = targetref.fetch(reference=self.target, start=alignseq_tstart, end=alignseq_tend).upper()
                    tpos = alignseq_tstart + 1
                else:
                    alignseq_tstart = intvl.begin + intvl.data[0]
                    alignseq_tend = intvl.begin + intvl.data[0] + intvl.data[2] + 1
                    alignseq_t = targetref.fetch(reference=self.target, start=alignseq_tstart, end=alignseq_tend).upper()
                    alignseq_t = reverse_complement(alignseq_t)
                    tpos = alignseq_tend

                # end position of indel allele is position to left of matching segment's start, start has a base appended to the left:
                alignseq_send = intvl.begin
                alignseq_sstart = intvl.begin - intvl.data[1] - 1
                alignseq_s = sourceref.fetch(reference=self.source, start=alignseq_sstart, end=alignseq_send).upper()
                spos = alignseq_sstart + 1

                msg += (f'{self.source}\t{spos}\t.\t{alignseq_s}\t{alignseq_t}\t.\tAUTO\tALN_SCORE={self.score};' +
                        f'ALN_TARGET={self.target};ALN_TPOS={tpos};ALN_STRAND={self.strand};ALN_DT={intvl.data[2]};ALN_DQ={intvl.data[1]}\n')

            # now print SNPs within matched segment:

            # zero-based, half-open start/end of matched segment
            if self.strand == "+":
                alignseq_tstart = intvl.begin + intvl.data[0]
                alignseq_tend = intvl.end + intvl.data[0]
                #alignseq_t = targetref[self.target][alignseq_tstart:alignseq_tend]
                alignseq_t = targetref.fetch(reference=self.target, start=alignseq_tstart, end=alignseq_tend).upper()
            else:
                alignseq_tend = intvl.begin + intvl.data[0]
                alignseq_tstart = alignseq_tend - (intvl.end - intvl.begin)
                #alignseq_t = targetref[self.target][alignseq_tstart:alignseq_tend]
                alignseq_t = targetref.fetch(reference=self.target, start=alignseq_tstart, end=alignseq_tend).upper()
                alignseq_t = reverse_complement(alignseq_t)

            alignseq_sstart = intvl.begin
            alignseq_send = intvl.end
            #alignseq_s = sourceref[self.source][alignseq_sstart:alignseq_send]
            alignseq_s = sourceref.fetch(reference=self.source, start=alignseq_sstart, end=alignseq_send).upper()

            segmentlength = intvl.end - intvl.begin

            #print(f'Tstart {alignseq_tstart} Tend {alignseq_tend} Sstart {alignseq_sstart} Send {alignseq_send} Intvl {intvl.begin}:{intvl.end} Toffset {intvl.data[1]} Seglength {segmentlength}')

            for i in range(segmentlength):
                if alignseq_t[i] != alignseq_s[i]:
                    # convert to one-based start position:
                    spos = alignseq_sstart + i + 1
                    if self.strand == "+":
                        tpos = alignseq_tstart + i + 1
                    else:
                        tpos = alignseq_tend - i
                    msg += (f'{self.source}\t{spos}\t.\t{alignseq_s[i]}\t{alignseq_t[i]}\t.\tAUTO\tALN_SCORE={self.score};ALN_TARGET={self.target};ALN_TPOS={tpos};ALN_STRAND={self.strand}\n')

        return msg
    
    def to_sam(self, sourceref, targetref) -> None:
        msg = ''
        cigarstring = ''
        nmtagval = 0

        # these are all correct for positive or negative strand (positions 1-based):
        qname = self.target
        rname = self.source
        pos = self.sstart + 1
        chainid = self.id
        segmentid = 1

        # hard clipping at start of alignment and 1-based target starts/ends:
        if self.strand == '+':
            lefthardclip = self.tstart
            targetstartpos = self.tstart + 1
            targetendpos = targetstartpos
            flag = 0
        else:
            lefthardclip = self.tstart
            targetstartpos = self.tlen - self.tstart
            targetendpos = targetstartpos
            flag = 16

        intervals = sorted(self.stree.all_intervals)
        for i, intvl in enumerate(intervals):
            # deal with indels first, because they are to the left of the matching segment:
            if i > 0:
                deltas = intvl.data[1]
                deltat = intvl.data[2]

                # split the alignment if there are unaligned bases in both the source and the target
                # might want to have an option to shut off this behavior
                if deltas > 0 and deltat > 0:
                    # write alignment line, increment segment number, reset start positions of alignment and left hard clipping
                    if self.strand == '+':
                        righthardclip = self.tlen - targetendpos + 1
                        seq = targetref.fetch(reference=qname, start=targetstartpos-1, end=targetendpos-1).upper()
                    else:
                        righthardclip = targetendpos
                        print(f'targetlimits {qname}:{targetendpos}-{targetstartpos} revcomp')
                        seq = targetref.fetch(reference=qname, start=targetendpos, end=targetstartpos).upper()
                        seq = reverse_complement(seq)
                    fullcigar = f'{lefthardclip}H' + cigarstring + f'{righthardclip}H'
                    msg += (f'{qname}.{chainid}.{segmentid}\t{flag}\t{rname}\t{pos}\t0\t{fullcigar}\t*\t0\t0\t{seq}\t*\tNM:i:{nmtagval}\n')
                    pos = intvl.begin + 1
                    cigarstring = ''
                    nmtagval = 0
                    if self.strand == '+':
                        lefthardclip = intvl.begin + intvl.data[0]
                        targetstartpos = lefthardclip + 1
                    else:
                        lefthardclip = self.tlen - intvl.begin - intvl.data[0]
                        targetstartpos = intvl.begin + intvl.data[0]
                    targetendpos = targetstartpos
                    segmentid += 1
                else:
                    if deltas > 0:
                        cigarstring += f'{deltas}D'
                        nmtagval += deltas
                    else:
                        cigarstring += f'{deltat}I'
                        nmtagval += deltat
                    if self.strand == '+':
                        targetendpos += deltat
                    else:
                        targetendpos -= deltat

            # now add X/= counts to cigar within matched segment:
            # zero-based, half-open start/end of matched segment

            segmentlength = intvl.end - intvl.begin
            if self.strand == "+":
                alignseq_tstart = intvl.begin + intvl.data[0]
                alignseq_tend = intvl.end + intvl.data[0]
                alignseq_t = targetref.fetch(reference=self.target, start=alignseq_tstart, end=alignseq_tend).upper()
                targetendpos += segmentlength
            else:
                alignseq_tend = intvl.begin + intvl.data[0]
                alignseq_tstart = alignseq_tend - (intvl.end - intvl.begin)
                alignseq_t = targetref.fetch(reference=self.target, start=alignseq_tstart, end=alignseq_tend).upper()
                alignseq_t = reverse_complement(alignseq_t)
                targetendpos -= segmentlength

            alignseq_sstart = intvl.begin
            alignseq_send = intvl.end
            alignseq_s = sourceref.fetch(reference=self.source, start=alignseq_sstart, end=alignseq_send).upper()

            #print(f'Tstart {alignseq_tstart} Tend {alignseq_tend} Sstart {alignseq_sstart} Send {alignseq_send} Intvl {intvl.begin}:{intvl.end} Toffset {intvl.data[1]} Seglength {segmentlength}')

            nummatches = 0
            nummismatches = 0
            for i in range(segmentlength):
                if alignseq_t[i] == alignseq_s[i]:
                    if nummismatches > 0:
                        cigarstring += f'{nummismatches}X'
                        nmtagval += nummismatches
                        nummismatches = 0
                    nummatches += 1
                else:
                    if nummatches > 0:
                        cigarstring += f'{nummatches}='
                        nummatches = 0
                    nummismatches += 1
            if nummatches > 0:
                cigarstring += f'{nummatches}='
            elif nummismatches > 0:
                cigarstring += f'{nummismatches}X'
                nmtagval += nummismatches

        # write last alignment line
        if self.strand == '+':
            righthardclip = self.tlen - targetendpos + 1
            seq = targetref.fetch(reference=qname, start=targetstartpos-1, end=targetendpos-1).upper()
        else:
            righthardclip = targetendpos
            seq = targetref.fetch(reference=qname, start=targetendpos, end=targetstartpos).upper()
            seq = reverse_complement(seq)
        fullcigar = f'{lefthardclip}H' + cigarstring + f'{righthardclip}H'
        msg += (f'{qname}.{chainid}.{segmentid}\t{flag}\t{rname}\t{pos}\t0\t{fullcigar}\t*\t0\t0\t{seq}\t*\tNM:i:{nmtagval}\n')

        return msg

def vcf_header(dict_contig_length) -> str:
    header = ''
    header += (f'##fileformat=VCFv4.3\n')
    header += (f'##FILTER=<ID=AUTO,Description="Generated automatically.">\n')

    for contig in sorted(dict_contig_length, key=lambda i: int(dict_contig_length[i]), reverse=True):
        header += (f'##contig=<ID={contig}, length={dict_contig_length[contig]}>\n')
    header += (f'##INFO=<ID=ALN_SCORE,Number=A,Type=Integer,Description="Score of chain alignment.">\n')
    header += (f'##INFO=<ID=ALN_TARGET,Number=A,Type=String,Description="Variant target contig.">\n')
    header += (f'##INFO=<ID=ALN_TPOS,Number=A,Type=Integer,Description="Variant position on target contig.">\n')
    header += (f'##INFO=<ID=ALN_STRAND,Number=A,Type=String,Description="Strand of chain alignment.">\n')
    header += (f'##INFO=<ID=ALN_DT,Number=A,Type=Integer,Description="Length of gap on SEQ_T.">\n')
    header += (f'##INFO=<ID=ALN_DQ,Number=A,Type=Integer,Description="Length of gap on SEQ_Q.">\n')
    header += (f'#CHROM       POS     ID      REF     ALT     QUAL    FILTER  INFO\n')

    return header

def sam_header(dict_contig_length) -> str:
    header = ''
    header += (f'@HD\tVN:1.6\tSO:unsorted\n')

    for contig in sorted(dict_contig_length, key=lambda i: int(dict_contig_length[i]), reverse=True):
        header += (f'@SQ\tSN:{contig}\tLN:{dict_contig_length[contig]}\n')

    return header
