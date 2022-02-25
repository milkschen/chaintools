'''
Utils for chain-related processing

Nae-Chyun Chen
Johns Hopkins University

Nancy Fisher Hansen
NIH/NHGRI

2021-2022
'''
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

def get_target_entries(fn: str)->dict:
    CC = ChainConst()
    dict_contig_length = {}
    with open(fn, 'r') as f:
        for line in f:
            line = line.split()
            if len(line) == 13:
                dict_contig_length[line[CC.CHDR_TARGET]] = line[CC.CHDR_TLEN]
    return dict_contig_length


class ChainConst():
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
        if self.qstrand == '+':
            self.qoffset = self.qstart
        else:
            self.qoffset = self.qlen - self.qstart
        self.qend = int(fields[self.CHDR_QEND])

        self.strand = '+' if (self.tstrand == '+' and self.qstrand == '+') \
                          else '-'
        self.id = fields[self.CHDR_ID]
        self.ds = 0
        self.dq = 0

        self.seglen = 0

    def add_record_three(self, fields):
        segment_size = int(fields[self.C_SIZE])
        if segment_size > 0:
            self.ttree[self.toffset : self.toffset + segment_size] = \
                (self.qoffset - self.toffset, self.ds, self.dq, False)
            if self.strand == '+':
                self.qtree[self.qoffset: self.qoffset + segment_size] = 1
            else:
                self.qtree[self.qoffset - segment_size: self.qoffset] = 1

        self.ds = int(fields[self.C_DT])
        self.dq = int(fields[self.C_DQ])
        self.toffset += (segment_size + self.ds)
        self.seglen += segment_size
        if self.strand == '+':
            self.qoffset += (segment_size + self.dq)
        else:
            self.qoffset -= (segment_size + self.dq)

    def add_record_one(self, fields):
        segment_size = int(fields[self.C_SIZE])
        self.seglen += segment_size
        if segment_size > 0:
            self.ttree[self.toffset : self.toffset + segment_size] = \
                (self.qoffset-self.toffset, self.ds, self.dq)
            if self.strand == '+':
                self.qtree[self.qoffset: self.qoffset + segment_size] = 1
            else:
                self.qtree[self.qoffset - segment_size: self.qoffset] = 1

    def print_chain(self) -> str:
        msg = (f'chain {self.score} '
               f'{self.target} {self.tlen} {self.tstrand} {self.tstart} {self.tend} '
               f'{self.query} {self.qlen} {self.qstrand} {self.qstart} {self.qend} '
               f'{self.id}\n')
        intervals = sorted(self.ttree.all_intervals)
        for i, intvl in enumerate(intervals):
            if i == 0:
                # If the size of the first interval is not zero
                if intvl.data[1:3] != (0, 0):
                    msg += (f'{intvl.begin - self.tstart - intvl.data[1]}\t'
                            f'{intvl.data[1]}\t{intvl.data[2]}\n')
            else:
                pintvl = intervals[i-1]
                msg += (f'{pintvl.end - pintvl.begin}\t{intvl.data[1]}\t{intvl.data[2]}\n')

        if intvl.end == self.tend:
            msg += (f'{intvl.end - intvl.begin}\n')
        else:
            msg += (f'{intvl.end - intvl.begin}\t'
                    f'{self.tend - intvl.end}\t'
                    f'{self.qend - (intvl.end+intvl.data[0])}\n0')
        return msg

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

    def to_paf(self) -> None:
        def update_cigar(msg, num_m, dt, dq) -> str:
            if dt != 0 and dq != 0:
                shared = min(dt, dq)
                num_m += shared
                dt -= shared
                dq -= shared
            if num_m > 0:
                msg += f'{num_m}M'
            if max(dt, dq) > 0:
                if dt > dq:
                    msg += f'{dt - dq}D'
                elif dt < dq:
                    msg += f'{dq - dt}I'
            return msg

        # aln_block_len = 0
        num_match = 0
        total_block_t = 0
        # total_block_q = 0
        msg = ''
        intervals = sorted(self.ttree.all_intervals)
        for i, intvl in enumerate(intervals):
            if i == 0:
                num_m = intvl.begin - self.tstart - intvl.data[1]
            else:
                pintvl = intervals[i-1]
                num_m = pintvl.end - pintvl.begin
            dt = intvl.data[1]
            dq = intvl.data[2]
            num_match += num_m
            # aln_block_len += (num_m + max([dt, dq]))
            total_block_t += (num_m + dt)
            # total_block_q += (num_m + dq)
            msg = update_cigar(msg, num_m, dt, dq)

        if intvl.end == self.tend:
            msg += f'{intvl.end - intvl.begin}M'
        else:
            num_m = intvl.end - intvl.begin
            dt = self.tend - intvl.end
            dq = self.qend - (intvl.end+intvl.data[0])
            num_match += num_m
            # aln_block_len += (num_m + max([dt, dq]))
            total_block_t += (num_m + dt)
            # total_block_q += (num_m + dq)
            msg = update_cigar(msg, num_m, dt, dq)

        if self.qstrand == '+':
            hdr = (f'{self.query}\t{self.qlen}\t{self.qstart}\t{self.qend}\t{self.qstrand}\t'
                   f'{self.target}\t{self.tlen}\t{self.tstart}\t{self.tend}\t{num_match}\t'
                   f'{total_block_t}\t255\tcg:Z:')
        else:
            hdr = (f'{self.query}\t{self.qlen}\t{self.qlen - self.qend}\t{self.qlen - self.qstart}\t{self.qstrand}\t'
                   f'{self.target}\t{self.tlen}\t{self.tstart}\t{self.tend}\t{num_match}\t'
                   f'{total_block_t}\t255\tcg:Z:')

        return hdr+msg
    
    def to_vcf(self, targetref, queryref) -> None:
        msg = ''
        rname = self.target
        qname = self.query

        intervals = sorted(self.ttree.all_intervals)
        for i, intvl in enumerate(intervals):
            # deal with indels first, because they are to the left of the matching segment:
            if i > 0:
                # calculate zero-based, half-open ref and alt start/end (subtracting one to add base prior to event)
                if self.strand == "+":
                    alignseq_qend = intvl.begin + intvl.data[0]
                    alignseq_qstart = intvl.begin + intvl.data[0] - intvl.data[2] - 1
                    alignseq_q = queryref.fetch(reference=qname, start=alignseq_qstart, end=alignseq_qend).upper()
                    qpos = alignseq_qstart + 1
                else:
                    alignseq_qstart = intvl.begin + intvl.data[0]
                    alignseq_qend = intvl.begin + intvl.data[0] + intvl.data[2] + 1
                    alignseq_q = queryref.fetch(reference=qname, start=alignseq_qstart, end=alignseq_qend).upper()
                    alignseq_q = reverse_complement(alignseq_q)
                    qpos = alignseq_qend

                # end position of indel allele is position to left of matching segment's start, start has a base appended to the left:
                alignseq_tend = intvl.begin
                alignseq_tstart = intvl.begin - intvl.data[1] - 1
                alignseq_t = targetref.fetch(reference=rname, start=alignseq_tstart, end=alignseq_tend).upper()
                rpos = alignseq_tstart + 1

                msg += (f'{rname}\t{rpos}\t.\t{alignseq_t}\t{alignseq_q}\t.\tAUTO\tALN_SCORE={self.score};' +
                        f'ALN_QUERY={qname};ALN_QPOS={qpos};ALN_STRAND={self.strand};ALN_DQ={intvl.data[2]};ALN_DT={intvl.data[1]}\n')

            # now print SNPs within matched segment:

            # zero-based, half-open start/end of matched segment
            if self.strand == "+":
                alignseq_qstart = intvl.begin + intvl.data[0]
                alignseq_qend = intvl.end + intvl.data[0]
                alignseq_q = queryref.fetch(reference=qname, start=alignseq_qstart, end=alignseq_qend).upper()
            else:
                alignseq_qend = intvl.begin + intvl.data[0]
                alignseq_qstart = alignseq_qend - (intvl.end - intvl.begin)
                alignseq_q = queryref.fetch(reference=qname, start=alignseq_qstart, end=alignseq_qend).upper()
                alignseq_q = reverse_complement(alignseq_q)

            alignseq_tstart = intvl.begin
            alignseq_tend = intvl.end
            alignseq_t = targetref.fetch(reference=rname, start=alignseq_tstart, end=alignseq_tend).upper()

            segmentlength = intvl.end - intvl.begin

            #print(f'Tstart {alignseq_tstart} Tend {alignseq_tend} Sstart {alignseq_sstart} Send {alignseq_send} Intvl {intvl.begin}:{intvl.end} Toffset {intvl.data[1]} Seglength {segmentlength}')

            for i in range(segmentlength):
                if alignseq_q[i] != alignseq_t[i]:
                    # convert to one-based start position:
                    rpos = alignseq_tstart + i + 1
                    if self.strand == "+":
                        qpos = alignseq_qstart + i + 1
                    else:
                        qpos = alignseq_qend - i
                    msg += (f'{rname}\t{rpos}\t.\t{alignseq_t[i]}\t{alignseq_q[i]}\t.\tAUTO\tALN_SCORE={self.score};ALN_QUERY={qname};ALN_QPOS={qpos};ALN_STRAND={self.strand}\n')

        return msg
    
    def to_sam(self, targetref, queryref) -> None:
        msg = ''
        cigarstring = ''
        nmtagval = 0

        # these are all correct for positive or negative strand (positions 1-based):
        qname = self.query
        rname = self.target
        pos = self.tstart + 1
        qlen = self.qlen
        qstart = self.qstart
        chainid = self.id
        segmentid = 1

        # hard clipping at start of alignment and 1-based target starts/ends:
        if self.strand == '+':
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
                    # write alignment line, increment segment number, reset start positions of alignment and left hard clipping
                    if self.strand == '+':
                        righthardclip = qlen - queryendpos + 1
                        seq = queryref.fetch(reference=qname, start=querystartpos-1, end=queryendpos-1).upper()
                    else:
                        righthardclip = queryendpos
                        seq = queryref.fetch(reference=qname, start=queryendpos, end=querystartpos).upper()
                        seq = reverse_complement(seq)
                    fullcigar = f'{lefthardclip}H' + cigarstring + f'{righthardclip}H'
                    msg += (f'{qname}.{chainid}.{segmentid}\t{flag}\t{rname}\t{pos}\t0\t{fullcigar}\t*\t0\t0\t{seq}\t*\tNM:i:{nmtagval}\n')
                    pos = intvl.begin + 1
                    cigarstring = ''
                    nmtagval = 0
                    if self.strand == '+':
                        lefthardclip = intvl.begin + intvl.data[0]
                        querystartpos = lefthardclip + 1
                    else:
                        lefthardclip = qlen - intvl.begin - intvl.data[0]
                        querystartpos = intvl.begin + intvl.data[0]
                    queryendpos = querystartpos
                    segmentid += 1
                else:
                    if deltaq > 0:
                        cigarstring += f'{deltaq}I'
                        nmtagval += deltaq
                    else:
                        cigarstring += f'{deltat}D'
                        nmtagval += deltat
                    if self.strand == '+':
                        queryendpos += deltaq
                    else:
                        queryendpos -= deltaq

            # now add X/= counts to cigar within matched segment:
            # zero-based, half-open start/end of matched segment

            segmentlength = intvl.end - intvl.begin
            if self.strand == "+":
                alignseq_qstart = intvl.begin + intvl.data[0]
                alignseq_qend = intvl.end + intvl.data[0]
                alignseq_q = queryref.fetch(reference=qname, start=alignseq_qstart, end=alignseq_qend).upper()
                queryendpos += segmentlength
            else:
                alignseq_qend = intvl.begin + intvl.data[0]
                alignseq_qstart = alignseq_qend - (intvl.end - intvl.begin)
                alignseq_q = queryref.fetch(reference=qname, start=alignseq_qstart, end=alignseq_qend).upper()
                alignseq_q = reverse_complement(alignseq_q)
                queryendpos -= segmentlength

            alignseq_tstart = intvl.begin
            alignseq_tend = intvl.end
            alignseq_t = targetref.fetch(reference=rname, start=alignseq_tstart, end=alignseq_tend).upper()

            #print(f'Tstart {alignseq_tstart} Tend {alignseq_tend} Sstart {alignseq_sstart} Send {alignseq_send} Intvl {intvl.begin}:{intvl.end} Toffset {intvl.data[1]} Seglength {segmentlength}')

            nummatches = 0
            nummismatches = 0
            for i in range(segmentlength):
                if alignseq_q[i] == alignseq_t[i]:
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
            righthardclip = qlen - queryendpos + 1
            seq = queryref.fetch(reference=qname, start=querystartpos-1, end=queryendpos-1).upper()
        else:
            righthardclip = queryendpos
            seq = queryref.fetch(reference=qname, start=queryendpos, end=querystartpos).upper()
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
    header += (f'##INFO=<ID=ALN_QUERY,Number=A,Type=String,Description="Variant query contig.">\n')
    header += (f'##INFO=<ID=ALN_QPOS,Number=A,Type=Integer,Description="Variant position on query contig.">\n')
    header += (f'##INFO=<ID=ALN_STRAND,Number=A,Type=String,Description="Strand of chain alignment.">\n')
    header += (f'##INFO=<ID=ALN_DT,Number=A,Type=Integer,Description="Length of gap on target sequence.">\n')
    header += (f'##INFO=<ID=ALN_DQ,Number=A,Type=Integer,Description="Length of gap on query sequence.">\n')
    header += (f'#CHROM       POS     ID      REF     ALT     QUAL    FILTER  INFO\n')

    return header

def sam_header(dict_contig_length) -> str:
    header = ''
    header += (f'@HD\tVN:1.6\tSO:unsorted\n')

    for contig in sorted(dict_contig_length, key=lambda i: int(dict_contig_length[i]), reverse=True):
        header += (f'@SQ\tSN:{contig}\tLN:{dict_contig_length[contig]}\n')

    return header
