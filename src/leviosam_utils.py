'''
Utils for levioSAM

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import pysam

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


def reverse_complement(seq):
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

