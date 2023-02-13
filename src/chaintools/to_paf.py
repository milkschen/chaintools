#!/usr/bin/env python3
'''
Convert a chain file to the PAF format

Nae-Chyun Chen
Johns Hopkins University
2022
'''
import argparse
import sys
from typing import TextIO

import pysam

from chaintools import utils


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c',
                        '--chain',
                        default='-',
                        help='Path to the chain file')
    parser.add_argument('-t',
                        '--targetfasta',
                        default='',
                        help=('Path to the fasta file for the target '
                              '(reference) genome of the chain file'))
    parser.add_argument('-q',
                        '--queryfasta',
                        default='',
                        help=('Path to the fasta file for the query '
                              'genome of the chain file'))
    parser.add_argument('-o',
                        '--output',
                        default='',
                        help='Path to the output PAF file.')
    args = parser.parse_args()
    return args


def write_to_paf(f: TextIO,
                 targetref: pysam.FastaFile = None,
                 queryref: pysam.FastaFile = None) -> str:
    for line in f:
        fields = line.split()
        if len(fields) == 0:
            pass
        elif fields[0] == 'chain':
            c = utils.Chain(fields)
        else:
            c.add_record(fields)
            if len(fields) == 1:
                yield c.to_paf(targetref=targetref, queryref=queryref)
                c = None


def write_to_paf_io(fn_chain: str,
                    fn_paf: str,
                    fn_targetfasta: str = '',
                    fn_queryfasta: str = '') -> None:
    if fn_chain == '-':
        f = sys.stdin
    else:
        f = open(fn_chain, 'r')
    if fn_paf:
        fo = open(fn_paf, 'w')
    else:
        fo = sys.stdout

    targetref = utils.fasta_reader(fn_targetfasta)
    queryref = utils.fasta_reader(fn_queryfasta)

    out = write_to_paf(f=f, targetref=targetref, queryref=queryref)
    while True:
        try:
            print(next(out), file=fo)
        except StopIteration:
            break


def main(argv=sys.argv):
    args = parse_args()
    write_to_paf_io(fn_chain=args.chain,
                    fn_paf=args.output,
                    fn_targetfasta=args.targetfasta,
                    fn_queryfasta=args.queryfasta)


if __name__ == '__main__':
    main()
