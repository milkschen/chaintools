'''
Convert a chain file to the PAF format

Nae-Chyun Chen
Johns Hopkins University
2022
'''
import argparse
import utils
import sys


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--chain', default='-',
        help='Path to the chain file'
    )
    parser.add_argument(
        '-t', '--targetfasta', default='',
        help='Path to the fasta file for the target (reference) genome of the chain file'
    )
    parser.add_argument(
        '-q', '--queryfasta', default='',
        help='Path to the fasta file for the query genome of the chain file'
    )
    parser.add_argument(
        '-o', '--output', default='',
        help='Path to the output PAF file.'
    )
    args = parser.parse_args()
    return args


def write_to_paf(
    fn_chain: str, fn_paf: str,
    fn_targetfasta: str='', fn_queryfasta: str=''
) -> None:
    if fn_chain == '-':
        f = sys.stdin
    else:
        f = open(fn_chain, 'r')
    if fn_paf:
        fo = open(fn_paf, 'w')
    else:
        fo = sys.stdout

    if fn_targetfasta:
        targetref = utils.fasta_reader(fn_targetfasta)
    if fn_queryfasta:
        queryref = utils.fasta_reader(fn_queryfasta)

    for line in f:
        fields = line.split()
        if len(fields) == 0:
            continue
        elif line.startswith('chain'):
            c = utils.Chain(fields)
        elif len(fields) == 3:
            c.add_record_three(fields)
        elif len(fields) == 1:
            c.add_record_one(fields)
            print(c.to_paf(targetref=targetref, queryref=queryref),
                  file=fo)
            c = None


if __name__ == '__main__':
    args = parse_args()
    write_to_paf(
        fn_chain=args.chain, fn_paf=args.output,
        fn_targetfasta=args.targetfasta, fn_queryfasta=args.queryfasta)
