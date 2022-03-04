'''
Convert a chain file to the SAM format

Nancy Fisher Hansen
NIH/NHGRI

Nae-Chyun Chen
Johns Hopkins University

2021-2022
'''
import argparse
import utils
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--chain', default='',
        help='Path to the chain file'
    )
    parser.add_argument(
        '-t', '--targetfasta', required=True,
        help='Path to the fasta file for the target (reference) genome of the chain file'
    )
    parser.add_argument(
        '-q', '--queryfasta', required=True,
        help='Path to the fasta file for the query genome of the chain file'
    )
    parser.add_argument(
        '-o', '--output', default='',
        help='Path to the output SAM-formatted file.'
    )
    args = parser.parse_args()
    return args

def write_to_sam(fn_chain: str, fn_sam: str, fn_targetfasta: str, fn_queryfasta: str):
    if fn_chain == '-':
        f = sys.stdin
    else:
        f = open(fn_chain, 'r')
    if fn_sam:
        fo = open(fn_sam, 'w')
    else:
        fo = sys.stdout

    print(utils.sam_header(utils.get_target_entries(fn_chain)), file=fo, end='')

    targetref = utils.fasta_reader(fn_targetfasta)
    queryref = utils.fasta_reader(fn_queryfasta)

    for line in f:
        fields = line.split()
        if len(fields) == 0:
            continue
        elif line.startswith('chain'):
            c = utils.Chain(fields)
        else:
            c.add_record(fields)
            if len(fields) == 1:
                print(c.to_sam(targetref, queryref), file=fo, end='')
                c = None


if __name__ == '__main__':
    args = parse_args()
    write_to_sam(fn_chain=args.chain, fn_sam=args.output, fn_targetfasta=args.targetfasta, fn_queryfasta=args.queryfasta)
