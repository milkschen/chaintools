import argparse
import utils
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--chain', required=True,
        help='Path to the chain file'
    )
    parser.add_argument(
        '-t', '--targetfasta', required=True,
        help='Path to the fasta file for the target genome of the chain file'
    )
    parser.add_argument(
        '-s', '--sourcefasta', required=True,
        help='Path to the fasta file for the source genome of the chain file'
    )
    parser.add_argument(
        '-o', '--output', default='',
        help='Path to the output SAM-formatted file.'
    )
    args = parser.parse_args()
    return args

def write_to_sam(fn_chain: str, fn_sam: str, fn_targetfasta: str, fn_sourcefasta: str):
    f = open(fn_chain, 'r')
    if fn_sam:
        fo = open(fn_sam, 'w')
    else:
        fo = sys.stdout

    print(utils.sam_header(utils.get_source_entries(fn_chain)), file=fo, end='')

    targetref = utils.fasta_reader(fn_targetfasta)
    sourceref = utils.fasta_reader(fn_sourcefasta)

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
            print(c.to_sam(targetref, sourceref), file=fo, end='')
            c = None

if __name__ == '__main__':
    args = parse_args()
    write_to_sam(fn_chain=args.chain, fn_sam=args.output, fn_targetfasta=args.targetfasta, fn_sourcefasta=args.sourcefasta)
