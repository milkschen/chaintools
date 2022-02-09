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
        '-o', '--output', default='',
        help='Path to the output PAF file.'
    )
    args = parser.parse_args()
    return args

def write_to_paf(fn_chain: str, fn_paf: str):
    f = open(fn_chain, 'r')
    if fn_paf:
        fo = open(fn_paf, 'w')
    else:
        fo = sys.stdout

    for line in f:
        fields = line.split()
        if len(fields) == 0:
            continue
        elif line.startswith('chain'):
            c = utils.Chain(fields)
        elif len(fields) == 3:
            c.add_record_three(fields)
            pass
        elif len(fields) == 1:
            c.add_record_one(fields)
            print(c.to_paf(), file=fo)
            c = None


if __name__ == '__main__':
    args = parse_args()
    write_to_paf(fn_chain=args.chain, fn_paf=args.output)
