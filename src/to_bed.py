'''
Convert a chain file to the BED format

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
        '-c', '--chain', required=True,
        help='Path to the chain file'
    )
    parser.add_argument(
        '-o', '--output', default='',
        help='Path to the output BED file.'
    )
    parser.add_argument(
        '--coord', default='target',
        help='Output coordinate system. Options: target, query. ["target"]'
    )
    args = parser.parse_args()
    return args


def to_bed(fn_chain: str, fn_paf: str, coord: str):
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
            segment_size = int(fields[0])
            if segment_size > 0:
                if coord == 'target':
                    print(f'{c.target}\t{c.toffset}\t{c.toffset + segment_size}', file=fo)
                elif coord == 'query':
                    if c.strand == '+':
                        print(f'{c.query}\t{c.qoffset}\t{c.qoffset + segment_size}', file=fo)
                    else:
                        print(f'{c.query}\t{c.qoffset - segment_size}\t{c.qoffset}', file=fo)
            c.add_record_three(fields)
        elif len(fields) == 1:
            segment_size = int(fields[0])
            if segment_size > 0:
                if coord == 'target':
                    print(f'{c.target}\t{c.toffset}\t{c.toffset + segment_size}', file=fo)
                elif coord == 'query':
                    if c.strand == '+':
                        print(f'{c.query}\t{c.qoffset}\t{c.qoffset + segment_size}', file=fo)
                    else:
                        print(f'{c.query}\t{c.qoffset - segment_size}\t{c.qoffset}', file=fo)
            c.add_record_one(fields)
            c = None


if __name__ == '__main__':
    args = parse_args()
    if args.coord not in ['target', 'query']:
        raise(ValueError, 'Illegal `-coord` value. Should be in ["target", "query"]')

    to_bed(fn_chain=args.chain, fn_paf=args.output, coord=args.coord)
