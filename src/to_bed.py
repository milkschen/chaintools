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
        '-c', '--chain', default='-',
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


def write_to_bed(fn_chain: str, fn_bed: str, coord: str):
    if fn_chain == '-':
        f = sys.stdin
    else:
        f = open(fn_chain, 'r')
    if fn_bed:
        fo = open(fn_bed, 'w')
    else:
        fo = sys.stdout

    for line in f:
        fields = line.split()
        if len(fields) == 0:
            continue
        elif line.startswith('chain'):
            c = utils.Chain(fields)
        else:
            bed_str = c.record_to_bed(fields=fields, coord=coord)
            if bed_str != '':
                print(bed_str, file=fo)
            c.add_record(fields)


if __name__ == '__main__':
    args = parse_args()
    if args.coord not in ['target', 'query']:
        raise(ValueError, 'Illegal `-coord` value. Should be in ["target", "query"]')

    write_to_bed(fn_chain=args.chain, fn_bed=args.output, coord=args.coord)
