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
        help='Path to the output merged chain file.'
    )
    args = parser.parse_args()
    return args

def filter(fn_chain: str, fn_out: str):
    # tr_dict = {}
    chains = []
    f = open(fn_chain, 'r')
    if fn_out:
        fo = open(fn_out, 'w')
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
            # chains.append(c)
            print(c.print_chain(), file=fo)
            # if c.source in tr_dict:
            #     for tr in tr_dict[c.source]:
            #         if not tr.try_merge(c):
            #             tr_dict[c.source].append(c)
            #             break
            # else:
            #     tr_dict[c.source] = [c]
            c = None

    # for c in chains:
    #     c.print_chain(fo)
    # for contig in tr_dict.keys():
    #     for tr in tr_dict[contig]:
    #         tr.print_chain(fo)


if __name__ == '__main__':
    args = parse_args()
    filter(fn_chain=args.chain, fn_out=args.output)
