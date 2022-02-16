import argparse
import intervaltree
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
    parser.add_argument(
        '-u', '--unique', action='store_true',
        help='Activate to remove mappings that are not 1-1. Chains with smaller scores are excluded.'
    )
    parser.add_argument(
        '-oc', '--overlapped_chain', default='',
        help='Path to the chains that overlap with higher-scoring chains in `-c`. Leave this field empty to not write. Must with `-u`. [None]'
    )
    parser.add_argument(
        '-s', '--segment_size', default=0, type=int,
        help='Minimal segment size (the sum of chain segment sizes) allowed. [0]'
    )
    args = parser.parse_args()
    return args

def filter(
    fn_chain: str, fn_out: str, unique: bool, segment_size: int,
    fn_overlapped_chain: str=''
):
    stree_dict = {}
    ttree_dict = {}

    f = open(fn_chain, 'r')
    if fn_out:
        fo = open(fn_out, 'w')
    else:
        fo = sys.stdout
    if fn_overlapped_chain:
        foc = open(fn_overlapped_chain, 'w')

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
            
            # Filter by segment size
            if c.seglen < segment_size:
                continue

            print(
                f'score={c.score}\tsource={c.source}:{c.sstart}-{c.send}\ttarget={c.target}:{c.tstart}-{c.tend}\t{c.strand}',
                file=sys.stderr)

            s_overlap = False
            t_overlap = False
            if unique:
                # If two strees intersect
                if (c.source in stree_dict) and (len(c.stree & stree_dict[c.source]) > 0):
                    s_overlap = True
                    print('overlap (source)')
                    if fn_overlapped_chain:
                        print(c.print_chain(), file=foc)
                if (c.target in ttree_dict) and (len(c.ttree & ttree_dict[c.target]) > 0):
                    t_overlap = True
                    print('overlap (target)')
                    # print(c.ttree & ttree_dict[c.target])
                    if fn_overlapped_chain:
                        print(c.print_chain(), file=foc)
            
            # If no overlap
            if (not s_overlap) and (not t_overlap):
                # Update source tree dict
                if c.source in stree_dict:
                    stree_dict[c.source] = stree_dict[c.source] | c.stree
                else:
                    stree_dict[c.source] = c.stree
                # Update target tree dict
                if c.target in ttree_dict:
                    ttree_dict[c.target] = ttree_dict[c.target] | c.ttree
                else:
                    ttree_dict[c.target] = c.ttree
                print(c.print_chain(), file=fo)

            c = None


if __name__ == '__main__':
    args = parse_args()
    if args.overlapped_chain:
        if not args.unique:
            print('[E::filter] `-u` must be set if `-oc` is not empty', file=sys.stderr)
            exit(1)

    filter(
        fn_chain=args.chain, fn_out=args.output,
        unique=args.unique, segment_size=args.segment_size,
        fn_overlapped_chain=args.overlapped_chain)
