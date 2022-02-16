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


def filter_core(
    c: utils.Chain, segment_size: int, unique: bool,
    stree_dict: dict, ttree_dict: dict
) -> list:
    filter_size = False
    filter_overlap_source = False
    filter_overlap_target = False

    # Check segment size
    if c.seglen < segment_size:
        filter_size = True

    # Check overlap
    if unique:
        # If two source trees intersect
        if c.source in stree_dict:
            for s_intvl in c.stree:
                if stree_dict[c.source].overlaps(s_intvl):
                    filter_overlap_source = True
                    break
        # If two target trees intersect
        if c.target in ttree_dict:
            for t_intvl in c.ttree:
                if ttree_dict[c.target].overlaps(t_intvl):
                    filter_overlap_target = True
                    break
    return filter_size, filter_overlap_source, filter_overlap_target


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
            
            filter_size, filter_overlap_source, filter_overlap_target = filter_core(
                c=c, segment_size=segment_size, unique=unique,
                stree_dict=stree_dict, ttree_dict=ttree_dict)

            msg = f'score={c.score}\tsource={c.source}:{c.sstart}-{c.send}\ttarget={c.target}:{c.tstart}-{c.tend}\t{c.strand}\t'

            if filter_size:
                msg += 'SIZE'
            elif filter_overlap_source and fn_overlapped_chain:
                msg += 'OVERLAP_SOURCE'
                print(c.print_chain(), file=foc)
            elif filter_overlap_target and fn_overlapped_chain:
                msg += 'OVERLAP_TARGET'
                print(c.print_chain(), file=foc)
            
            # if pass
            if not any([filter_size, filter_overlap_source, filter_overlap_target]):
                msg += 'PASS'
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

            print(msg, file=sys.stderr)
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
