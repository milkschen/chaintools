'''
Split a chain at alignment breakpoints and large gaps

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
        help='Path to the output BED file. Leave empty to write to stdout.'
    )
    parser.add_argument(
        '--min_gap', default=10000, type=int,
        help='Min size of chain gaps to split. Set to a negative value to diable. Example of a 2000-bp gap `100 2000 0`. [10000]'
    )
    parser.add_argument(
        '--min_bp', default=1000, type=int,
        help='Min size of chain break point to split. Set to a negative value to diable. Example of a 341-bp breakpoint: `149 341 2894`. [1000]'
    )
    args = parser.parse_args()
    return args


def check_split(
    dt: int, dq: int, min_bp: int, min_gap: int
) -> bool:
    if min_bp >= 0:
        if min([dt, dq]) > min_bp:
            return True
    if min_gap >= 0:
        if max([dt, dq]) > min_gap:
            return True
    return False


# TODO: naechyun -- add tests
def split_chain(
    fn_chain: str, fn_out: str, min_bp: int, min_gap: int
) -> None:
    if fn_chain == '-':
        f = sys.stdin
    else:
        f = open(fn_chain, 'r')

    if fn_out:
        fo = open(fn_out, 'w')
    else:
        fo = sys.stdout

    num_orig_chain = 0
    num_split_chain = 0
    for i, line in enumerate(f):
        fields = line.split()
        if len(fields) == 0:
            continue
        elif line.startswith('chain'):
            c = utils.Chain(fields)
            num_orig_chain += 1
        elif len(fields) == 3:
            dt = int(fields[1])
            dq = int(fields[2])
            # Split a chain object if encountering an alignment breakpoint. 
            # An alignment breakpoint refers to a chain segment gap that is not simple indelic, 
            # such as `149 341 2894`
            if check_split(dt=dt, dq=dq, min_bp=min_bp, min_gap=min_gap):
                # print(f'  Split at line {i}:\t{line.rstrip()}', file=sys.stderr)
                c.add_record([fields[0]])
                tmp_tend = c.tend
                tmp_qend = c.qend
                c.tend = c.toffset
                c.qend = c.qoffset
                print(c.print_chain(), file=fo)
                num_split_chain += 1
                # Reset
                c.reset_at_break(
                    dt=dt, dq=dq, tend=tmp_tend, qend=tmp_qend)
            else:
                c.add_record(fields)
        elif len(fields) == 1:
            c.add_record(fields)
            print(c.print_chain(), file=fo)
            num_split_chain += 1

    print(f'Num original chains: {num_orig_chain}', file=sys.stderr)
    print(f'Num split chains   : {num_split_chain}', file=sys.stderr)


if __name__ == '__main__':
    args = parse_args()
    print(f'Split parameters:')
    print(f' * min_bp : {args.min_bp}', file=sys.stderr)
    print(f' * min_gap: {args.min_gap}', file=sys.stderr)

    split_chain(
        fn_chain=args.chain, fn_out=args.output,
        min_bp=args.min_bp, min_gap=args.min_gap)
