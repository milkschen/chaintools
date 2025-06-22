#!/usr/bin/env python3
"""
Split a chain at alignment breakpoints and large gaps

Nae-Chyun Chen
Johns Hopkins University
2022
"""
import argparse
import sys
from typing import TextIO

from chaintools_bio import utils


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--chain", required=True, help="Path to the chain file"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="",
        help="Path to the output BED file. Leave empty to write to stdout.",
    )
    parser.add_argument(
        "--min_gap",
        default=10000,
        type=int,
        help=(
            "Min size of chain gaps to split. Set to a negative value to "
            "diable. Example of a 2000-bp gap `100 2000 0`."
        ),
    )
    parser.add_argument(
        "--min_bp",
        default=1000,
        type=int,
        help=(
            "Min size of chain break point to split. Set to a negative value "
            "to diable. Example of a 341-bp breakpoint: `149 341 2894`."
        ),
    )
    args = parser.parse_args()
    return args


""" Check if a chain line (3-field) should be split

Inputs:
    - dt: gap size wrt target
    - dq: gap size wrt query
    - min_bp: min breakpoint size that we split. An alignment 
              breakpoint refers to a chain segment gap that is
              not simply indelic, such as `149 341 2894`
    - min_gap: min gap size that we split
"""


def check_split(dt: int, dq: int, min_bp: int, min_gap: int) -> bool:
    if min_bp >= 0:
        if min([dt, dq]) > min_bp:
            return True
    if min_gap >= 0:
        if max([dt, dq]) > min_gap:
            return True
    return False


def split_chain(f: TextIO, min_bp: int, min_gap: int) -> str:
    for line in f:
        fields = line.split()
        if len(fields) == 0:
            pass
        elif fields[0] == "chain":
            c = utils.Chain(fields)
        elif len(fields) == 3:
            dt = int(fields[1])
            dq = int(fields[2])
            # Split a chain object if encountering an alignment breakpoint.
            if check_split(dt=dt, dq=dq, min_bp=min_bp, min_gap=min_gap):
                c.add_record([fields[0]])
                tmp_tend = c.tend
                tmp_qend = c.qend
                c.tend = c.toffset
                if c.strand == "+":
                    c.qend = c.qoffset
                else:
                    c.qend = c.qlen - c.qoffset
                yield (c.print_chain())
                # Reset
                c.reset_at_break(dt=dt, dq=dq, tend=tmp_tend, qend=tmp_qend)
            else:
                c.add_record(fields)
        elif len(fields) == 1:
            c.add_record(fields)
            yield (c.print_chain())
            c = None


def split_chain_io(
    fn_chain: str, fn_out: str, min_bp: int, min_gap: int
) -> None:
    if fn_chain == "-":
        f = sys.stdin
    else:
        f = open(fn_chain, "r")
    if fn_out:
        fo = open(fn_out, "w")
    else:
        fo = sys.stdout

    out = split_chain(f=f, min_bp=min_bp, min_gap=min_gap)
    while True:
        try:
            print(next(out), file=fo)
        except StopIteration:
            break


def main(argv=sys.argv):
    args = parse_args()
    print(f"Split parameters:", file=sys.stderr)
    print(f" * min_bp : {args.min_bp}", file=sys.stderr)
    print(f" * min_gap: {args.min_gap}", file=sys.stderr)

    split_chain_io(
        fn_chain=args.chain,
        fn_out=args.output,
        min_bp=args.min_bp,
        min_gap=args.min_gap,
    )


if __name__ == "__main__":
    main()
