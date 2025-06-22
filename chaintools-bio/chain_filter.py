#!/usr/bin/env python3
"""
Filter a chain file

Nae-Chyun Chen
Johns Hopkins University
2022
"""
import argparse
import sys

from chaintools_bio import utils


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--chain", default="", help="Path to the chain file"
    )
    parser.add_argument(
        "-o",
        "--output",
        default="",
        help="Path to the output merged chain file.",
    )
    parser.add_argument(
        "-u",
        "--unique",
        action="store_true",
        help="Activate to remove mappings that are not 1-1. Chains with smaller scores are excluded.",
    )
    parser.add_argument(
        "-oc",
        "--overlapped_chain",
        default="",
        help="Path to the chains that overlap with higher-scoring chains in `-c`. Leave this field empty to not write. Must with `-u`. [None]",
    )
    parser.add_argument(
        "-s",
        "--segment_size",
        default=0,
        type=int,
        help="Minimal segment size (the sum of chain segment sizes) allowed. [0]",
    )
    args = parser.parse_args()
    return args


def filter_core(
    c: utils.Chain,
    segment_size: int,
    unique: bool,
    ttree_dict: dict,
    qtree_dict: dict,
) -> list:
    filter_size = False
    filter_overlap_target = False
    filter_overlap_query = False

    # Check segment size
    if c.seglen < segment_size:
        filter_size = True

    # Check overlap
    if unique:
        # If two target trees intersect
        if c.target in ttree_dict:
            for t_intvl in c.ttree:
                if ttree_dict[c.target].overlaps(t_intvl):
                    filter_overlap_target = True
                    break
        # If two query trees intersect
        if c.query in qtree_dict:
            for q_intvl in c.qtree:
                if qtree_dict[c.query].overlaps(q_intvl):
                    filter_overlap_query = True
                    break
    return filter_size, filter_overlap_target, filter_overlap_query


def chain_filter(
    fn_chain: str,
    fn_out: str,
    unique: bool,
    segment_size: int,
    fn_overlapped_chain: str = "",
):
    stree_dict = {}
    qtree_dict = {}

    if fn_chain == "-":
        f = sys.stdin
    else:
        f = open(fn_chain, "r")
    if fn_out:
        fo = open(fn_out, "w")
    else:
        fo = sys.stdout
    if fn_overlapped_chain:
        foc = open(fn_overlapped_chain, "w")

    for line in f:
        fields = line.split()
        if len(fields) == 0:
            continue
        elif line.startswith("chain"):
            c = utils.Chain(fields)
        elif len(fields) == 3:
            c.add_record(fields)
        elif len(fields) == 1:
            c.add_record(fields)

            (
                filter_size,
                filter_overlap_target,
                filter_overlap_query,
            ) = filter_core(
                c=c,
                segment_size=segment_size,
                unique=unique,
                ttree_dict=stree_dict,
                qtree_dict=qtree_dict,
            )

            msg = (
                f"score={c.score}\ttarget={c.target}:{c.tstart}-{c.tend}\t"
                f"query={c.query}:{c.qstart}-{c.qend}\t{c.strand}\t"
            )

            if filter_size:
                msg += "SIZE"
            elif filter_overlap_target and fn_overlapped_chain:
                msg += "OVERLAP_TARGET"
                print(c.print_chain(), file=foc)
            elif filter_overlap_query and fn_overlapped_chain:
                msg += "OVERLAP_QUERY"
                print(c.print_chain(), file=foc)

            # if pass
            if not any(
                [filter_size, filter_overlap_target, filter_overlap_query]
            ):
                msg += "PASS"
                # Update target tree dict
                if c.target in stree_dict:
                    stree_dict[c.target] = stree_dict[c.target] | c.ttree
                else:
                    stree_dict[c.target] = c.ttree
                # Update query tree dict
                if c.query in qtree_dict:
                    qtree_dict[c.query] = qtree_dict[c.query] | c.qtree
                else:
                    qtree_dict[c.query] = c.qtree
                print(c.print_chain(), file=fo)

            print(msg, file=sys.stderr)
            c = None


def main(argv=sys.argv):
    args = parse_args()
    if args.overlapped_chain:
        if not args.unique:
            print(
                "[E::chain_filter] `-u` must be set if `-oc` is not empty",
                file=sys.stderr,
            )
            exit(1)

    chain_filter(
        fn_chain=args.chain,
        fn_out=args.output,
        unique=args.unique,
        segment_size=args.segment_size,
        fn_overlapped_chain=args.overlapped_chain,
    )


if __name__ == "__main__":
    main()
