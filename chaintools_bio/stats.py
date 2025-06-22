#!/usr/bin/env python3
"""
Filtering a chain file

Nae-Chyun Chen
Johns Hopkins University
2022
"""
import argparse
import sys

import pandas as pd

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
        help="Path to the output merged chain file.",
    )
    args = parser.parse_args()
    return args


def stats(fn_chain: str, fn_out: str):
    f = open(fn_chain, "r")

    print(f"TARGET\tQUERY\tSCORE\tSTRAND\tSEGLEN", file=sys.stderr)
    num_forward = 0
    num_reversed = 0
    seglen_forward = 0
    seglen_reversed = 0

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

            if c.strand == "+":
                num_forward += 1
                seglen_forward += c.seglen
            else:
                num_reversed += 1
                seglen_reversed += c.seglen

            msg = (
                f"{c.target}:{c.tstart}-{c.tend}\t"
                f"{c.query}:{c.qstart}-{c.qend}"
                f"\t{c.score}\t{c.strand}\t{c.seglen}"
            )
            print(msg, file=sys.stderr)

            c = None
    data = [
        [num_forward, seglen_forward],
        [num_reversed, seglen_reversed],
        [num_forward + num_reversed, seglen_forward + seglen_reversed],
    ]
    idx = ["Forward", "Reversed", "Total"]
    hdr = ["NumChains", "SegmentSize"]
    df = pd.DataFrame(data, idx, hdr)
    if fn_out:
        df.to_csv(fn_out, sep="\t")
    else:
        print(df)


def main(argv=sys.argv):
    args = parse_args()
    stats(fn_chain=args.chain, fn_out=args.output)


if __name__ == "__main__":
    main()
