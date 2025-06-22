#!/usr/bin/env python3
"""
Convert a chain file to the SAM format

Nancy Fisher Hansen
NIH/NHGRI

Nae-Chyun Chen
Johns Hopkins University

2021-2022
"""
import argparse
import sys
from typing import TextIO

import pysam

from chaintools_bio import utils


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--chain", default="", help="Path to the chain file"
    )
    parser.add_argument(
        "-t",
        "--targetfasta",
        required=True,
        help=(
            "Path to the fasta file for the "
            "target (reference) genome of the chain file"
        ),
    )
    parser.add_argument(
        "-q",
        "--queryfasta",
        required=True,
        help=(
            "Path to the fasta file for the query " "genome of the chain file"
        ),
    )
    parser.add_argument(
        "-o", "--output", default="", help="Path to the output VCF file."
    )
    args = parser.parse_args()
    return args


def write_to_vcf(
    f: TextIO,
    targetref: pysam.FastaFile = None,
    queryref: pysam.FastaFile = None,
) -> str:
    for line in f:
        fields = line.split()
        if len(fields) == 0:
            pass
        elif fields[0] == "chain":
            c = utils.Chain(fields)
        else:
            c.add_record(fields)
            if len(fields) == 1:
                yield c.to_vcf(targetref=targetref, queryref=queryref)
                c = None


def write_to_vcf_io(
    fn_chain: str,
    fn_vcf: str,
    fn_targetfasta: str = "",
    fn_queryfasta: str = "",
) -> None:
    if fn_chain == "-":
        f = sys.stdin
    else:
        f = open(fn_chain, "r")
    if fn_vcf:
        fo = open(fn_vcf, "w")
    else:
        fo = sys.stdout

    targetref = utils.fasta_reader(fn_targetfasta)
    queryref = utils.fasta_reader(fn_queryfasta)

    print(utils.vcf_header(utils.get_target_entries(fn_chain)), file=fo)

    out = write_to_vcf(f=f, targetref=targetref, queryref=queryref)
    while True:
        try:
            nextline = next(out)
            if nextline:
                print(nextline, file=fo)
        except StopIteration:
            break


def main(argv=sys.argv):
    args = parse_args()
    write_to_vcf_io(
        fn_chain=args.chain,
        fn_vcf=args.output,
        fn_targetfasta=args.targetfasta,
        fn_queryfasta=args.queryfasta,
    )


if __name__ == "__main__":
    main()
