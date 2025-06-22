#!/usr/bin/env python3
"""
Annotate a chain file

Nae-Chyun Chen
Johns Hopkins University
2021-2022
"""
import argparse
import re
import sys

from chaintools_bio import utils


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-c", "--chain", required=True, help="Path to the input chain file."
    )
    parser.add_argument(
        "-o",
        "--out",
        default="",
        help=("Path to the output verbose chain " "file. [empty string]"),
    )
    parser.add_argument(
        "-b",
        "--bed_prefix",
        default="",
        help=(
            "Prefix to the BED files that specify "
            "the aligned regions in `-c`. [empty string]"
        ),
    )
    parser.add_argument(
        "-s",
        "--summary",
        default="",
        help=(
            "Path to the output summary. "
            "Leave empty for no output. [empty string]"
        ),
    )
    parser.add_argument(
        "-fs",
        "--s_ref",
        default="",
        help="Path to the source reference (optional).",
    )
    parser.add_argument(
        "-ft",
        "--t_ref",
        default="",
        help="Path to the destination reference (optional).",
    )
    args = parser.parse_args()
    return args


""" Write a chain record in the BED format """


def write_to_summary(
    fs_fn, fs, strand, l, hd, source, s_start, target, t_start
):
    if fs_fn:
        if strand == "+":
            print(
                (
                    f"{l}\t{hd:.6f}\t{source}\t{s_start}\t{s_start+l}\t+"
                    f"\t{target}\t{t_start}\t{t_start+l}"
                ),
                file=fs,
            )
        else:
            print(
                (
                    f"{l}\t{hd:.6f}\t{source}\t{s_start}\t{s_start+l}\t+"
                    f"\t{target}\t{t_start-l}\t{t_start}"
                ),
                file=fs,
            )


def annotate(
    chain: str,
    out: str,
    bed_prefix: str,
    summary: str,
    s_ref: dict,
    t_ref: dict,
) -> None:
    f = open(chain, "r")
    if out == "":
        fo = sys.stderr
    else:
        fo = open(out, "w")
    if bed_prefix != "":
        f_sbed = open(bed_prefix + "-source.bed", "w")
        f_dbed = open(bed_prefix + "-target.bed", "w")

    check_hdist = True if (s_ref != {} and t_ref != {}) else False
    if not check_hdist:
        hd = None
    if summary:
        assert check_hdist == True
        fs = open(summary, "w")
        # Write header
        print(
            "SIZE\tHDIST\tSOURCE\tS_START\tS_END\t"
            "STRAND\tTARGET\tT_START\tT_END",
            file=fs,
        )
    else:
        fs = None

    total_bases = 0
    total_bases_idy = 0
    for line in f:
        if not line:
            continue
        line = line.rstrip()
        fields = re.split(r"[\s\t]+", line)
        if fields[0] == "chain":
            print(line, file=fo)
            source = fields[2]
            s_start = int(fields[5])
            s_end = int(fields[6])
            target = fields[7]
            target_len = int(fields[8])
            strand = fields[9]
            if strand == "+":
                t_start = int(fields[10])
                # t_end = int(fields[11])
            else:
                t_start = target_len - int(fields[10])
                # t_end = target_len - int(fields[11])
        elif len(fields) == 3:
            l = int(fields[0])
            total_bases += l
            ds = int(fields[1])
            dd = int(fields[2])
            if strand == "+":
                msg = (
                    f"\t{source}:{s_start}-{s_start+l}=>{target}:"
                    f"{t_start}-{t_start+l} ({t_start-s_start})"
                )
                if bed_prefix != "":
                    f_sbed.write(f"{source}\t{s_start}\t{s_start+l}\n")
                    f_dbed.write(f"{target}\t{t_start}\t{t_start+l}\n")
                if check_hdist:
                    hd = utils.compute_hamming_dist(
                        True,
                        s_ref,
                        source,
                        s_start,
                        s_start + l,
                        t_ref,
                        target,
                        t_start,
                        t_start + l,
                    )
                    msg += f"\t{hd:.2f}"
            else:
                msg = (
                    f"\t{source}:{s_start}-{s_start+l}=>{target}:"
                    f"{t_start}-{t_start-l} ({t_start-s_start})"
                )
                if bed_prefix != "":
                    f_sbed.write(f"{source}\t{s_start}\t{s_start+l}\n")
                    f_dbed.write(f"{target}\t{t_start-l}\t{t_start}\n")
                if check_hdist:
                    hd = utils.compute_hamming_dist(
                        False,
                        s_ref,
                        source,
                        s_start,
                        s_start + l,
                        t_ref,
                        target,
                        t_start - l,
                        t_start,
                    )
                    msg += f"\t{hd:.2f}"
            print(line + msg, file=fo)
            write_to_summary(
                summary, fs, strand, l, hd, source, s_start, target, t_start
            )
            if strand == "+":
                s_start += l + ds
                t_start += l + dd
            else:
                s_start += l + ds
                t_start -= l + dd
            if check_hdist:
                total_bases_idy += l * hd
        elif len(fields) == 1 and fields[0] != "":
            l = int(fields[0])
            total_bases += l
            if strand == "+":
                msg = f"\t\t\t{source}:{s_start}-{s_start+l}=>{target}:{t_start}-{t_start+l}"
                if bed_prefix != "":
                    f_sbed.write(f"{source}\t{s_start}\t{s_start+l}\n")
                    f_dbed.write(f"{target}\t{t_start}\t{t_start+l}\n")
                if check_hdist:
                    hd = utils.compute_hamming_dist(
                        True,
                        s_ref,
                        source,
                        s_start,
                        s_start + l,
                        t_ref,
                        target,
                        t_start,
                        t_start + l,
                    )
                    msg += f"\t{hd:.2f}"
            else:
                msg = f"\t\t\t{source}:{s_start}-{s_start+l}=>{target}:{t_start}-{t_start-l}"
                if bed_prefix != "":
                    f_sbed.write(f"{source}\t{s_start}\t{s_start+l}\n")
                    f_dbed.write(f"{target}\t{t_start-l}\t{t_start}\n")
                if check_hdist:
                    hd = utils.compute_hamming_dist(
                        False,
                        s_ref,
                        source,
                        s_start,
                        s_start + l,
                        t_ref,
                        target,
                        t_start - l,
                        t_start,
                    )
                    msg += f"\t{hd:.2f}"
            write_to_summary(
                summary, fs, strand, l, hd, source, s_start, target, t_start
            )
            if check_hdist:
                total_bases_idy += l * hd
            print(line + msg + "\n", file=fo)

    print(
        f"Total number of gapless aligned bases = {total_bases}",
        file=sys.stderr,
    )
    print(
        f"Total number of gapless matched bases = {total_bases_idy}",
        file=sys.stderr,
    )


def main(argv=sys.argv):
    args = parse_args()

    print("Input chain             :", args.chain, file=sys.stderr)
    print("Output chain (annotated):", args.out, file=sys.stderr)

    assert (args.s_ref != "" and args.t_ref != "") or (
        args.s_ref == "" and args.t_ref == ""
    )
    print("Source reference        :", args.s_ref, file=sys.stderr)
    print("Target reference        :", args.t_ref, file=sys.stderr)

    s_ref = utils.read_fasta(args.s_ref)
    t_ref = utils.read_fasta(args.t_ref)

    annotate(
        chain=args.chain,
        out=args.out,
        bed_prefix=args.bed_prefix,
        summary=args.summary,
        s_ref=s_ref,
        t_ref=t_ref,
    )


if __name__ == "__main__":
    main()
