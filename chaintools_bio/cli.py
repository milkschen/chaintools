#!/usr/bin/env python3
"""
Main CLI interface for chaintools_bio using Typer
"""
import typer
from enum import Enum

from chaintools_bio import (
    __version__,
    to_bed,
    to_paf,
    to_sam,
    to_vcf,
    chain_filter,
    annotate,
    invert,
    split,
    stats,
)

app = typer.Typer(
    help="Utilities for the genomic chain format",
    no_args_is_help=True,
)


class CoordSystem(str, Enum):
    target = "target"
    query = "query"


@app.command("to-bed")
def to_bed_cmd(
    chain: str = typer.Option(
        "-", "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output BED file"
    ),
    coord: CoordSystem = typer.Option(
        CoordSystem.target, "--coord", help="Output coordinate system"
    ),
):
    """Convert a chain file to the BED format"""
    to_bed.write_to_bed_io(fn_chain=chain, fn_bed=output, coord=coord.value)


@app.command("to-paf")
def to_paf_cmd(
    chain: str = typer.Option(
        "-", "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output PAF file"
    ),
    target_fasta: str = typer.Option(
        "", "-t", "--target-fasta", help="Path to the target FASTA file"
    ),
    query_fasta: str = typer.Option(
        "", "-q", "--query-fasta", help="Path to the query FASTA file"
    ),
    preserve_breakpoint: bool = typer.Option(
        False,
        "--preserve-breakpoint",
        help="Preserve breakpoint information in PAF output",
    ),
):
    """Convert a chain file to the PAF format"""
    to_paf.write_to_paf_io(
        fn_chain=chain,
        fn_paf=output,
        fn_targetfasta=target_fasta,
        fn_queryfasta=query_fasta,
        preserve_breakpoint=preserve_breakpoint,
    )


@app.command("to-sam")
def to_sam_cmd(
    chain: str = typer.Option(
        "-", "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output SAM file"
    ),
    target_fasta: str = typer.Option(
        ..., "-t", "--target-fasta", help="Path to the target FASTA file"
    ),
    query_fasta: str = typer.Option(
        ..., "-q", "--query-fasta", help="Path to the query FASTA file"
    ),
):
    """Convert a chain file to the SAM format"""
    to_sam.write_to_sam_io(
        fn_chain=chain,
        fn_paf=output,
        fn_targetfasta=target_fasta,
        fn_queryfasta=query_fasta,
    )


@app.command("to-vcf")
def to_vcf_cmd(
    chain: str = typer.Option(
        "-", "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output VCF file"
    ),
    target_fasta: str = typer.Option(
        ..., "-t", "--target-fasta", help="Path to the target FASTA file"
    ),
    query_fasta: str = typer.Option(
        ..., "-q", "--query-fasta", help="Path to the query FASTA file"
    ),
):
    """Convert a chain file to the VCF format"""
    to_vcf.write_to_vcf_io(
        fn_chain=chain,
        fn_vcf=output,
        fn_targetfasta=target_fasta,
        fn_queryfasta=query_fasta,
    )


@app.command("filter")
def filter_cmd(
    chain: str = typer.Option(
        "", "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output chain file"
    ),
    unique: bool = typer.Option(
        False, "-u", "--unique", help="Remove mappings that are not 1-1"
    ),
    segment_size: int = typer.Option(
        0, "-s", "--segment-size", help="Minimal segment size allowed"
    ),
    overlapped_chain: str = typer.Option(
        "", "--overlapped-chain", help="Path to output overlapped chains file"
    ),
):
    """Filter chain files based on various criteria"""
    chain_filter.chain_filter(
        fn_chain=chain,
        fn_out=output,
        unique=unique,
        segment_size=segment_size,
        fn_overlapped_chain=overlapped_chain,
    )


@app.command("annotate")
def annotate_cmd(
    chain: str = typer.Option(
        ..., "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output annotated chain file"
    ),
    bed_prefix: str = typer.Option(
        "", "-b", "--bed-prefix", help="Prefix for BED output files"
    ),
    summary: str = typer.Option(
        "", "-s", "--summary", help="Path to summary output file"
    ),
    source_fasta: str = typer.Option(
        "", "--source-fasta", help="Path to source FASTA file"
    ),
    target_fasta: str = typer.Option(
        "", "--target-fasta", help="Path to target FASTA file"
    ),
):
    """Annotate chain files with BED file information"""
    import pysam

    s_ref = pysam.FastaFile(source_fasta) if source_fasta else {}
    t_ref = pysam.FastaFile(target_fasta) if target_fasta else {}
    annotate.annotate(
        chain=chain,
        out=output,
        bed_prefix=bed_prefix,
        summary=summary,
        s_ref=s_ref,
        t_ref=t_ref,
    )


@app.command("invert")
def invert_cmd(
    chain: str = typer.Option(
        "-", "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output inverted chain file"
    ),
):
    """Invert a chain file (swap target and query)"""
    invert.invert(in_fn=chain, out_fn=output)


@app.command("split")
def split_cmd(
    chain: str = typer.Option(
        ..., "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output split chain file"
    ),
    min_bp: int = typer.Option(
        1000, "--min-bp", help="Minimum breakpoint size to split"
    ),
    min_gap: int = typer.Option(
        10000, "--min-gap", help="Minimum gap size to split"
    ),
):
    """Split a chain file by breakpoints"""
    split.split_chain_io(
        fn_chain=chain, fn_out=output, min_bp=min_bp, min_gap=min_gap
    )


@app.command("stats")
def stats_cmd(
    chain: str = typer.Option(
        ..., "-c", "--chain", help="Path to the chain file"
    ),
    output: str = typer.Option(
        "", "-o", "--output", help="Path to the output statistics file"
    ),
):
    """Generate statistics for a chain file"""
    stats.stats(fn_chain=chain, fn_out=output)


def _version_callback(value: bool):
    if value:
        print(f"chaintools_bio version: {__version__}")
        raise typer.Exit()


@app.callback(invoke_without_command=True)
def version_callback(
    version: bool = typer.Option(
        None, "--version", callback=_version_callback, is_eager=True
    ),
):
    pass


def main():
    """Primary entry point for the CLI."""
    version_callback()
    app()


if __name__ == "__main__":
    main()
