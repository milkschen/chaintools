# chaintools: utilities for the genomic chain format

This toolkit provides utilities to process whole-genome maps in the ([chain format](https://genome.ucsc.edu/goldenPath/help/chain.html)).

A chain can be used to convert genomic information from the "target" coordinate system to the "query" coordinate system.
For example, in the `hg38ToHg19.over.chain` file hg38 uses the target fields and hg19 uses the query fields.
Lift-over software such as [UCSC LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver), [CrossMap](https://github.com/liguowang/CrossMap), or [levioSAM](https://github.com/alshai/levioSAM) can be used to convert a locus in the hg38 coordinates to hg19's using the chain file.

## Utilities supported

- [Annotate](#annotate)
- [Convert to BED](#convert-to-bed)
- [Convert to PAF](#convert-to-paf)
- [Convert to SAM](#convert-to-sam)
- [Convert to VCF](#convert-to-vcf)
- [Filter](#filter)
- [Invert](#invert)
- [Split](#split)
- [Stats](#stats)

## Install

```shell
git clone git@github.com:milkschen/chaintools.git
```

- Python 3.7+
- Dependencies: [intervaltree](https://github.com/chaimleib/intervaltree), [pandas](https://pandas.pydata.org), and [pysam](https://pysam.readthedocs.io/en/latest/). See [INSTALL.md](INSTALL.md) for instructions.

## Usage

### Set up $PYTHONPATH

```shell
export PYTHONPATH=$(pwd)/:${PYTHONPATH}
```

### Annotate

Annotate a chain file:

- Specify the contig and start/end positions of each segment
- Calculate the identity of each segment (optional)
- Write liftable regions to a pair of BED files (one for target and one for query) (optional)

```shell
# Annotate contig and positions
python chaintools/annotate.py -c <in.chain> -o <out.chain>
# Add identity
python chaintools/annotate.py -c <in.chain> -o <out.chain> -fs <target.fasta> -ft <query.fasta>
# Also write liftable regions to BED files
python chaintools/annotate.py -c <in.chain> -o <out.chain> -fs <target.fasta> -ft <query.fasta> -b <bed_prefix>
```

### Convert to BED

Convert a chain file to the BED format using either target or query coordinates

```shell
# Report using the target coordinates
python chaintools/to_bed.py -c <in.chain> -o <out.bed> --coord target
# Report using the query coordinates
python chaintools/to_bed.py -c <in.chain> -o <out.bed> --coord query
```

### Convert to PAF

Convert a chain file to the [PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md).

The target chain is converted as the target sequence, and the query chain is converted as the query sequence.

If both `target.fa` and `query.fa` are provided, this script checks the reference sequences and updates the cigar (`cg:Z` tag) using `[=XID]+` operators.
Otherwise, it uses `[MID]+` and `[X]+` at chain break points. A breakpoint is a gap wrt both target and query, e.g., `149 341 2894`.

```shell
python chaintools/to_paf.py -c <in.chain> -o <out.paf> [-t <target.fa> -q <query.fa>]
```

### Convert to SAM

Convert a chain file to the [SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf),
using the target fasta file for the genome _from_ which
the chain lifts, and the query fasta file for the genome _to_ which the chain lifts.

```shell
python chaintools/to_sam.py -c <in.chain> -t <target.fa> -q <query.fa> -o <out.sam>
```

Note: For a chain file used to convert from a target genome's coordinates to a query
genome's coordinates, the chain header lines have target data in the second through
sixth fields, and query data in the seventh through eleventh fields.

### Convert to VCF

Convert a chain file to the [VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf),
using the target fasta file for the genome _from_ which
the chain lifts, and the query fasta file for the genome _to_ which the chain lifts.

```shell
python chaintools/to_vcf.py -c <in.chain> -t <target.fa> -q <query.fa> -o <out.vcf>
```

### Filter

Filter a chain file by critera including chain sizes and overlap status.
The size of a chain is the sum of all its segments, including matches (`=`) and mismatches (`X`).
The overlap filter makes sure no chains overlap wrt either target or query references. If two chains overlap, the smaller one is removed.

```shell
# Filter by chain size
python chaintools/chain_filter.py -c <in.chain> -o <out.filtered.chain> -s <size>
# Filter by both chain size and overlap status
python chaintools/chain_filter.py -c <in.chain> -o <out.filtered.chain> -u -oc <out.overlapped.chain> -s <size>
```

### Invert

Invert a chain file by switching the target and query references

```shell
python chaintools/invert.py -c <a_to_b.chain> -o <b_to_a.chain>
```

### Split

Split a chain at large gaps or breakpoints. A breakpoint is a gap wrt both target and query, e.g., `149 341 2894`.

```shell
python chaintools/split.py -c <in.chain> -o <split.chain> [--min_gap <INT> --min_bp <INT>]
```

### Stats

Calculate summary statistics of a chain file

```shell
python chaintools/stats.py -c <in.chain> -o <stats.tsv>
```
