# chaintools: utilities for the genomic chain format

This toolkit provides utilities to process whole-genome maps in the ([chain format](https://genome.ucsc.edu/goldenPath/help/chain.html)).

A chain can be used to convert genomic information from the "target" coordinate system to the "query" coordinate system.
For example, in the `hg38ToHg19.over.chain` file hg38 uses the target fields and hg19 uses the query fields.
Lift-over software such as [UCSC LiftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver), [CrossMap](https://github.com/liguowang/CrossMap), or [levioSAM](https://github.com/alshai/levioSAM) can be used to  convert a locus in the hg38 coordinates to hg19's using the chain file.


## Utilities supported
* [Annotate](#annotate)
* [Convert to BED](#to_bed)
* [Convert to PAF](#to_paf)
* [Convert to SAM](#to_sam)
* [Convert to VCF](#to_vcf)
* [Filter](#filter)
* [Invert](#invert)
* [Split](#split)
* [Stats](#stats)


## Install

```
git clone git@github.com:milkschen/chaintools.git
```

* Dependencies: [intervaltree](https://github.com/chaimleib/intervaltree), [pandas](https://pandas.pydata.org), and [pysam](https://pysam.readthedocs.io/en/latest/). See [INSTALL.md](INSTALL.md) for instructions.



## Usage

<a name="annotate"></a>
### Annotate 
Annotate a chain file:
* Specify the contig and start/end positions of each segment
* Calculate the identity of each segment (optional)
* Write liftable regions to a pair of BED files (one for target and one for query) (optional)

```
# Annotate contig and positions
python src/annotate.py -c <in.chain> -o <out.chain>
# Add identity
python src/annotate.py -c <in.chain> -o <out.chain> -fs <target.fasta> -ft <query.fasta>
# Also write liftable regions to BED files
python src/annotate.py -c <in.chain> -o <out.chain> -fs <target.fasta> -ft <query.fasta> -b <bed_prefix>
```


<a name="to_bed"></a>
### Convert to BED
Convert a chain file to the BED format using either target or query coordinates

```
# Report using the target coordinates
python src/to_bed.py -c <in.chain> -o <out.bed> --coord target
# Report using the query coordinates
python src/to_bed.py -c <in.chain> -o <out.bed> --coord query
```


<a name="to_paf"></a>
### Convert to PAF
Convert a chain file to the ([PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md)). 

The target chain is converted as the target sequence, and the query chain is converted as the query sequence.

If both `target.fa` and `query.fa` are provided, this script checks the reference sequences and updates the cigar (`cg:Z` tag) using `[=XID]+` operators.
Otherwise, it uses `[MID]+` and `[X]+` at chain break points. A breakpoint is a gap wrt both target and query, e.g., `149 341 2894`.

```
python src/to_paf.py -c <in.chain> -o <out.paf> [-t <target.fa> -q <query.fa>]
```


<a name="to_sam"></a>
### Convert to SAM
Convert a chain file to the ([SAM format](https://samtools.github.io/hts-specs/SAMv1.pdf)), 
using the target fasta file for the genome *from* which
the chain lifts, and the query fasta file for the genome *to* which the chain lifts.

```
python src/to_sam.py -c <in.chain> -t <target.fa> -q <query.fa> -o <out.sam> 
```

Note: For a chain file used to convert from a target genome's coordinates to a query
genome's coordinates, the chain header lines have target data in the second through
sixth fields, and query data in the seventh through eleventh fields.


<a name="to_vcf"></a>
### Convert to VCF
Convert a chain file to the ([VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf)),
using the target fasta file for the genome *from* which
the chain lifts, and the query fasta file for the genome *to* which the chain lifts.

```
python src/to_vcf.py -c <in.chain> -t <target.fa> -q <query.fa> -o <out.vcf>
```


<a name="filter"></a>
### Filter
Filter a chain file by critera including chain sizes and overlap status. 
The size of a chain is the sum of all its segments, including matches (`=`) and mismatches (`X`). 
The overlap filter makes sure no chains overlap wrt either target or query references. If two chains overlap, the smaller one is removed.

```
# Filter by chain size
python src/filter.py -c <in.chain> -o <out.filtered.chain> -s <size>
# Filter by both chain size and overlap status
python src/filter.py -c <in.chain> -o <out.filtered.chain> -u -oc <out.overlapped.chain> -s <size>
```


<a name="invert"></a>
### Invert
Invert a chain file by switching the target and query references

```
python src/invert.py -c <a_to_b.chain> -o <b_to_a.chain>
```


<a name="split"></a>
### Split 
Split a chain at large gaps or breakpoints. A breakpoint is a gap wrt both target and query, e.g., `149 341 2894`.

```
python src/split.py -c <in.chain> -o <split.chain> [--min_gap <INT> --min_bp <INT>]
```


<a name="stats"></a>
### Stats
Calculate summary statistics of a chain file

```
python src/stats.py -c <in.chain> -o <stats.tsv>
```
