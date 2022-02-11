# chaintools: utilities for the genomic chain format

## Usage

### Annotate
Annotate a chain file:
* Specify the contig and start/end positions of each segment
* Calculate the identity of each segment (optional)
* Write liftable regions to a pair of BED files (one for source and one for target) (optional)

```
# Annotate contig and positions
python src/annotate.py -c <in.chain> -o <out.chain>
# Add identity
python src/annotate.py -c <in.chain> -o <out.chain> -fs <source.fasta> -ft <target.fasta>
# Also write liftable regions to BED files
python src/annotate.py -c <in.chain> -o <out.chain> -fs <source.fasta> -ft <target.fasta> -b <bed_prefix>
```


### Invert
Invert a chain file by switching the source and dest references

```
python src/invert.py -c <in.chain> -o <out.paf>
```

### Convert to PAF
Convert a chain file to the ([PAF format](https://github.com/lh3/miniasm/blob/master/PAF.md)). 

The source chain is converted as the query sequence and the target chain is converted as the target sequence.

```
python src/to_paf.py -c <in.chain> -o <out.paf>
```

