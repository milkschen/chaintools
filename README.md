# chaintools: utilities for the genomic chain format

## Usage

### Annotate
Annotate a chain file:
* Specify the contig and start/end positions of each segment
* Calculate the identity of each segment (optional)
* Write liftable regions to a pair of BED files (one for source and one for target) (optional)

```
python src/annotate.py -c <chain> -o <out>
python src/annotate.py -c <chain> -o <out> -f1 <source.fasta> -f2 <target.fasta>
python src/annotate.py -c <chain> -o <out> -f1 <source.fasta> -f2 <target.fasta> -b <bed_prefix>
```


### Invert
Invert a chain file by switching the source and dest references

```
python src/invert.py -c <chain> -o <out>
```
