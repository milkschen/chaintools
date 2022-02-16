import argparse
import utils
import sys

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-c', '--chain', required=True,
        help='Path to the chain file'
    )
    parser.add_argument(
        '-s', '--sourcefasta', required=True,
        help='Path to the fasta file for the source genome'
    )
    parser.add_argument(
        '-t', '--targetfasta', required=True,
        help='Path to the fasta file for the target genome'
    )
    parser.add_argument(
        '-o', '--output', default='',
        help='Path to the output VCF file.'
    )
    args = parser.parse_args()
    return args

def write_to_vcf(fn_chain: str, fn_vcf: str, fn_sourcefasta: str, fn_targetfasta: str):
    f = open(fn_chain, 'r')
    if fn_vcf:
        fo = open(fn_vcf, 'w')
    else:
        fo = sys.stdout

    print(utils.vcf_header(utils.get_source_entries(fn_chain)), file=fo, end='')

    sourceref = utils.fasta_reader(fn_sourcefasta)
    targetref = utils.fasta_reader(fn_targetfasta)

    for line in f:
        fields = line.split()
        if len(fields) == 0:
            continue
        elif line.startswith('chain'):
            c = utils.Chain(fields)
        elif len(fields) == 3:
            c.add_record_three(fields)
        elif len(fields) == 1:
            c.add_record_one(fields)
            print(c.to_vcf(sourceref, targetref), file=fo, end='')
            c = None

if __name__ == '__main__':
    args = parse_args()
    write_to_vcf(fn_chain=args.chain, fn_vcf=args.output, fn_sourcefasta=args.sourcefasta, fn_targetfasta=args.targetfasta)
