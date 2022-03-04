'''
Tests for chaintools utilities

Nae-Chyun Chen
Johns Hopkins University

Nancy Fisher Hansen
NIH/NHGRI

2021-2022
'''
import intervaltree
import unittest
import sys
# chaintools
import filter
import to_paf
import utils


class TestReadingChain(unittest.TestCase):
    def read_and_compare(self, fn):
        input_txt = ''
        output_txt = ''
        with open(fn, 'r') as f:
            for line in f:
                input_txt += line
                fields = line.split()
                if len(fields) == 0:
                    continue
                elif line.startswith('chain'):
                    c = utils.Chain(fields)
                elif len(fields) == 3:
                    c.add_record_three(fields)
                elif len(fields) == 1:
                    c.add_record_one(fields)
                    output_txt += c.print_chain()
                    output_txt += '\n'
                    c = None
        
        # A tab and a space are considered to be indifferent
        input_txt = input_txt.replace('\t', ' ').rstrip()
        output_txt = output_txt.replace('\t', ' ').rstrip()
        ii = input_txt.split('\n')
        oo = output_txt.split('\n')
        if len(ii) != len(oo):
            print(f'[E::read_and_compare] Lengths differ ({len(ii)} vs. {len(oo)})', file=sys.stderr)
            return False
        for i in range(min(len(ii), len(oo))):
            if ii[i] != oo[i]:
                print(f'[E::read_and_compare] {ii[i]} != {oo[i]}', file=sys.stderr)
                return False
        return input_txt == output_txt

    def test_read_forward(self):
        fn = 'testdata/forward.chain'
        self.assertTrue(self.read_and_compare(fn),
                        f'Failed when reading {fn}')

    def test_read_reversed(self):
        fn = 'testdata/reversed.chain'
        self.assertTrue(self.read_and_compare(fn),
                        f'Failed when reading {fn}')

    def test_read_forward_end_with_zero(self):
        fn = 'testdata/end_with_zero.chain'
        self.assertTrue(self.read_and_compare(fn),
                        f'Failed when reading {fn}')
        
    def test_read_forward_small(self):
        fn = 'testdata/small.chain'
        self.assertTrue(self.read_and_compare(fn),
                        f'Failed when reading {fn}')


class TestGenerateVcf(unittest.TestCase):
    def generate_and_check(self, chainfn, targetfn, queryfn, vcffn):
        targetref = utils.fasta_reader(targetfn)
        queryref = utils.fasta_reader(queryfn)
        output_txt = utils.vcf_header(utils.get_target_entries(chainfn))
        with open(chainfn, 'r') as f:
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
                    output_txt += c.to_vcf(targetref, queryref)
                    c = None
        with open(vcffn, 'r') as f:
            vcf_txt = ''
            output_lines = output_txt.split('\n')
            for i, line in enumerate(f):
                if line.rstrip() != output_lines[i]:
                    print('\n   truth vcf line:', line.rstrip())
                    print('reported vcf line:', output_lines[i])
                    return False
                vcf_txt += line

        return output_txt == vcf_txt

    def test_generate_vcf_from_small(self):
        fn = 'testdata/target-query.chain'
        targetfn = 'testdata/target.fasta'
        queryfn = 'testdata/query.fasta'
        vcffn = 'testdata/target-query.vcf'
        
        self.assertTrue(self.generate_and_check(fn, targetfn, queryfn, vcffn),
                        f'Failed when generating vcf from {fn}')


class TestGenerateSAM(unittest.TestCase):
    def generate_and_check(self, chainfn, targetfn, queryfn, samfn):
        output_txt = ''
        targetref = utils.fasta_reader(targetfn)
        queryref = utils.fasta_reader(queryfn)
        output_txt += utils.sam_header(utils.get_target_entries(chainfn))
        with open(chainfn, 'r') as f:
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
                    output_txt += c.to_sam(targetref, queryref)
                    c = None
        with open(samfn, 'r') as f:
            sam_txt = ''
            for line in f:
                sam_txt += line

        return output_txt == sam_txt

    def test_generate_sam_from_small(self):
        fn = 'testdata/target-query.chain'
        targetfn = 'testdata/target.fasta'
        queryfn = 'testdata/query.fasta'
        samfn = 'testdata/target-query.sam'
        
        self.assertTrue(self.generate_and_check(fn, targetfn, queryfn, samfn),
                        f'Failed when generating sam from {fn}')


class TestGeneratePAF(unittest.TestCase):
    def generate_and_check(self, chainfn, targetfn, queryfn, samfn):
        f = open(chainfn, 'r')
        targetref = utils.fasta_reader(targetfn)
        queryref = utils.fasta_reader(queryfn)
        output_txt = ''
        out = to_paf.write_to_paf(f, targetref, queryref)
        while True:
            try:
                output_txt += (next(out) + '\n')
            except StopIteration:
                break
        f.close()
        with open(samfn, 'r') as f:
            paf_txt = ''
            for line in f:
                paf_txt += line

        return output_txt == paf_txt

    def test_generate_paf_from_small(self):
        fn = 'testdata/target-query.chain'
        targetfn = 'testdata/target.fasta'
        queryfn = 'testdata/query.fasta'
        paffn = 'testdata/target-query.paf'
        
        self.assertTrue(self.generate_and_check(fn, targetfn, queryfn, paffn),
                        f'Failed when generating PAF from {fn}')
    
    def test_generate_paf_from_small_no_ref(self):
        fn = 'testdata/target-query.chain'
        paffn = 'testdata/target-query.no_ref.paf'
        
        self.assertTrue(self.generate_and_check(fn, '', '', paffn),
                        f'Failed when generating PAF from {fn}')

    def test_generate_paf_from_forward_no_ref(self):
        fn = 'testdata/forward.chain'
        paffn = 'testdata/forward.no_ref.paf'
        
        self.assertTrue(self.generate_and_check(fn, '', '', paffn),
                        f'Failed when generating PAF from {fn}')

    def test_generate_paf_from_reversed_no_ref(self):
        fn = 'testdata/reversed.chain'
        paffn = 'testdata/reversed.no_ref.paf'
        
        self.assertTrue(self.generate_and_check(fn, '', '', paffn),
                        f'Failed when generating PAF from {fn}')


class TestGenerateBED(unittest.TestCase):
    def generate_and_check(self, chainfn, bedfn, coord):
        output_txt = ''
        with open(chainfn, 'r') as f:
            for line in f:
                fields = line.split()
                if len(fields) == 0:
                    continue
                elif line.startswith('chain'):
                    c = utils.Chain(fields)
                else:
                    bed_str = c.record_to_bed(fields=fields, coord=coord)
                    if bed_str != '':
                        output_txt += (bed_str + '\n')
                    c.add_record(fields)
        with open(bedfn, 'r') as f:
            bed_txt = ''
            for line in f:
                bed_txt += line
        return output_txt == bed_txt

    def test_generate_bed_from_small_target(self):
        fn = 'testdata/target-query.chain'
        bedfn = 'testdata/target-query.target.bed'
        self.assertTrue(self.generate_and_check(fn, bedfn, 'target'),
                        f'Failed when generating BED (target) from {fn}')
    
    def test_generate_bed_from_small_query(self):
        fn = 'testdata/target-query.chain'
        bedfn = 'testdata/target-query.query.bed'
        self.assertTrue(self.generate_and_check(fn, bedfn, 'query'),
                        f'Failed when generating BED (query) from {fn}')


class TestFilter(unittest.TestCase):
    def read_filter_compare(
        self, fn: str, segment_size: int, unique: bool,
        ttree_dict: dict, qtree_dict: dict
    ) -> list:
        filter_results = []
        with open(fn, 'r') as f:
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
                    filter_results.append(filter.filter_core(
                        c=c, segment_size=segment_size, unique=unique,
                        ttree_dict=ttree_dict, qtree_dict=qtree_dict))
                    c = None
        return filter_results

    def test_read_forward_no_filter(self):
        fn = 'testdata/forward.chain'
        filter_results = self.read_filter_compare(
            fn=fn, segment_size=0, unique=False,
            ttree_dict={}, qtree_dict={})
        self.assertTrue(not any(filter_results[0]), f'Failed when reading {fn} [0]')
        self.assertTrue(not any(filter_results[1]), f'Failed when reading {fn} [1]')
    
    def test_read_forward_filter_size(self):
        fn = 'testdata/forward.chain'
        filter_results = self.read_filter_compare(
            fn=fn, segment_size=360000, unique=False,
            ttree_dict={}, qtree_dict={})
        self.assertTrue(filter_results[0] == (True, False, False), 
                        f'Failed when reading {fn} [0]')
        self.assertTrue(not any(filter_results[1]), f'Failed when reading {fn} [1]')
    
    def test_read_forward_filter_unique_source(self):
        fn = 'testdata/forward.chain'
        stree_dict = {'chr7': intervaltree.IntervalTree(),
                      'chr19': intervaltree.IntervalTree()}
        stree_dict['chr7'][60195160: 60195161] = 1
        filter_results = self.read_filter_compare(
            fn=fn, segment_size=0, unique=True,
            ttree_dict=stree_dict, qtree_dict={})
        self.assertTrue(filter_results[0] == (False, True, False), 
                        f'Failed when reading {fn} [0]')
        self.assertTrue(not any(filter_results[1]), f'Failed when reading {fn} [1]')
    
    def test_read_forward_filter_unique_target(self):
        fn = 'testdata/forward.chain'
        ttree_dict = {'chr7': intervaltree.IntervalTree(),
                      'chr5': intervaltree.IntervalTree()}
        ttree_dict['chr5'][50042645: 50042646] = 1
        filter_results = self.read_filter_compare(
            fn=fn, segment_size=0, unique=True,
            ttree_dict={}, qtree_dict=ttree_dict)
        self.assertTrue(not any(filter_results[0]), f'Failed when reading {fn} [0]')
        self.assertTrue(filter_results[1] == (False, False, True), 
                        f'Failed when reading {fn} [1]')


if __name__ == '__main__':
    unittest.main()
