'''
Tests for chaintools utilities

Nae-Chyun Chen
Johns Hopkins University

Nancy Fisher Hansen
NIH/NHGRI

2021-2022
'''
import sys
import unittest

import intervaltree

from chaintools import (chain_filter, split, to_bed, to_paf, to_sam, to_vcf,
                        utils)


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
                else:
                    c.add_record(fields)
                    if len(fields) == 1:
                        output_txt += (c.print_chain() + '\n')
                        c = None

        # A tab and a space are considered to be indifferent
        input_txt = input_txt.replace('\t', ' ').rstrip()
        output_txt = output_txt.replace('\t', ' ').rstrip()
        ii = input_txt.split('\n')
        oo = output_txt.split('\n')
        if len(ii) != len(oo):
            print(
                f'[E::read_and_compare] Lengths differ ({len(ii)} vs. {len(oo)})',
                file=sys.stderr)
            return False
        for i in range(min(len(ii), len(oo))):
            if ii[i] != oo[i]:
                print(f'[E::read_and_compare] {ii[i]} != {oo[i]}',
                      file=sys.stderr)
                return False
        return input_txt == output_txt

    def test_read_forward(self):
        fn = 'testdata/forward.chain'
        self.assertTrue(self.read_and_compare(fn), f'Failed when reading {fn}')

    def test_read_reversed(self):
        fn = 'testdata/reversed.chain'
        self.assertTrue(self.read_and_compare(fn), f'Failed when reading {fn}')

    def test_read_forward_end_with_zero(self):
        fn = 'testdata/end_with_zero.chain'
        self.assertTrue(self.read_and_compare(fn), f'Failed when reading {fn}')

    def test_read_forward_small(self):
        fn = 'testdata/small.chain'
        self.assertTrue(self.read_and_compare(fn), f'Failed when reading {fn}')


class TestGenerateVcf(unittest.TestCase):

    def generate_and_check(self, chainfn, targetfn, queryfn, vcffn):
        targetref = utils.fasta_reader(targetfn)
        queryref = utils.fasta_reader(queryfn)
        output_txt = (utils.vcf_header(utils.get_target_entries(chainfn)) +
                      '\n')
        with open(chainfn, 'r') as f:
            out = to_vcf.write_to_vcf(f, targetref, queryref)
            while True:
                try:
                    output_txt += (next(out) + '\n')
                except StopIteration:
                    break
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
        output_txt += (utils.sam_header(utils.get_target_entries(chainfn)) +
                       '\n')
        with open(chainfn, 'r') as f:
            out = to_sam.write_to_sam(f, targetref, queryref)
            while True:
                try:
                    output_txt += (next(out) + '\n')
                except StopIteration:
                    break
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
        targetref = utils.fasta_reader(targetfn)
        queryref = utils.fasta_reader(queryfn)
        output_txt = ''
        with open(chainfn, 'r') as f:
            out = to_paf.write_to_paf(f, targetref, queryref)
            while True:
                try:
                    output_txt += (next(out) + '\n')
                except StopIteration:
                    break
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
            out = to_bed.write_to_bed(f, coord)
            while True:
                try:
                    output_txt += (next(out) + '\n')
                except StopIteration:
                    break
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

    def test_generate_bed_from_forward_target(self):
        fn = 'testdata/forward.chain'
        bedfn = 'testdata/forward.target.bed'
        self.assertTrue(self.generate_and_check(fn, bedfn, 'target'),
                        f'Failed when generating BED (target) from {fn}')

    def test_generate_bed_from_forward_query(self):
        fn = 'testdata/forward.chain'
        bedfn = 'testdata/forward.query.bed'
        self.assertTrue(self.generate_and_check(fn, bedfn, 'query'),
                        f'Failed when generating BED (query) from {fn}')

    def test_generate_bed_from_reversed_target(self):
        fn = 'testdata/reversed.chain'
        bedfn = 'testdata/reversed.target.bed'
        self.assertTrue(self.generate_and_check(fn, bedfn, 'target'),
                        f'Failed when generating BED (target) from {fn}')

    def test_generate_bed_from_reversed_query(self):
        fn = 'testdata/reversed.chain'
        bedfn = 'testdata/reversed.query.bed'
        self.assertTrue(self.generate_and_check(fn, bedfn, 'query'),
                        f'Failed when generating BED (query) from {fn}')


class TestFilter(unittest.TestCase):

    def read_filter_compare(self, fn: str, segment_size: int, unique: bool,
                            ttree_dict: dict, qtree_dict: dict) -> list:
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
                    filter_results.append(
                        chain_filter.filter_core(c=c,
                                                 segment_size=segment_size,
                                                 unique=unique,
                                                 ttree_dict=ttree_dict,
                                                 qtree_dict=qtree_dict))
                    c = None
        return filter_results

    def test_read_forward_no_filter(self):
        fn = 'testdata/forward.chain'
        filter_results = self.read_filter_compare(fn=fn,
                                                  segment_size=0,
                                                  unique=False,
                                                  ttree_dict={},
                                                  qtree_dict={})
        self.assertTrue(not any(filter_results[0]),
                        f'Failed when reading {fn} [0]')
        self.assertTrue(not any(filter_results[1]),
                        f'Failed when reading {fn} [1]')

    def test_read_forward_filter_size(self):
        fn = 'testdata/forward.chain'
        filter_results = self.read_filter_compare(fn=fn,
                                                  segment_size=360000,
                                                  unique=False,
                                                  ttree_dict={},
                                                  qtree_dict={})
        self.assertTrue(filter_results[0] == (True, False, False),
                        f'Failed when reading {fn} [0]')
        self.assertTrue(not any(filter_results[1]),
                        f'Failed when reading {fn} [1]')

    def test_read_forward_filter_unique_source(self):
        fn = 'testdata/forward.chain'
        stree_dict = {
            'chr7': intervaltree.IntervalTree(),
            'chr19': intervaltree.IntervalTree()
        }
        stree_dict['chr7'][60195160:60195161] = 1
        filter_results = self.read_filter_compare(fn=fn,
                                                  segment_size=0,
                                                  unique=True,
                                                  ttree_dict=stree_dict,
                                                  qtree_dict={})
        self.assertTrue(filter_results[0] == (False, True, False),
                        f'Failed when reading {fn} [0]')
        self.assertTrue(not any(filter_results[1]),
                        f'Failed when reading {fn} [1]')

    def test_read_forward_filter_unique_target(self):
        fn = 'testdata/forward.chain'
        ttree_dict = {
            'chr7': intervaltree.IntervalTree(),
            'chr5': intervaltree.IntervalTree()
        }
        ttree_dict['chr5'][50042645:50042646] = 1
        filter_results = self.read_filter_compare(fn=fn,
                                                  segment_size=0,
                                                  unique=True,
                                                  ttree_dict={},
                                                  qtree_dict=ttree_dict)
        self.assertTrue(not any(filter_results[0]),
                        f'Failed when reading {fn} [0]')
        self.assertTrue(filter_results[1] == (False, False, True),
                        f'Failed when reading {fn} [1]')


class TestSplit(unittest.TestCase):

    def generate_and_check(self, chainfn, splitfn, min_bp, min_gap):
        output_txt = ''
        with open(chainfn, 'r') as f:
            out = split.split_chain(f=f, min_bp=min_bp, min_gap=min_gap)
            while True:
                try:
                    output_txt += (next(out) + '\n')
                except StopIteration:
                    break
        with open(splitfn, 'r') as f:
            split_txt = ''
            for line in f:
                split_txt += line
        oo = output_txt.split('\n')
        for i, l in enumerate(split_txt.split('\n')):
            if l != oo[i]:
                print(l)
                print(oo[i])
                return False
        return output_txt == split_txt

    def test_split_from_forward(self):
        fn = 'testdata/forward.chain'
        splitfn = 'testdata/forward-split.chain'
        self.assertTrue(
            self.generate_and_check(chainfn=fn,
                                    splitfn=splitfn,
                                    min_bp=1000,
                                    min_gap=10000),
            f'Failed when splitting {fn}')

    def test_split_from_reversed(self):
        fn = 'testdata/reversed.chain'
        splitfn = 'testdata/reversed-split.chain'
        self.assertTrue(
            self.generate_and_check(chainfn=fn,
                                    splitfn=splitfn,
                                    min_bp=1000,
                                    min_gap=10000),
            f'Failed when splitting {fn}')

    def test_split_from_reverse2(self):
        fn = 'testdata/reversed2.chain'
        splitfn = 'testdata/reversed2-split.chain'
        self.assertTrue(
            self.generate_and_check(chainfn=fn,
                                    splitfn=splitfn,
                                    min_bp=1000,
                                    min_gap=10000),
            f'Failed when splitting {fn}')

    def test_check_split_min_bp(self):
        self.assertTrue(
            split.check_split(dt=101, dq=90, min_bp=89, min_gap=1000),
            f'Failed: chain_filter::check_split()')
        self.assertFalse(
            split.check_split(dt=101, dq=90, min_bp=90, min_gap=1000),
            f'Failed: chain_filter::check_split()')

    def test_check_split_min_gap(self):
        self.assertTrue(
            split.check_split(dt=1001, dq=0, min_bp=89, min_gap=1000),
            f'Failed: chain_filter::check_split()')
        self.assertFalse(
            split.check_split(dt=10, dq=0, min_bp=0, min_gap=1000),
            f'Failed: chain_filter::check_split()')


if __name__ == '__main__':
    unittest.main()
