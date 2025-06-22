"""
Tests for chaintools_bio utilities

Nae-Chyun Chen
Johns Hopkins University

Nancy Fisher Hansen
NIH/NHGRI

2021-2022
"""

import unittest

import intervaltree

from chaintools_bio import (
    chain_filter,
    split,
    to_bed,
    to_paf,
    to_sam,
    to_vcf,
    utils,
)


class TestUtils(unittest.TestCase):
    def test_reverse_complement(self):
        self.assertEqual(utils.reverse_complement("acg"), "CGT")
        self.assertEqual(utils.reverse_complement("ACG"), "CGT")
        self.assertEqual(utils.reverse_complement("TCGAN"), "NTCGA")


class TestChain(unittest.TestCase):
    """Tests of the Chain class."""

    def test_update_cigar_indel_simple(self):
        c = utils.Chain()
        c.cigarstring = ""
        c.nm_tag_val = 0
        c.strand = "+"
        queryendpos = c.update_cigar_indel(deltaq=0, deltat=1)
        self.assertEqual(c.cigarstring, "1D")
        self.assertEqual(queryendpos, 0)
        queryendpos = c.update_cigar_indel(deltaq=2, deltat=0)
        self.assertEqual(c.cigarstring, "1D2I")
        self.assertEqual(queryendpos, 2)

    def test_update_cigar_indel_breakpoint(self):
        c = utils.Chain()
        c.cigarstring = ""
        c.nm_tag_val = 0
        c.strand = "+"
        queryendpos = c.update_cigar_indel(deltaq=2000, deltat=300)
        self.assertEqual(c.cigarstring, "2000I300D")
        self.assertEqual(queryendpos, 2000)


class TestReadingChain(unittest.TestCase):
    """Tests reading a chain file."""

    def read_and_compare(self, fn):
        """Reads a chain file using utils.Chain().

        The processed text (`output_txt`) should equal to the input text
        (`input_txt`).
        """
        input_txt = ""
        output_txt = ""
        with open(fn, "r") as f:
            for line in f:
                input_txt += line
                fields = line.split()
                if len(fields) == 0:
                    continue
                elif line.startswith("chain"):
                    c = utils.Chain(fields)
                else:
                    c.add_record(fields)
                    if len(fields) == 1:
                        output_txt += c.print_chain() + "\n"
                        c = None

        # A tab and a space are considered to be indifferent
        input_txt = input_txt.replace("\t", " ").rstrip()
        output_txt = output_txt.replace("\t", " ").rstrip()
        ii = input_txt.split("\n")
        oo = output_txt.split("\n")
        if len(ii) != len(oo):
            raise ValueError(f"Lengths differ ({len(ii)} vs. {len(oo)})")
        for i in range(len(ii)):
            if ii[i] != oo[i]:
                raise ValueError(f"Context differ ({ii[i]} != {oo[i]})")
        self.assertSequenceEqual(input_txt, output_txt)

    def test_read_chain(self):
        for fn in [
            "testdata/forward.chain",
            "testdata/reversed.chain",
            "testdata/end_with_zero.chain",
            "testdata/small.chain",
        ]:
            with self.subTest(msg=f"Failed when reading {fn}"):
                self.read_and_compare(fn)


class TestGenerateVcf(unittest.TestCase):
    def generate_and_check(self, chainfn, targetfn, queryfn, vcffn):
        targetref = utils.fasta_reader(targetfn)
        queryref = utils.fasta_reader(queryfn)
        output_txt = utils.vcf_header(utils.get_target_entries(chainfn)) + "\n"
        with open(chainfn, "r") as f:
            out = to_vcf.write_to_vcf(f, targetref, queryref)
            while True:
                try:
                    output_txt += next(out) + "\n"
                except StopIteration:
                    break
        with open(vcffn, "r") as f:
            vcf_txt = ""
            output_lines = output_txt.split("\n")
            for i, line in enumerate(f):
                if line.rstrip() != output_lines[i]:
                    print("\n   truth vcf line:", line.rstrip())
                    print("reported vcf line:", output_lines[i])
                    return False
                vcf_txt += line

        return output_txt == vcf_txt

    def test_generate_vcf_from_small(self):
        fn = "testdata/target-query.chain"
        targetfn = "testdata/target.fasta"
        queryfn = "testdata/query.fasta"
        vcffn = "testdata/target-query.vcf"

        self.assertTrue(
            self.generate_and_check(fn, targetfn, queryfn, vcffn),
            f"Failed when generating vcf from {fn}",
        )


class TestGenerateSAM(unittest.TestCase):
    def generate_and_check(self, chainfn, targetfn, queryfn, samfn):
        output_txt = ""
        targetref = utils.fasta_reader(targetfn)
        queryref = utils.fasta_reader(queryfn)
        output_txt += utils.sam_header(utils.get_target_entries(chainfn)) + "\n"
        with open(chainfn, "r") as f:
            out = to_sam.write_to_sam(f, targetref, queryref)
            while True:
                try:
                    output_txt += next(out) + "\n"
                except StopIteration:
                    break
        with open(samfn, "r") as f:
            sam_txt = ""
            for line in f:
                sam_txt += line
        return output_txt == sam_txt

    def test_generate_sam_from_small(self):
        fn = "testdata/target-query.chain"
        targetfn = "testdata/target.fasta"
        queryfn = "testdata/query.fasta"
        samfn = "testdata/target-query.sam"

        self.assertTrue(
            self.generate_and_check(fn, targetfn, queryfn, samfn),
            f"Failed when generating sam from {fn}",
        )


class TestPaf(unittest.TestCase):
    def generate_and_check(self, chainfn, targetfn, queryfn, paffn):
        targetref = utils.fasta_reader(targetfn)
        queryref = utils.fasta_reader(queryfn)
        output_txt = ""
        with open(chainfn, "r") as f:
            out = to_paf.write_to_paf(f, targetref, queryref)
            while True:
                try:
                    output_txt += next(out) + "\n"
                except StopIteration:
                    break
        with open(paffn, "r") as f:
            paf_txt = ""
            for line in f:
                paf_txt += line
        self.assertSequenceEqual(output_txt, paf_txt)

    def test_paf_read_write_integrity(self):
        for test_case in [
            (
                "testdata/target-query.chain",
                "testdata/target.fasta",
                "testdata/query.fasta",
                "testdata/target-query.paf",
            ),
            (
                "testdata/target-query.chain",
                "",
                "",
                "testdata/target-query.no_ref.paf",
            ),
            ("testdata/forward.chain", "", "", "testdata/forward.no_ref.paf"),
            ("testdata/reversed.chain", "", "", "testdata/reversed.no_ref.paf"),
        ]:
            with self.subTest(msg=f"Failed when reading {test_case[0]}"):
                self.generate_and_check(
                    chainfn=test_case[0],
                    targetfn=test_case[1],
                    queryfn=test_case[2],
                    paffn=test_case[3],
                )


class TestGenerateBED(unittest.TestCase):
    def generate_and_check(self, chainfn, bedfn, coord):
        output_txt = ""
        with open(chainfn, "r") as f:
            out = to_bed.write_to_bed(f, coord)
            while True:
                try:
                    output_txt += next(out) + "\n"
                except StopIteration:
                    break
        with open(bedfn, "r") as f:
            bed_txt = ""
            for line in f:
                bed_txt += line
        self.assertSequenceEqual(output_txt, bed_txt)

    def test_bed_read_write_integrity(self):
        for coord in ["target", "query"]:
            for test_case_name in ["target-query", "forward", "reversed"]:
                with self.subTest(msg=f"Failed: {test_case_name}"):
                    self.generate_and_check(
                        chainfn=f"testdata/{test_case_name}.chain",
                        bedfn=f"testdata/{test_case_name}.{coord}.bed",
                        coord=coord,
                    )


class TestFilter(unittest.TestCase):
    def read_filter_compare(
        self,
        fn: str,
        segment_size: int,
        unique: bool,
        ttree_dict: dict,
        qtree_dict: dict,
    ) -> list:
        filter_results = []
        with open(fn, "r") as f:
            for line in f:
                fields = line.split()
                if len(fields) == 0:
                    continue
                elif line.startswith("chain"):
                    c = utils.Chain(fields)
                elif len(fields) == 3:
                    c.add_record(fields)
                elif len(fields) == 1:
                    c.add_record(fields)
                    filter_results.append(
                        chain_filter.filter_core(
                            c=c,
                            segment_size=segment_size,
                            unique=unique,
                            ttree_dict=ttree_dict,
                            qtree_dict=qtree_dict,
                        )
                    )
                    c = None
        return filter_results

    def test_read_forward_no_filter(self):
        fn = "testdata/forward.chain"
        filter_results = self.read_filter_compare(
            fn=fn, segment_size=0, unique=False, ttree_dict={}, qtree_dict={}
        )
        self.assertTrue(
            not any(filter_results[0]), f"Failed when reading {fn} [0]"
        )
        self.assertTrue(
            not any(filter_results[1]), f"Failed when reading {fn} [1]"
        )

    def test_read_forward_filter_size(self):
        fn = "testdata/forward.chain"
        filter_results = self.read_filter_compare(
            fn=fn,
            segment_size=360000,
            unique=False,
            ttree_dict={},
            qtree_dict={},
        )
        self.assertTrue(
            filter_results[0] == (True, False, False),
            f"Failed when reading {fn} [0]",
        )
        self.assertTrue(
            not any(filter_results[1]), f"Failed when reading {fn} [1]"
        )

    def test_read_forward_filter_unique_source(self):
        fn = "testdata/forward.chain"
        stree_dict = {
            "chr7": intervaltree.IntervalTree(),
            "chr19": intervaltree.IntervalTree(),
        }
        stree_dict["chr7"][60195160:60195161] = 1
        filter_results = self.read_filter_compare(
            fn=fn,
            segment_size=0,
            unique=True,
            ttree_dict=stree_dict,
            qtree_dict={},
        )
        self.assertTrue(
            filter_results[0] == (False, True, False),
            f"Failed when reading {fn} [0]",
        )
        self.assertTrue(
            not any(filter_results[1]), f"Failed when reading {fn} [1]"
        )

    def test_read_forward_filter_unique_target(self):
        fn = "testdata/forward.chain"
        ttree_dict = {
            "chr7": intervaltree.IntervalTree(),
            "chr5": intervaltree.IntervalTree(),
        }
        ttree_dict["chr5"][50042645:50042646] = 1
        filter_results = self.read_filter_compare(
            fn=fn,
            segment_size=0,
            unique=True,
            ttree_dict={},
            qtree_dict=ttree_dict,
        )
        self.assertTrue(
            not any(filter_results[0]), f"Failed when reading {fn} [0]"
        )
        self.assertTrue(
            filter_results[1] == (False, False, True),
            f"Failed when reading {fn} [1]",
        )


class TestSplit(unittest.TestCase):
    def generate_and_check(self, chainfn, splitfn, min_bp, min_gap):
        output_txt = ""
        with open(chainfn, "r") as f:
            out = split.split_chain(f=f, min_bp=min_bp, min_gap=min_gap)
            while True:
                try:
                    output_txt += next(out) + "\n"
                except StopIteration:
                    break
        with open(splitfn, "r") as f:
            split_txt = ""
            for line in f:
                split_txt += line
        oo = output_txt.split("\n")
        for i, l in enumerate(split_txt.split("\n")):
            if l != oo[i]:
                raise ValueError(
                    f"Lines are different: input={l}, output={oo[i]}"
                )
        self.assertSequenceEqual(output_txt, split_txt)
        return output_txt == split_txt

    def test_bed_read_write_integrity(self):
        for test_case_name in ["forward", "reversed", "reversed2"]:
            with self.subTest(msg=f"Failed: {test_case_name}"):
                self.generate_and_check(
                    chainfn=f"testdata/{test_case_name}.chain",
                    splitfn=f"testdata/{test_case_name}-split.chain",
                    min_bp=1000,
                    min_gap=10000,
                )

    def test_check_split_min_bp(self):
        self.assertTrue(
            split.check_split(dt=101, dq=90, min_bp=89, min_gap=1000),
            f"Failed: chain_filter::check_split()",
        )
        self.assertFalse(
            split.check_split(dt=101, dq=90, min_bp=90, min_gap=1000),
            f"Failed: chain_filter::check_split()",
        )

    def test_check_split_min_gap(self):
        self.assertTrue(
            split.check_split(dt=1001, dq=0, min_bp=89, min_gap=1000),
            f"Failed: chain_filter::check_split()",
        )
        self.assertFalse(
            split.check_split(dt=10, dq=0, min_bp=0, min_gap=1000),
            f"Failed: chain_filter::check_split()",
        )


if __name__ == "__main__":
    unittest.main()
