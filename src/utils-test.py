'''
Tests for levioSAM utilities

Nae-Chyun Chen
Johns Hopkins University
2021
'''
import unittest
import subprocess
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
                    pass
                elif len(fields) == 1:
                    c.add_record_one(fields)
                    output_txt += c.print_chain()
                    c = None
        
        # A tab and a space are considered indifferent
        input_txt = input_txt.replace('\t', ' ').rstrip()
        output_txt = output_txt.replace('\t', ' ').rstrip()
        ii = input_txt.split('\n')
        oo = output_txt.split('\n')
        if len(ii) != len(oo):
            return False
        for i in range(min(len(ii), len(oo))):
            if ii[i] != oo[i]:
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


if __name__ == '__main__':
    unittest.main()
