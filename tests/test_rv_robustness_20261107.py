"""Synthetic limiting-case tests; not observational validation."""
import unittest
from research.rv_robustness_20261107 import weighted_mean,huber_location,leave_one_group_out

class RobustRVTests(unittest.TestCase):
    def test_constant(self):
        self.assertAlmostEqual(huber_location([4,4,4],[1,2,3]),4)
    def test_outlier_downweight(self):
        self.assertLess(abs(huber_location([0,0,0,20],[1]*4)),abs(weighted_mean([0,0,0,20],[1]*4)))
    def test_group_deletion(self):
        out=leave_one_group_out([1,1,5,5],[1]*4,['a','a','b','b'])
        self.assertAlmostEqual(out['a'],5)
        self.assertAlmostEqual(out['b'],1)
    def test_invalid(self):
        with self.assertRaises(ValueError): weighted_mean([1],[0])
        with self.assertRaises(ValueError): leave_one_group_out([1],[1],['a'])
if __name__ == "__main__": unittest.main()
