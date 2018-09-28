#!/usr/bin/env python
import os
import sys
import json
import unittest

from utils import load_results

class TestAssemblyJson(unittest.TestCase):
   
    def setUp(self):
        """dummy pass"""
        benchmark_json = "ecoli/benchmark/asm_stats.json"
        result_json = "ecoli/benchmark/asm_stats.json"
        self.benchmark = load_results(benchmark_json)
        self.result = load_results(result_json)
        

    def test_asm_contigs(self):
        self.assertEqual(self.benchmark['asm_contigs'], self.result['asm_contigs'])

    def test_asm_esize(self):
        self.assertEqual(self.benchmark['asm_esize'], self.result['asm_esize'])

    def test_asm_mean(self):
        self.assertEqual(self.benchmark['asm_mean'], self.result['asm_mean'])

    def test_asm_median(self):
        self.assertEqual(self.benchmark['asm_median'], self.result['asm_median'])

    def test_asm_max(self):
        self.assertEqual(self.benchmark['asm_max'], self.result['asm_max'])

    def test_asm_min(self):
        self.assertEqual(self.benchmark['asm_min'], self.result['asm_min'])

    def test_asm_n50(self):
        self.assertEqual(self.benchmark['asm_n50'], self.result['asm_n50'])

    def test_asm_n90(self):
        self.assertEqual(self.benchmark['asm_n90'], self.result['asm_n90'])

    def test_asm_n95(self):
        self.assertEqual(self.benchmark['asm_n95'], self.result['asm_n95'])

    def test_asm_total_bp(self):
        self.assertEqual(self.benchmark['asm_total_bp'], self.result['asm_total_bp'])

if __name__ == "__main__":
    unittest.main()

