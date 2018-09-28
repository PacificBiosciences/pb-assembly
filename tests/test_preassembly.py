#!/usr/bin/env python
import os
import sys
import json
import unittest
from utils import load_results

class TestPreassemblyJson(unittest.TestCase):
   
    def setUp(self):
        """dummy pass"""
        benchmark_json = "ecoli/benchmark/pre_assembly_stats.json"
        result_json = "ecoli/benchmark/pre_assembly_stats.json"
        self.benchmark = load_results(benchmark_json)
        self.result = load_results(result_json)

    def test_genome_length(self):
        self.assertEqual(self.benchmark['genome_length'], self.result['genome_length'])

    def test_preassembled_bases(self):
        self.assertEqual(self.benchmark['preassembled_bases'], self.result['preassembled_bases'])

    def test_preassembled_coverage(self):
        self.assertEqual(self.benchmark['preassembled_coverage'], self.result['preassembled_coverage'])

    def test_preassembled_esize(self):
        self.assertEqual(self.benchmark['preassembled_esize'], self.result['preassembled_esize'])
 
    def test_preassembled_mean(self):
        self.assertEqual(self.benchmark['preassembled_mean'], self.result['preassembled_mean'])

    def test_preassembled_n50(self):
        self.assertEqual(self.benchmark['preassembled_n50'], self.result['preassembled_n50'])

    def test_preassembled_p95(self):
        self.assertEqual(self.benchmark['preassembled_p95'], self.result['preassembled_p95'])

    def test_preassembled_reads(self):
        self.assertEqual(self.benchmark['preassembled_reads'], self.result['preassembled_reads'])

    def test_preassembled_seed_fragmentation(self):
        self.assertEqual(self.benchmark['preassembled_seed_fragmentation'], self.result['preassembled_seed_fragmentation'])

    def test_preassembled_seed_truncation(self):
        self.assertEqual(self.benchmark['preassembled_seed_truncation'], self.result['preassembled_seed_truncation'])

    def test_preassembled_yield(self):
        self.assertEqual(self.benchmark['preassembled_yield'], self.result['preassembled_yield'])

    def test_raw_bases(self):
        self.assertEqual(self.benchmark['raw_bases'], self.result['raw_bases'])

    def test_raw_coverage(self):
        self.assertEqual(self.benchmark['raw_coverage'], self.result['raw_coverage'])

    def test_raw_esize(self):
        self.assertEqual(self.benchmark['raw_esize'], self.result['raw_esize'])

    def test_raw_mean(self):
        self.assertEqual(self.benchmark['raw_mean'], self.result['raw_mean'])

    def test_raw_n50(self):
        self.assertEqual(self.benchmark['raw_n50'], self.result['raw_n50'])

    def test_raw_p95(self):
        self.assertEqual(self.benchmark['raw_p95'], self.result['raw_p95'])

    def test_raw_reads(self):
        self.assertEqual(self.benchmark['raw_reads'], self.result['raw_reads'])

    def test_seed_bases(self):
        self.assertEqual(self.benchmark['seed_bases'], self.result['seed_bases'])

    def test_seed_coverage(self):
        self.assertEqual(self.benchmark['seed_coverage'], self.result['seed_coverage'])

    def test_seed_esize(self):
        self.assertEqual(self.benchmark['seed_esize'], self.result['seed_esize'])

    def test_seed_mean(self):
        self.assertEqual(self.benchmark['seed_mean'], self.result['seed_mean'])

    def test_seed_n50(self):
        self.assertEqual(self.benchmark['seed_n50'], self.result['seed_n50'])

    def test_seed_p95(self):
        self.assertEqual(self.benchmark['seed_p95'], self.result['seed_p95'])

    def test_seed_reads(self):
        self.assertEqual(self.benchmark['seed_reads'], self.result['seed_reads'])

if __name__ == "__main__":
    unittest.main()

