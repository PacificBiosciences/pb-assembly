#!/usr/bin/env python
"""
@gconcepcion
Get FALCON assembly stats in JSON format
"""
import os
import sys
import json
import math
import argparse
from pbcore.io import FastaReader


def get_fasta_stats(fasta, genome_size):
    """Calculate basic fasta stats"""
    lengths = [len(record.sequence) for record in FastaReader(fasta)]
    lengths.sort(reverse=True)
    asm_contigs = len(lengths)
    asm_total_bp = sum(lengths)

    def get_nstat(lens, stat, genome_size=None):
        """Calculate all N* stats"""
        lens.sort(reverse=True)
        if genome_size is not None:
            total = genome_size
        else:
            total = sum(lens)
        limit = total * stat
        for num in lens:
            total -= num
            if total <= limit:
                return num

    asm_n50 = get_nstat(lengths, 0.50)
    asm_n90 = get_nstat(lengths, 0.10)
    asm_n95 = get_nstat(lengths, 0.05)
    
    asm_min = lengths[-1]
    asm_max = lengths[0]
    asm_mean = asm_total_bp / asm_contigs
    asm_median = int((lengths[int(math.floor(asm_contigs * .5))] +
                      lengths[int(math.floor(asm_contigs * .5))]) / 2)
    asm_esize = sum([x*x for x in lengths]) / asm_total_bp


    fasta_stats = {'asm_contigs': asm_contigs,
                   'asm_total_bp': asm_total_bp,
                   'asm_esize': asm_esize,
                   'asm_min': asm_min,
                   'asm_max': asm_max,
                   'asm_mean': asm_mean,
                   'asm_median': asm_median,
                   'asm_n50': asm_n50,
                   'asm_n90': asm_n90,
                   'asm_n95': asm_n95}

    if genome_size is not None:
        asm_ng50 = get_nstat(lengths, 0.50, genome_size)
        asm_ng90 = get_nstat(lengths, 0.10, genome_size)
        asm_ng95 = get_nstat(lengths, 0.05, genome_size)
    

        fasta_stats.update({'asm_ng50': asm_ng50,
                            'asm_ng90': asm_ng90,
                            'asm_ng95': asm_ng95})


    return fasta_stats


def main():
    """Run fasta stat script"""

    parser = get_parser()
    args = parser.parse_args()

    if os.path.isdir(args.asm_dir):
        fasta = os.path.join(args.asm_dir, 'p_ctg.fa')
    else:
        if args.asm_dir.endswith(('.fa', '.fasta')):
            fasta = args.asm_dir
        else:
            print "Please specify a fasta or a 2-asm-falcon dir!"
    if args.genome_size is not None:
        genome_size = args.genome_size
    else:
        genome_size = None

    fasta_stats = get_fasta_stats(fasta, genome_size)

    if args.output:
        if os.path.isdir(args.output):
            outfile = os.path.join(args.output, 'asm_stats.json')
        else:
            outfile = args.output
        with open(outfile, 'w') as outfile:

            json.dump(fasta_stats, outfile, indent=1, sort_keys=True)

    print json.dumps(fasta_stats, indent=1, sort_keys=True)


def get_parser():
    """Get options"""

    parser = argparse.ArgumentParser()
    parser.add_argument('--debug', action='store_true',
                        help="Print debug logging to stdout")
    parser.add_argument('asm_dir', type=str, default='./',
                        help="path to 2-asm-falcon directory or fasta file")
    parser.add_argument('--output', type=str,
                        help="path to outdir or filename")
    parser.add_argument('--genome_size', type=int,
                        help="if known; genome_size for NG50 calculation")

    return parser


if __name__ == "__main__":
    sys.exit(main())
