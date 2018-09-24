#!/usr/bin/env python
"""Save current falcon asm results"""
import os
import sys
import shutil
import argparse
import falcon_kit
import falcon_unzip


def main():
    """Set benchmark data"""
    parser = get_parser()
    args = parser.parse_args()

    absolute_path = os.path.abspath(args.asm_dir)
    version_slug = "falcon_kit-{f}-falcon_unzip-{u}".format(
        f=falcon_kit.__version__, u=falcon_unzip.__version__)
    benchmark_dir = os.path.join(absolute_path, 'benchmark')

    if not os.path.exists(benchmark_dir):
        os.mkdir(benchmark_dir)

    new_pre_assembly_json = os.path.join(
        absolute_path, '0-rawreads/report/pre_assembly_stats.json')
    new_asm_json = os.path.join(absolute_path, '2-asm-falcon/asm_stats.json')

    link_pre_assembly_json = os.path.join(
        absolute_path, 'benchmark', 'pre_assembly_stats.json')
    link_asm_json = os.path.join(absolute_path, 'benchmark', 'asm_stats.json')

    archive_pre_assembly_json = os.path.join(
        absolute_path, 'benchmark', 'pre_assembly_stats_{d}.json'.format(d=version_slug))
    archive_asm_json = os.path.join(
        absolute_path, 'benchmark', 'asm_stats_{d}.json'.format(d=version_slug))

    if os.path.exists(link_pre_assembly_json):
        os.unlink(link_pre_assembly_json)
        print "Preassembly results for falcon_kit {d} already exist; overwriting...".format(
            d=version_slug)

    if os.path.exists(link_asm_json):
        print "Assembly results for falcon_kit {d} already exist; overwriting...".format(
            d=version_slug)
        os.unlink(link_asm_json)

    shutil.copyfile(new_pre_assembly_json, archive_pre_assembly_json)
    shutil.copyfile(new_asm_json, archive_asm_json)

    os.symlink(archive_pre_assembly_json, link_pre_assembly_json)
    os.symlink(archive_asm_json, link_asm_json)


def get_parser():
    """Get options"""

    parser = argparse.ArgumentParser()
    parser.add_argument('asm_dir', type=str, default='./',
                        help="path to a completed FALCON job")

    return parser


if __name__ == "__main__":
    sys.exit(main())
