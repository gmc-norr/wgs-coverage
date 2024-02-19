#!/usr/bin/env python

import argparse
import sys

import pysam


def get_bam_seqs(bam):
    bamfile = pysam.AlignmentFile(bam, "rb")
    seqs = dict(zip(bamfile.header.references, bamfile.header.lengths))
    return seqs


def guess_genome(bam_seqs):
    chr1_len = bam_seqs.get("chr1")

    if chr1_len == 249250621:
        return "hg19"
    elif chr1_len == 248956422:
        return "hg38"

    raise ValueError("could not determine reference genome")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("bam")
    args = parser.parse_args()

    bam_seqs = get_bam_seqs(args.bam)
    genome_build = guess_genome(bam_seqs)

    print(genome_build, end="")


if __name__ == "__main__":
    main()
