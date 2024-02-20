# -*- coding: utf-8 -*-
#########################################################################
# File Name: Part-2-Generate_PRE_fa.py
# Created on : 2022-12-26 13:49:24
# Author: JFF
# Last Modified: 2024-02-20 14:01:00
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import time
import re
import os
import sys
# sys.path.insert(0, ".")
from eQTac.generate_mut_fa import generate_mut_fa
from eQTac.generate_PRE import generate_PRE
from eQTac.generate_snp_dict import generate_snp_dict
import textwrap
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawTextHelpFormatter


# https://nibes.cn/blog/5845 showed line breaks in help text


class CustomArgumentFormatter(ArgumentDefaultsHelpFormatter, RawTextHelpFormatter):
    """Formats argument help which maintains line length restrictions as well as appends default value if present."""

    def _split_lines(self, text, width):
        text = super()._split_lines(text, width)
        new_text = []

        # loop through all the lines to create the correct wrapping for each line segment.
        for line in text:
            if not line:
                # this would be a new line.
                new_text.append(line)
                continue

            # wrap the line's help segment which preserves new lines but ensures line lengths are
            # honored
            new_text.extend(textwrap.wrap(line, width))

        return new_text


parser = ArgumentParser(
    formatter_class=CustomArgumentFormatter,
    description="Generate PRE fasta file.\nInput: PRE bed, plink genotype prefix (to estimate ld), SNP list file, whole genome fasta file\nOutput: xxx.ld_info, xxx.ld_info.snplist.bed.mutate.fa"
)
parser.add_argument('-pre',
                    '--PRE_bed',
                    help="Potential regulatory region .bed file. chr{1-22}\\tpos1\\tpos2",
                    type=str,
                    required=True)
parser.add_argument('-geno',
                    '--geno_prefix',
                    help="Genotype file prefix in plink format (prefix.bed, prefix.bim, prefix.fam).",
                    type=str,
                    required=True)
parser.add_argument('-snp',
                    '--snp_list',
                    help="SNP list file used in eQTac analysis. Note: only single nucleotide mutations.",
                    type=str,
                    required=True)
parser.add_argument(
    '-fa', '--fasta', help="Genome fasta file. Must with .fa.fai (e.g. xxx.fa & xxx.fa.fai)", type=str, required=True)
parser.add_argument('-r2', '--r2_max', help="Maximum ld between PRE SNPs",
                    type=float, required=False, default=0.3)
parser.add_argument('-d', '--distance_min',
                    help="Minimum distance bewtween PRE SNPs", type=int, required=False, default=10)
parser.add_argument('-o',
                    '--outfolder',
                    help="Output folder. Default eQTac_train_outfolder",
                    type=str,
                    default="eQTac_train_outfolder",
                    required=False)

args = parser.parse_args()
pre_bed = args.PRE_bed
geno_prefix = args.geno_prefix
snp_list = args.snp_list
fasta = args.fasta
r2_max = args.r2_max
distance_min = args.distance_min
outfolder = args.outfolder

t0 = time.time()
# Part2. Generate PRE fasta file
# 4. Generate PRE regions
print("#---- EQTac STEP4: Generate PRE regions START. ----#")

d_snp = generate_snp_dict(geno_prefix + ".bim")
ld_info = generate_PRE(geno_prefix, pre_bed, snp_list,
                       d_snp, outfolder, ld_max=r2_max, distance_min=distance_min)

print("#---- EQTac STEP4: Generate PRE regions FINISHED: %.6f. ----#" %
      (time.time() - t0))

# 5. Generate mutated fasta file
print("#---- EQTac STEP5: Generate mutated fasta file START. ----#")

mutate_fa = generate_mut_fa(ld_info, fasta, d_snp)

print("#---- EQTac STEP5: Generate mutated fasta file FINISHED: %.6f. ----#" %
      (time.time() - t0))
