# -*- coding: utf-8 -*-
#########################################################################
# File Name: Part-3-Predict_PRE_score.py
# Created on : 2022-12-26 13:56:46
# Author: JFF
# Last Modified: 2024-02-20 14:01:31
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
# sys.path.insert(0, ".")
import re
import time
import textwrap
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawTextHelpFormatter
from eQTac.geno2score import geno2score


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


#parser = argparse.ArgumentParser(description="A pipeline for calculating eQTac.")
parser = ArgumentParser(
    formatter_class=CustomArgumentFormatter,
    description="Predcit PRE score for each individual.\nInput: plink genotype prefix (to estimate PRE score), svm model from Part1, SNP list file from Part2 (xxx.pre_snplist), mutate fasta file from Part2, ld info file from Part2.\nOutput: PRE score file (xxx.PRE_score)"
)
parser.add_argument('-m', '--model', help="SVM model from Part1. (xxx.svmmodel.{t}_{l}_{k}_{e}.model.txt)", type=str, required=True)
parser.add_argument('-l', '--ld_info', help="PRE ld infomation file. (xxx.ld_info)", type=str, required=True)
parser.add_argument('-mfa', '--mutate_fa', help="Mutation fa file. (xxx.ld_info.snplist.bed.mutate.fa)", type=str, required=True)
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
parser.add_argument('-n', '--normalize', help="If normalized PRE score to mean=0, std=1", default=True, type=bool, required=False)
parser.add_argument('-r2', '--r2_max', help="Maximum ld between PRE SNPs",
                    type=float, required=False, default=0.3)
parser.add_argument('-d', '--distance_min',
                    help="Minimum distance bewtween PRE SNPs", type=int, required=False, default=10)
parser.add_argument('-T',
                    required=False,
                    help="set the number of threads for parallel calculation, 1, 4, or 16\n",
                    default=1,
                    type=int,
                    choices=[1, 4, 16])
parser.add_argument('-o',
                    '--outfolder',
                    help="Output folder. Default eQTac_train_outfolder",
                    type=str,
                    default="eQTac_train_outfolder",
                    required=False)

args = parser.parse_args()
model = args.model
mutate_fa = args.mutate_fa
snp_list = args.snp_list
ld_info = args.ld_info
geno_prefix = args.geno_prefix
normalize = args.normalize
r2_max = args.r2_max
distance_min = args.distance_min
outfolder = args.outfolder
T = args.T

t0 = time.time()
# make dir
if not os.path.exists(outfolder):
    os.mkdir(outfolder)
# Part3. Predict PRE fasta score for individuals
# 6.Predict svm weights for each SNPs
print("#---- EQTac STEP6: Predict svm weights START. ----#")

pred_out = outfolder + "/" + re.sub("\.fa$", "", os.path.basename(mutate_fa)) + ".pred_out"
os.system(f"gkmpredict {mutate_fa} {model} {pred_out} -T {T}")

print("#---- EQTac STEP6: Predict svm weights FINISHED: %.6f. ----#" % (time.time() - t0))

# 7. Calculate scores for each individuals in each PRE
print("#---- EQTac STEP7: Calculate scores for each individuals in each PRE START. ----#")

geno_prefix_name = os.path.basename(geno_prefix)
if not os.path.exists(outfolder + "/" + geno_prefix_name + ".vcf.gz"):
    os.system(
        f"plink --allow-no-sex --bfile {geno_prefix} --extract {snp_list} --recode vcf-iid bgz --output-chr chr26 --out {outfolder}/{geno_prefix_name} && tabix -p vcf {outfolder}/{geno_prefix_name}.vcf.gz"
    )
PRE_scorefile = geno2score(f"{outfolder}/{geno_prefix_name}.vcf.gz", pred_out, ld_info, normalize, ld_max=r2_max, distance_min=distance_min)

print("#---- EQTac STEP7: Calculate scores for each individuals in each PRE FINISHED: %.6f. ----#" % (time.time() - t0))
