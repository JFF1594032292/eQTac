# -*- coding: utf-8 -*-
#########################################################################
# File Name: Part-1-Train_model.py
# Created on : 2022-12-26 13:39:41
# Author: JFF
# Last Modified: 2023-07-28 11:45:17
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
import time
import textwrap
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, RawTextHelpFormatter
from eQTac.get_nullseq import get_nullseq
from eQTac.filter_bkg import filter_bkg


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


#parser = argparse.ArgumentParser(description="Training SVM model to get weights for each kmers.")
parser = ArgumentParser(
    formatter_class=CustomArgumentFormatter,
    description="Training SVM model to get weights for each kmers.\nInput: positive bed file, exclude bed file\nOutput: xxx.svmmodel.{t}_{l}_{k}_{e}.model.txt.\nNote: Parameters after -o was from gkmtrain software."
)
parser.add_argument('-p',
                    '--positive_bed',
                    help="Positive .bed file for gkSVM trainning. chr{1-22}\\tpos1\\tpos2",
                    type=str,
                    required=True)
parser.add_argument(
    '-ex',
    '--exclude_bed',
    help="Exclude regions from negtive dataset. Usually be set to peaks called with P<0.05, or just set same with --positive_bed",
    type=str,
    required=True)
parser.add_argument('-o',
                    '--outfolder',
                    help="Output folder. Default eQTac_train_outfolder",
                    type=str,
                    default="eQTac_train_outfolder",
                    required=False)
parser.add_argument('-f', '--xfold', help="Background fold.", type=float, default=2.5, required=False)
parser.add_argument(
    '-t',
    required=False,
    help=textwrap.dedent(
        "set kernel function (default: 2 gkm).\nNOTE: RBF kernels (3 and 5) work best with -c 10 -g 2.\n0 -- gapped-kmer.\n1 -- estimated l-mer with full filter.\n2 -- estimated l-mer with truncated filter (gkm).\n3 -- gkm + RBF (gkmrbf).\n4 -- gkm + center weighted (wgkm).\n[weight = max(M, floor(M*exp(-ln(2)*D/H)+1))].\n5 -- gkm + center weighted + RBF (wgkmrbf)."
    ),
    default="2",
    type=str,
    choices=["1", "2", "3", "4", "5"])
parser.add_argument('-l', required=False, help="set word length, 3<=l<=12", default=11, type=int)
parser.add_argument('-k', required=False, help="set number of informative column, k<=l", default=7, type=int)
parser.add_argument('-d', required=False, help="set maximum number of mismatches to consider, d<=4", default=3, type=int)
parser.add_argument('-g', required=False, help="set gamma for RBF kernel. -t 3 or 5 only", default=1.0, type=float)
parser.add_argument('-M',
                    required=False,
                    help="set the initial value (M) of the exponential decay function\nfor wgkm-kernels. max=255, -t 4 or 5 only ",
                    default=50,
                    type=int)
parser.add_argument(
    '-H',
    required=False,
    help="set the half-life parameter (H) that is the distance (D) required\nto fall to half of its initial value in the exponential decay\nfunction for wgkm-kernels. -t 4 or 5 only ",
    default=50.0,
    type=float)
parser.add_argument('-c', required=False, help="set the regularization parameter SVM-C ", default=1.0, type=float)
parser.add_argument('-e', required=False, help="set the precision parameter epsilon ", default=0.001, type=float)
parser.add_argument('-w', required=False, help="set the parameter SVM-C to w*C for the positive set ", default=1.0, type=float)
parser.add_argument('-m',
                    required=False,
                    help="set cache memory size in MB\nNOTE: Large cache signifcantly reduces runtime. >4Gb is recommended",
                    default=100.0,
                    type=float)
#parser.add_argument('-x', required=False,help="set N-fold cross validation mode (default: no cross validation)",default=0,type=int)
#parser.add_argument('-i', required=False,help="run i-th cross validation only 1<=i<=ncv (default: all)",default=0,type=int)
parser.add_argument('-r', required=False, help="set random seed for shuffling in cross validation mode ", default=1, type=int)
parser.add_argument(
    '-v',
    required=False,
    help="set the level of verbosity (default: 2)\n0 -- error msgs only (ERROR)\n1 -- warning msgs (WARN)\n2 -- progress msgs at coarse-grained level (INFO)\n3 -- progress msgs at fine-grained level (DEBUG)\n4 -- progress msgs at finer-grained level (TRACE)",
    default="2",
    type=str,
    choices=["0", "1", "2", "3", "4"])
parser.add_argument('-T',
                    required=False,
                    help="set the number of threads for parallel calculation, 1, 4, or 16\n",
                    default=1,
                    type=int,
                    choices=[1, 4, 16])

args = parser.parse_args()
pos_bed_path = args.positive_bed
exclude_bed = args.exclude_bed
outfolder = args.outfolder
f = args.xfold
t = args.t
l = args.l
k = args.k
d = args.d
g = args.g
M = args.M
H = args.H
c = args.c
e = args.e
w = args.w
m = args.m
#x = args.x
#i = args.i
r = args.r
v = args.v
T = args.T

t0 = time.time()
# Part1. Train model

# 1. Get negative datasets
print("#---- EQTac STEP1: Get negative datasets START. ----#")

pos_fa, pos_bed, negraw_fa, negraw_bed = get_nullseq(pos_bed_path, outfolder, xfold=f)

print("#---- EQTac STEP1: Get negative datasets FINISHED: %.6f. ----#" % (time.time() - t0))

# 2. Fiter negative datasets
print("---- EQTac STEP2: Fiter negative datasets START. ----#")

neg_fa, neg_bed = filter_bkg(negraw_bed, negraw_fa, exclude_bed)

print("#---- EQTac STEP2: Fiter negative datasets FINISHED: %.6f. ----#" % (time.time() - t0))

# 3. Train model
print("#---- EQTac STEP3: Train model START. ----#")

out_model = re.sub("\.bed$", "", pos_bed) + f".svmmodel.{t}_{l}_{k}_{e}"
os.system(
    f"gkmtrain {pos_fa} {neg_fa} {out_model} -t {t} -l {l} -k {k} -d {d} -g {g} -M {M} -H {H} -c {c} -e {e} -w {w} -m {m}  -r {r} -v {v} -T {T}"
)

print("#---- EQTac STEP3: Train model FINISHED: %.6f. ----#" % (time.time() - t0))
