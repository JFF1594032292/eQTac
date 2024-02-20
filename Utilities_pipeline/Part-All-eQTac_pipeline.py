# -*- coding: utf-8 -*-
#########################################################################
# File Name: All-eQTac_pipeline.py
# Created on : 2022-12-17 11:56:36
# Author: JFF
# Last Modified: 2024-02-20 14:46:13
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
from eQTac.generate_snp_dict import generate_snp_dict
from eQTac.generate_PRE import generate_PRE
from eQTac.generate_mut_fa import generate_mut_fa
from eQTac.geno2score import geno2score
from eQTac.eQTac_correlation import eQTac_correlation
from eQTac.eQTac_permutation import eQTac_permutation
from eQTac.control_FDR import control_FDR


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
parser = ArgumentParser(formatter_class=CustomArgumentFormatter,
                        description="A pipeline for calculating eQTac.\nNote: Parameters after -o was from gkmtrain software.")
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
parser.add_argument('-fa', '--fasta', help="Genome fasta file. Must with .fa.fai (e.g. xxx.fa & xxx.fa.fai)", type=str, required=True)
parser.add_argument('-n', '--n_permutation', help="Permutation counts.", type=int, required=True)
parser.add_argument('-norm', '--normalize', help="If normalized PRE score to mean=0, std=1", default=True, type=bool, required=False)
parser.add_argument('-r2', '--r2_max', help="Maximum ld between PRE SNPs",
                            type = float, required = False, default = 0.3)
parser.add_argument('-dis', '--distance_min',
                            help="Minimum distance bewtween PRE SNPs", type=int, required=False, default=10)
parser.add_argument('-exp', '--expression_file', help="Gene expression file", type=str, required=True)
parser.add_argument('-o',
                    '--outfolder',
                    help="Output folder. Default eQTac_train_outfolder",
                    type=str,
                    default="eQTac_train_outfolder",
                    required=False)
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
pre_bed = args.PRE_bed
geno_prefix = args.geno_prefix
exclude_bed = args.exclude_bed
snp_list = args.snp_list
fasta = args.fasta
exp_file = args.expression_file
n_permutation = args.n_permutation
normalize = args.normalize
r2_max = args.r2_max
distance_min = args.distance_min
outfolder = args.outfolder
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
print("#---- EQTac STEP1: Get negative datasets START. ----#", flush=True)

pos_fa, pos_bed, negraw_fa, negraw_bed = get_nullseq(pos_bed_path, outfolder, xfold=2.5)

print("#---- EQTac STEP1: Get negative datasets FINISHED: %.6f. ----#" % (time.time() - t0), flush=True)

# 2. Fiter negative datasets
print("---- EQTac STEP2: Fiter negative datasets START. ----#", flush=True)

neg_fa, neg_bed = filter_bkg(negraw_bed, negraw_fa, exclude_bed)

print("#---- EQTac STEP2: Fiter negative datasets FINISHED: %.6f. ----#" % (time.time() - t0), flush=True)

# 3. Train model
print("#---- EQTac STEP3: Train model START. ----#", flush=True)

out_model = re.sub("\.bed$", "", pos_bed) + f".svmmodel.{t}_{l}_{k}_{e}"
os.system(
    f"gkmtrain {pos_fa} {neg_fa} {out_model} -t {t} -l {l} -k {k} -d {d} -g {g} -M {M} -H {H} -c {c} -e {e} -w {w} -m {m}  -r {r} -v {v} -T {T}"
)

print("#---- EQTac STEP3: Train model FINISHED: %.6f. ----#" % (time.time() - t0), flush=True)

# Part2. Generate PRE fasta file
# 4. Generate PRE regions
print("#---- EQTac STEP4: Generate PRE regions START. ----#", flush=True)

d_snp = generate_snp_dict(geno_prefix + ".bim")
ld_info = generate_PRE(geno_prefix, pre_bed, snp_list, d_snp, outfolder, ld_max=r2_max, distance_min=distance_min)

print("#---- EQTac STEP4: Generate PRE regions FINISHED: %.6f. ----#" % (time.time() - t0), flush=True)

# 5. Generate mutated fasta file
print("#---- EQTac STEP5: Generate mutated fasta file START. ----#", flush=True)

mutate_fa = generate_mut_fa(ld_info, fasta, d_snp)

print("#---- EQTac STEP5: Generate mutated fasta file FINISHED: %.6f. ----#" % (time.time() - t0), flush=True)

# Part3. Predict PRE fasta score for individuals
# 6.Predict svm weights for each SNPs
print("#---- EQTac STEP6: Predict svm weights START. ----#", flush=True)

pred_out = outfolder + "/" + re.sub("\.fa$", "", os.path.basename(mutate_fa)) + ".pred_out"
os.system(f"gkmpredict {mutate_fa} {out_model}.model.txt {pred_out} -T {T}")

print("#---- EQTac STEP6: Predict svm weights FINISHED: %.6f. ----#" % (time.time() - t0), flush=True)

# 7. Calculate scores for each individuals in each PRE
print("#---- EQTac STEP7: Calculate scores for each individuals in each PRE START. ----#", flush=True)

geno_prefix_name = os.path.basename(geno_prefix)
if not os.path.exists(outfolder + "/" + geno_prefix_name + ".vcf.gz"):
    os.system(
        f"plink --allow-no-sex --bfile {geno_prefix} --extract {snp_list} --recode vcf-iid bgz --output-chr chr26 --out {outfolder}/{geno_prefix_name} && tabix -p vcf {outfolder}/{geno_prefix_name}.vcf.gz"
    )
PRE_scorefile = geno2score(f"{outfolder}/{geno_prefix_name}.vcf.gz", pred_out, ld_info, normalize, ld_max=r2_max, distance_min=distance_min)

print("#---- EQTac STEP7: Calculate scores for each individuals in each PRE FINISHED: %.6f. ----#" % (time.time() - t0), flush=True)

# Part4. Calculate correlation
# 8. Calculate eQTac correlation
print("#---- EQTac STEP8: Calculate eQTac correlation. ----#", flush=True)

result_file = outfolder + "/" + os.path.basename(PRE_scorefile) + ".eQTac_result"
eQTac_correlation(PRE_scorefile, exp_file, result_file)

print("#---- EQTac STEP8: Calculate eQTac correlation: %.6f. ----#" % (time.time() - t0), flush=True)

# 9. Permutation P values
print("#---- EQTac STEP9: Permutation P values. ----#", flush=True)

permutation_plist_file = outfolder + "/" + os.path.basename(result_file) + ".permutation_plist"
eQTac_permutation(PRE_scorefile, exp_file, n_permutation, permutation_plist_file)

print("#---- EQTac STEP9: Permutation P values: %.6f. ----#" % (time.time() - t0), flush=True)

# 10. Control FDR
print("#---- EQTac STEP10: Control FDR. ----#", flush=True)

FDR_result_file = outfolder + "/" + os.path.basename(result_file) + ".FDR.txt"
control_FDR(result_file, permutation_plist_file, FDR_result_file)

print("#---- EQTac STEP10: Control FDR: %.6f. ----#" % (time.time() - t0), flush=True)

print("", flush=True)
print("#---- EQTac STEP ALL FINISHED: %.6f. ----#" % (time.time() - t0), flush=True)
