# -*- coding: utf-8 -*-
#########################################################################
# File Name: Part-4-Calculate_eQTac_correlation.py
# Created on : 2022-12-26 13:57:20
# Author: JFF
# Last Modified: 2023-09-22 22:01:53
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
parser = ArgumentParser(formatter_class=CustomArgumentFormatter, description="Calculates eQTac and controls FDR.")

parser.add_argument('-pre', '--PRE_score', help="PRE score file (xxx.PRE_score). From Part3", type=str, required=True)
parser.add_argument('-exp', '--expression_file', help="Gene expression file", type=str, required=True)
parser.add_argument('-u', '--upstream_bp', help="Upstream range of gene", type=int, default=1000000, required=False)
parser.add_argument('-d', '--downstream_bp', help="Downstream range of gene", type=int, default=1000000, required=False)
parser.add_argument('-n', '--n_permutation', help="Permutation counts.", type=int, required=True)
parser.add_argument('-o',
                    '--outfolder',
                    help="Output folder. Default eQTac_train_outfolder",
                    type=str,
                    default="eQTac_train_outfolder",
                    required=False)

args = parser.parse_args()
pre_score = args.PRE_score
exp_file = args.expression_file
upstream_bp = args.upstream_bp
downstream_bp = args.downstream_bp
n_permutation = args.n_permutation
outfolder = args.outfolder

t0 = time.time()
# make dir
if not os.path.exists(outfolder):
    os.mkdir(outfolder)
# Part4. Calculate correlation
# 8. Calculate eQTac correlation
print("#---- EQTac STEP8: Calculate eQTac correlation. ----#")

result_file = outfolder + "/" + os.path.basename(pre_score) + ".eQTac_result"
eQTac_correlation(pre_score, exp_file, result_file, upstream_bp, downstream_bp)

print("#---- EQTac STEP8: Calculate eQTac correlation: %.6f. ----#" % (time.time() - t0))

# 9. Permutation P values
print("#---- EQTac STEP9: Permutation P values. ----#")

permutation_plist_file = outfolder + "/" + os.path.basename(result_file) + ".permutation_plist"
eQTac_permutation(pre_score, exp_file, n_permutation, permutation_plist_file, upstream_bp, downstream_bp)

print("#---- EQTac STEP9: Permutation P values: %.6f. ----#" % (time.time() - t0))

# 10. Control FDR
print("#---- EQTac STEP10: Control FDR. ----#")

FDR_result_file = outfolder + "/" + os.path.basename(result_file) + ".FDR.txt"
control_FDR(result_file, permutation_plist_file, FDR_result_file)

print("#---- EQTac STEP10: Control FDR: %.6f. ----#" % (time.time() - t0))

print("")
print("#---- EQTac STEP ALL FINISHED: %.6f. ----#" % (time.time() - t0))
