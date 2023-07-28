# -*- coding: utf-8 -*-
#########################################################################
# File Name: filter_bkg.py
# Created on : 2022-12-17 16:43:57
# Author: JFF
# Last Modified: 2022-12-17 16:43:59
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
from pysam import FastxFile


def filter_bkg(negraw_bed, negraw_fa, exclude_bed):
    """
    Filter out the duplicate sequences and remove potential positive sequences from negative sequences.
    
    ## Parameters
    negraw_bed
        negraw_bed file from get_nullseq.
    negraw_fa
        negraw_fa file from get_nullseq.
    exclude_bed
        bed files contained potential positive regions. (e.g. peaks with P<0.05 from macs2).
    
    ## Returns
    out_neg_fa
        out_neg_fa can be used for trainning
    out_neg_bed
        out_neg_bed
    
    ## Notes
    Recommended at least setting ``exclude_bed`` as pos_bed to remove positive regions.
    
    """
    # out filename
    out_neg_bed = re.sub("\.negraw\.bed$", "", negraw_bed) + ".neg.bed"
    out_neg_fa = re.sub("\.negraw\.fa$", "", negraw_fa) + ".neg.fa"
    #used chr
    used_chrset = ["chr" + str(i) for i in list(range(1, 24))]
    pos_set = set([])
    name_set = set([])
    # find negraw bed intersect with exclude_bed
    with os.popen("awk 'BEGIN{OFS=\"\t\"}{print $0,$1\"_\"$2\"_\"$3\"_neg_\"NR}' %a|bedtools intersect -v -a - -b %s -wo " %
                  (negraw_bed, exclude_bed)) as f1:
        with open(out_neg_bed, 'w') as ff1:
            for i in f1:
                i = i.strip().split("\t")
                if i[0] in used_chrset:
                    name = i[3]  # chr3_159811111_159811310_neg_11
                    pos = "_".join(i[:3])  # chr3_159811111_159811310
                    if pos not in pos_set:  # 用位置判断是否重复
                        pos_set.add(pos)
                        name_set.add(name)
                        ff1.write("\t".join(i[:3]) + "\n")
    # write filtered fa to new file
    with FastxFile(negraw_fa) as f1:
        with open(out_neg_fa, 'w') as ff2:
            for i in f1:
                name = i.name
                seq = i.sequence
                if re.sub("^>", "", name) in name_set:
                    ff2.write(">" + name + "\n")
                    ff2.write(seq + "\n")
    return out_neg_fa, out_neg_bed


if __name__ == '__main__':
    filter_bkg("../Utilities_pipeline/example.output_eQTac/test.positive.negraw.bed", 
               "../Utilities_pipeline/example.output_eQTac/test.positive.negraw.fa",
               "../Utilities_pipeline/test_data/test.exclude.bed")
