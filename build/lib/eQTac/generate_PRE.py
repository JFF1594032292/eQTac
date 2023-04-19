# -*- coding: utf-8 -*-
#########################################################################
# File Name: generate_PRE.py
# Created on : 2022-12-19 10:19:00
# Author: JFF
# Last Modified: 2022-12-19 10:19:02
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
import time
import math
import numpy as np
from itertools import combinations
from shutil import copyfile


def generate_PRE(geno_prefix, pre_bed, snp_list, d_snp, outfolder):
    #make dir
    if not os.path.exists(outfolder):
        os.mkdir(outfolder)
    # get max pre length
    pre_length_set = set([])
    with open(pre_bed) as f1:
        for i in f1:
            i = i.strip().split("\t")
            pre_length_set.add(abs(int(i[2]) - int(i[1])))
    pre_max_length_kb = math.ceil(max(pre_length_set) / 1000)
    # generate snp bed
    snp_bed = outfolder + "/" + os.path.basename(snp_list) + ".bed"
    with open(snp_list) as f1:
        with open(snp_bed, 'w') as ff1:
            for i in f1:
                snp = i.strip()
                ch, pos, a1, a2 = d_snp[snp]
                if len(a1) == 1 and len(a2) == 1:
                    ff1.write("\t".join([ch, pos, pos, snp]) + "\n")
    # str format
    d_str = {
        "pre_bed": pre_bed,
        "pre_bed_name": os.path.basename(pre_bed),
        "snp_bed": snp_bed,
        "snp_bed_name": os.path.basename(snp_bed),
        "geno_prefix": geno_prefix,
        "outfolder": outfolder,
        "pre_max_length_kb": pre_max_length_kb
    }
    # 提取pre上的SNP，计算ld矩阵
    os.system(
        "bedtools intersect -a %(snp_bed)s -b %(pre_bed)s  |cut -f4|sort -u > %(outfolder)s/%(snp_bed_name)s--%(pre_bed_name)s.pre_snplist"
        % d_str)
    os.system(
        "plink --allow-no-sex --bfile %(geno_prefix)s --extract %(outfolder)s/%(snp_bed_name)s--%(pre_bed_name)s.pre_snplist --r2  --ld-window-r2 0 --ld-window-kb %(pre_max_length_kb)d --ld-window 999999 --out %(outfolder)s/%(snp_bed_name)s--%(pre_bed_name)s.pre_snplist"
        % d_str)
    # 提取ld信息
    d_ld = {}
    with open("%(outfolder)s/%(snp_bed_name)s--%(pre_bed_name)s.pre_snplist.ld" % d_str) as f1:
        next(f1)
        for i in f1:
            i = i.strip().split()
            snp1 = i[2]
            snp2 = i[5]
            r2 = float(i[6])
            d_ld.setdefault(snp1, {})[snp2] = r2
            d_ld.setdefault(snp2, {})[snp1] = r2  # snp2:{snp1:r2},snp1:{snp2:r2}
    # 获取SNP的位置
    d_pos = {}
    with open(snp_bed) as f1:
        for i in f1:
            i = i.strip().split("\t")
            d_pos[i[3]] = int(i[2])  # snp:pos
    # 计算平均r2并输出
    d_peak = {}
    with os.popen("bedtools intersect -a %(snp_bed)s -b %(pre_bed)s -wa -wb" % d_str) as f1:
        for i in f1:
            i = i.strip().split("\t")
            peak = i[4:8]
            if len(peak) == 3:  #没有名字的话加上名字
                peak = peak + ["PRE_" + "_".join(peak[:3])]
            peak = tuple(peak)
            snp = i[3]
            d_peak.setdefault(peak, set([])).add(snp)  # peak:set([snp,snp,])

    with open("%(outfolder)s/%(snp_bed_name)s--%(pre_bed_name)s.pre_snplist.ld_info" % d_str, 'w') as ff:
        ff.write("chr\tstart\tend\tPRE_region\tsnp_count\tmean_r2\tmin_r2\tr2\tmean_distance\tmin_distance\tdistance\tsnp_set\n")
        for peak in d_peak:
            snp_set = d_peak[peak]
            poss = [d_pos[j] for j in snp_set]
            distance_list = [abs(j[0] - j[1]) for j in combinations(poss, 2)]
            r2_list = []
            for snp_pair in combinations(snp_set, 2):
                #print(snp_pair[0], snp_pair[1])
                try:
                    r2 = d_ld[snp_pair[0]][snp_pair[1]]
                    r2_list.append(r2)
                except:
                    #print(snp_pair[0], snp_pair[1], "no ld information!")
                    continue
            # print(snp_set)
            # print(r2_list)
            if len(snp_set) > 1:
                r2_mean = np.mean(r2_list)
                r2_min = np.min(r2_list)
                min_dis = min(distance_list)
                mean_dis = np.mean(distance_list)
                if r2_min <= 0.3 and min_dis >= 10:
                    ff.write("\t".join(peak) + "\t" + "%d" % len(snp_set) + "\t" + "%.4g" % r2_mean + "\t" + "%.4g" % r2_min + "\t" +
                             ",".join(["%.3g" % i for i in r2_list]) + "\t" + "%.4g" % mean_dis + "\t" + "%.4g" % min_dis + "\t" +
                             ",".join(["%.3g" % i for i in distance_list]) + "\t" + ",".join(snp_set) + "\n")

    return "%(outfolder)s/%(snp_bed_name)s--%(pre_bed_name)s.pre_snplist.ld_info" % d_str
