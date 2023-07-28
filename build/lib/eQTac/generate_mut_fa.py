# -*- coding: utf-8 -*-
#########################################################################
# File Name: generate_mut_fa.py
# Created on : 2022-12-19 17:02:20
# Author: JFF
# Last Modified: 2022-12-19 17:02:21
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


def generate_mut_fa(ld_info, genome_fasta, d_snp):
    # Generate snp bed file
    snp_bed = ld_info + ".snplist.bed"
    snpset = set([])
    with open(ld_info) as f2:
        next(f2)
        for i in f2:
            i = set(i.strip().split("\t")[11].split(","))
            snpset = snpset | i
    with open(snp_bed, 'w') as ff1:
        for snp in snpset:
            ch, pos, a1, a2 = d_snp[snp]
            if len(a1) == 1 and len(a2) == 1:  # only SNPs
                ff1.write("\t".join([ch, str(int(pos) - 10), str(int(pos) + 9), snp, a1, a2]) + "\n")
    # Generate fa file
    snp_fa = snp_bed + ".fa"
    os.system("bedtools getfasta -fi %(genome_fasta)s -bed %(snp_bed)s -name |seqkit seq -u > %(snp_fa)s" % {
        "genome_fasta": genome_fasta,
        "snp_bed": snp_bed,
        "snp_fa": snp_fa
    })
    # Generate mutate fa file
    d = {}
    with open(snp_bed) as f1:
        for i in f1:
            i = i.strip().split("\t")
            d.setdefault(i[3], {})[i[4]] = i[5]
            d.setdefault(i[3], {})[i[5]] = i[4]  # rsxxx:{A:T,T:A}
    mutate_fa = re.sub("\.fa$", "", snp_fa) + ".mutate.fa"
    with open(snp_fa) as f1:
        with open(mutate_fa, 'w') as ff:
            l_save = []
            n = 0
            for i in f1:
                i = i.strip()
                n += 1
                if n % 2 == 1:
                    l_save.append([])
                    n_pair = n // 2
                    l_save[n_pair].append(i)
                elif n % 2 == 0:
                    n_pair = n // 2 - 1
                    l_save[n_pair].append(i)
            for i in l_save:
                rs = re.sub("^>", "", i[0].split("::")[0])
                raw_allele = i[1][9]
                mutate_allele = d[rs][raw_allele]
                raw_title = i[0] + ":" + raw_allele + "_" + mutate_allele + ":" + raw_allele
                mutate_title = i[0] + ":" + raw_allele + "_" + mutate_allele + ":" + mutate_allele
                raw_seq = i[1]
                mutate_seq = i[1][:9] + mutate_allele + i[1][10:]
                ff.write("\n".join([raw_title, raw_seq]) + "\n")
                ff.write("\n".join([mutate_title, mutate_seq]) + "\n")
    return mutate_fa