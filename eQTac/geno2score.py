# -*- coding: utf-8 -*-
#########################################################################
# File Name: geno2score.py
# Created on : 2022-12-19 20:30:37
# Author: JFF
# Last Modified: 2022-12-19 20:30:38
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
import numpy as np
import gzip


def geno2score(vcf_gz, pred_out, ld_info):
    PRE_scorefile = vcf_gz + ".PRE_score"
    # 增强子上的SNP
    d_pre = {}
    with open(ld_info) as f1:
        next(f1)
        for i in f1:
            i = i.strip().split("\t")
            if float(i[6]) < 0.3 and int(i[9]) > 10:  # enhancer上SNP之间最小的ld需要<0.3，SNP之间最小距离大于10bp
                enhancer = "_".join(i[:3])
                snp_set = i[11].split(",")
                d_pre[enhancer] = snp_set  # enhancer:set([snp,snp])
    # SNP得分
    d_score = {}
    with open(pred_out) as f2:
        for i in f2:
            i = i.strip().split("\t")
            snp = i[0].split(":")[0]
            allele = i[0].split(":")[-1]
            score = float(i[1])
            d_score.setdefault(snp, {})[allele] = score  # snp:{allele1:score,}

    # 将基因型转化为得分
    d_geno = {}
    miss_snp_set = set([])
    with gzip.open(vcf_gz, 'rt') as f3:
        for i in f3:
            i = i.strip().split("\t")
            if not i[0].startswith("##"):  # 注释行
                if i[0].startswith("#"):  # 标题行
                    header = i[9:]  # 保存标题行的样本顺序
                else:
                    snp = i[2]
                    if snp in d_score:  # snp有预测得分
                        geno_list = i[9:]
                        for j in range(len(header)):
                            geno = geno_list[j]  # 0/1,
                            alleles = re.sub("1", i[4], re.sub("0", i[3], geno)).split("/")  # 转成[A1,A2]
                            if "." in alleles:  # 基因型缺失，直接用“.”代替
                                score = "NA"
                                miss_snp_set.add(snp)  # 有缺失值的SNP
                            else:  # 基因型不缺失，计算得分
                                score = d_score[snp][alleles[0]] + d_score[snp][alleles[1]]
                            d_geno.setdefault(snp, {})[header[j]] = score  # snp:{OA-1:score,}
    # 填补缺失值（人群中的该SNP得分的均值）
    for snp in miss_snp_set:
        mean_score = np.mean([i for i in d_geno[snp].values() if i != "NA"])
        for sample in d_geno[snp]:
            if d_geno[snp][sample] == "NA":
                d_geno[snp][sample] = mean_score
    # 输出每个人每个enhancer的得分
    with open(PRE_scorefile, 'w') as ff:
        ff.write("\t".join(["chr", "start", "end", "PRE", "snp_count", 'snps'] + header) + "\n")
        for enhancer in d_pre:
            ch, start, end = enhancer.split("_")
            snp_set = [snp for snp in d_pre[enhancer] if snp in d_geno]
            snp_count = len(snp_set)
            score_list = ["%.6g" % sum([d_geno[snp][sample] for snp in snp_set]) for sample in header]
            snps = ",".join(snp_set)
            ff.write("\t".join([ch, start, end, enhancer, "%d" % snp_count, snps] + score_list) + "\n")
    return PRE_scorefile


if __name__ == "__main__":
    geno2score(
        "../Utilities_pipeline/example.output_eQTac/test.geno.vcf.gz",
        "../Utilities_pipeline/example.output_eQTac/test.geno.snplist.bed--test.pre.bed.pre_snplist.ld_info.snplist.bed.mutate.pred_out",
        "../Utilities_pipeline/example.output_eQTac/test.geno.snplist.bed--test.pre.bed.pre_snplist.ld_info"
    )
