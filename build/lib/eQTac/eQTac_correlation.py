# -*- coding: utf-8 -*-
#########################################################################
# File Name: eQTac_correlation.py
# Created on : 2022-12-19 21:39:53
# Author: JFF
# Last Modified: 2022-12-19 21:39:54
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
from pybedtools import BedTool
from scipy.stats import linregress
import numpy as np
import pandas as pd


def st_linear_reg(exp, score):  #这个速度更快
    #select = ~np.isnan(d_score[pre])  #得分不为nan（基因型不缺失）的pre
    #slope, intercept, r_value, p, std_err = stats.linregress(score[select], exp[select])
    slope, intercept, r_value, p, std_err = linregress(score, exp)
    return (slope, std_err, p, intercept, r_value)  #beta,se,p


def extract_data(PRE_scorefile, exp_file):
    score_data = pd.read_csv(PRE_scorefile, sep="\t", header=0, index_col=3).T
    exp_data = pd.read_csv(exp_file, sep="\t", header=0, index_col=3).T
    #提取基因cis区域（1MB）
    gene_range = exp_data.loc[["#chr", "start"], :].T
    d_gene_info = {i: list(gene_range.loc[i, :]) for i in gene_range.T}  #chr,tss,保存gene位置信息，用于后面输出
    gene_range.loc[:, "end"] = gene_range.loc[:, "start"] + 1000000
    gene_range.loc[:, "start"] = gene_range.loc[:, "start"] - 1000000
    gene_range.loc[:, "ENSG"] = gene_range.index
    gene_range.loc[gene_range.loc[:, "start"] < 0, "start"] = 0
    gene_range = BedTool(gene_range.values.tolist())  #chr,start,end, ensg
    #提取pre区域
    pre_range = score_data.loc[["chr", "start", "end", "snp_count", "snps"], :].T
    d_pre_info = {i: list(pre_range.loc[i, :]) for i in pre_range.T}  #chr,start,end,snp_count,snps,保存pre位置信息，用于后面输出
    pre_range.loc[:, "pre"] = pre_range.index
    pre_range = BedTool(pre_range.values.tolist())  #chr,start,end,snp_count,snps,pre
    #获取每个基因对应的SNP
    d_cis = {}
    for i in gene_range.intersect(pre_range, wo="-wo"):
        d_cis.setdefault(i[3], []).append(i[9])  #ensg:[pre,pre]
    score_data = score_data.iloc[5:, :]
    exp_data = exp_data.iloc[5:, :]
    #统一样本顺序
    samplelist = list(set(score_data.index) & set(exp_data.index))
    exp_data = exp_data.loc[samplelist, :]
    score_data = score_data.loc[samplelist, :]
    return score_data, exp_data, d_gene_info, d_cis, d_pre_info


def calculate_correlation(score_data, exp_data, d_gene_info, d_cis, d_pre_info, ff):
    # 基因表达和PRE信息存成字典
    #snp score
    d_score = {}
    for pre in score_data.columns:
        res = np.array(score_data.loc[:, pre]).astype(float)
        d_score[pre] = res  #pre:[score,score,]
    #exp data
    d_exp = {}
    for ensg in exp_data.columns:
        d_exp[ensg] = np.array(exp_data.loc[:, ensg]).astype(float)  #ensg:[exp,exp,]
    # 计算回归方程并输出结果
    ff.write("ENSG\tchr\tTSS\tPRE\tpre_start\tpre_end\tbeta\tse\tp\tpre_snps\n")
    for ensg in d_cis:
        ch, tss = d_gene_info[ensg][:]
        for pre in d_cis[ensg]:
            #if pre in d_score:  #有些SNP两个型没有差异，没存到字典中
            score = d_score[pre]
            exp = d_exp[ensg]
            beta, se, p, intercept, r_value = st_linear_reg(exp, score)
            pre_start, pre_end, snp_count, snps = d_pre_info[pre][1:]
            ff.write("\t".join(
                [ensg, ch, "%d" %
                 tss, pre, "%d" %
                 pre_start, "%d" %
                 pre_end, "%.5g" %
                 beta, "%.5g" %
                 se, "%.5g" % p, snps]) + "\n")
        #print(ensg, "is done")


def eQTac_correlation(PRE_scorefile, exp_file, result_file):
    score_data, exp_data, d_gene_info, d_cis, d_pre_info = extract_data(PRE_scorefile, exp_file)
    with open(result_file, 'w') as ff:
        calculate_correlation(score_data, exp_data, d_gene_info, d_cis, d_pre_info, ff)


if __name__ == "__main__":
    eQTac_correlation("../Utilities_pipeline/example.output_eQTac/test.geno.vcf.gz.PRE_score",
                      "../Utilities_pipeline/test_data/test.exp_residual",
                      "../Utilities_pipeline/example.output_eQTac/test.geno.vcf.gz.PRE_score.eQTac_result")
