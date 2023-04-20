# -*- coding: utf-8 -*-
#########################################################################
# File Name: eQTac_permutation.py
# Created on : 2022-12-20 15:58:46
# Author: JFF
# Last Modified: 2022-12-20 15:58:48
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
import numpy as np
from .eQTac_correlation import extract_data, st_linear_reg


def eQTac_permutation(PRE_scorefile, exp_file, permute_count, result_file):
    score_data, exp_data, d_gene_info, d_cis, d_pre_info = extract_data(PRE_scorefile, exp_file)
    with open(result_file, 'w') as ff:
        for seed in range(permute_count):
            exp_data_copy = exp_data.copy()
            # 随机置换表达值
            rng = np.random.default_rng(seed)
            exp_label = list(exp_data.index)
            rng.shuffle(exp_label)
            exp_data_copy = exp_data_copy.loc[exp_label, :]
            # 基因表达和PRE信息存成字典
            #snp score
            d_score = {}
            for pre in score_data.columns:
                res = np.array(score_data.loc[:, pre]).astype(float)
                d_score[pre] = res  #pre:[score,score,]
            #exp data
            d_exp = {}
            for ensg in exp_data_copy.columns:
                d_exp[ensg] = np.array(exp_data_copy.loc[:, ensg]).astype(float)  #ensg:[exp,exp,]
            #计算回归方程
            for ensg in d_cis:
                for pre in d_cis[ensg]:
                    #if pre in d_score:  #有些SNP两个型没有差异，没存到字典中
                    score = d_score[pre]
                    exp = d_exp[ensg]
                    beta, se, p, intercept, r_value = st_linear_reg(exp, score)
                    ff.write("%.5g" % p + "\n")
            print("Permutation", seed + 1, "finished")


if __name__ == "__main__":
    eQTac_permutation(
        "../Utilities_pipeline/example.output_eQTac/test.geno.vcf.gz.PRE_scoree",
        "../Utilities_pipeline/test_data/test.exp_residual", 100,
        "../Utilities_pipeline/example.output_eQTac/test.geno.vcf.gz.PRE_score.eQTac_result.permutation_plist")
