# -*- coding: utf-8 -*-
#########################################################################
# File Name: control_FDR.py
# Created on : 2022-12-20 19:51:25
# Author: JFF
# Last Modified: 2023-04-18 17:28:48
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
import time
import pandas as pd
import numpy as np


def control_FDR(PRE_scorefile, permute_file, output_file):
    permute_plist = [float(i.strip()) for i in open(permute_file)]
    # observe_data 处理
    observe_data = pd.read_table(PRE_scorefile, header=0, index_col=None, sep="\t")
    observe_data = observe_data.sort_values(by="p")
    observe_data.loc[:, "FDR"] = 1
    observe_plist = list(observe_data.loc[:, "p"])
    # permute plist处理
    n_fold = len(permute_plist) / len(observe_plist)  # 背景集是真实值数量的多少倍
    permute_plist = np.array(permute_plist)
    # 计算FDR及阈值
    fdr_list = []
    for i in range(len(observe_plist)):
        n_lower = np.sum(permute_plist < observe_plist[i]) / n_fold
        fdr = n_lower / (i + 1)
        fdr_list.append(fdr)
    # FDR p值排序在后面的，FDR应大于等于前面的FDR值
    fdr_last = -1
    for i in range(len(fdr_list)):
        fdr_list[i] = max(fdr_list[i], fdr_last)
        fdr_last = fdr_list[i]
    observe_data.loc[:, "FDR"] = [i if i < 1 else 1 for i in fdr_list]
    observe_data = observe_data.sort_values(by="p")
    # 保存文件
    observe_data.to_csv(output_file, sep="\t", header=True, index=False)


if __name__ == "__main__":
    control_FDR(
        "../Utilities_pipeline/example.output_eQTac/test.geno.vcf.gz.PRE_score.eQTac_result",
        "../Utilities_pipeline/example.output_eQTac/test.geno.vcf.gz.PRE_score.eQTac_result.permutation_plist",
        "../Utilities_pipeline/example.output_eQTac/test.geno.vcf.gz.PRE_score.eQTac_result.FDR.txt"
    )
