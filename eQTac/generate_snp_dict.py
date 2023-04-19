# -*- coding: utf-8 -*-
#########################################################################
# File Name: generate_snp_dict.py
# Created on : 2022-12-21 15:43:12
# Author: JFF
# Last Modified: 2022-12-21 15:43:14
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re


def generate_snp_dict(bim_file):
    """
    Generate a dict from a plink .bim file. Used in generate_PRE and generate_mut_fa.
    
    ## Parameters
    bim_file: bim file path
        chr\tsnp\tcmol\tpos\tA1\tA2\n
        chr1\\trs202152658\\t0\\t751343\\tA\\tT
    
    ## Returns
    d_snp: dict
        {snp:[chr,pos,a1,a2],...}
    """
    d_snp = {}
    with open(bim_file, 'r') as f1:
        for i in f1:
            i = i.strip().split("\t")
            d_snp[i[1]] = [i[0], i[3], i[4], i[5]]  #snp:[chr,pos,a1,a2]
    return d_snp
