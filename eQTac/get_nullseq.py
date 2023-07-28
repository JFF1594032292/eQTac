# -*- coding: utf-8 -*-
#########################################################################
# File Name: get_nullseq.py
# Created on : 2022-12-17 16:06:59
# Author: JFF
<<<<<<< HEAD
# Last Modified: 2023-07-28 10:24:14
=======
# Last Modified: 2022-12-17 16:07:01
>>>>>>> 4b8bfa95564736c3bf45c48056ea656cf880d680
# Description:
# Usage:
# Input:
# Output:
#########################################################################
import os
import sys
import re
import rpy2.robjects as robjects
from shutil import copyfile


<<<<<<< HEAD
def get_nullseq(pos_bed_path, out_folder, xfold):
    """
    Generates null sequences (negative set) with matching repeat and GC content as the input bed file for positive set regions.\n
    A wrapper of genNullSeqs from gkmSVM R package.

=======
def get_nullseq(pos_bed_path, out_folder):
    """
    Generates null sequences (negative set) with matching repeat and GC content as the input bed file for positive set regions.\n
    A wrapper of genNullSeqs from gkmSVM R package.
    
>>>>>>> 4b8bfa95564736c3bf45c48056ea656cf880d680
    ## Parameters
    pos_bed_path
        positive training set bed file. e.g. ATAC peak bed file.
    out_folder
        output directory
<<<<<<< HEAD

=======
    
>>>>>>> 4b8bfa95564736c3bf45c48056ea656cf880d680
    ## Returns
    pos_fa
        positive sequence .fa file
    pos_bed_newpath
        positive sequence .bed file. (copy to output directory)
    negraw_fa
        raw negative sequence .fa file.
    negraw_bed
        raw negative sequence .bed file.
<<<<<<< HEAD

    ## Notes
    The output ``negraw_fa`` and ``negraw_bed`` file not remove the duplicate sequences, and raw negative sequences may also contain potential positive sequences.(e.g. peaks with P<0.05 but q>0.05)\n
    The ``negraw_fa`` and ``negraw_bed`` showed filtered by ``filter_bkg`` function.

    """
    pos_bed = os.path.basename(pos_bed_path)
    prefix = re.sub("\.bed$", "", pos_bed)
    # make dir
=======
    
    ## Notes
    The output ``negraw_fa`` and ``negraw_bed`` file not remove the duplicate sequences, and raw negative sequences may also contain potential positive sequences.(e.g. peaks with P<0.05 but q>0.05)\n
    The ``negraw_fa`` and ``negraw_bed`` showed filtered by ``filter_bkg`` function.
    
    """
    pos_bed = os.path.basename(pos_bed_path)
    prefix = re.sub("\.bed$", "", pos_bed)
    #make dir
>>>>>>> 4b8bfa95564736c3bf45c48056ea656cf880d680
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)
    # copy pos_bed to output folder & add a "pos" suffix
    pos_bed_newpath = out_folder + "/" + prefix + ".pos.bed"
    copyfile(pos_bed_path, pos_bed_newpath)
    # gkmSVM genNullSeqs
    pos_fa = out_folder + "/" + prefix + ".pos.fa"
    negraw_bed = out_folder + "/" + prefix + ".negraw.bed"
    negraw_fa = out_folder + "/" + prefix + ".negraw.fa"
<<<<<<< HEAD
    d = dict(pos_bed_newpath=pos_bed_newpath, pos_fa=pos_fa, negraw_bed=negraw_bed, negraw_fa=negraw_fa, xfold=xfold)
    robjects.r(
        "library(gkmSVM);genNullSeqs('%(pos_bed_newpath)s',nMaxTrials=20,xfold='%(xfold)f',GC_match_tol=0.05,repeat_match_tol=0.05,length_match_tol=0.05,batchsize=500000,genomeVersion='hg19', outputPosFastaFN='%(pos_fa)s', outputBedFN='%(negraw_bed)s', outputNegFastaFN='%(negraw_fa)s')"
=======
    d = dict(pos_bed_newpath=pos_bed_newpath, pos_fa=pos_fa, negraw_bed=negraw_bed, negraw_fa=negraw_fa)
    robjects.r(
        "library(gkmSVM);genNullSeqs('%(pos_bed_newpath)s',nMaxTrials=20,xfold=2.5,GC_match_tol=0.05,repeat_match_tol=0.05,length_match_tol=0.05,batchsize=500000,genomeVersion='hg19', outputPosFastaFN='%(pos_fa)s', outputBedFN='%(negraw_bed)s', outputNegFastaFN='%(negraw_fa)s')"
>>>>>>> 4b8bfa95564736c3bf45c48056ea656cf880d680
        % d)
    return pos_fa, pos_bed_newpath, negraw_fa, negraw_bed


if __name__ == '__main__':
<<<<<<< HEAD
    get_nullseq("../Utilities_pipeline/test_data/test.positive.bed", ".", xfold=2.5)
=======
    get_nullseq("../Utilities_pipeline/test_data/test.positive.bed", ".")
>>>>>>> 4b8bfa95564736c3bf45c48056ea656cf880d680
