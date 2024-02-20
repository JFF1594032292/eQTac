set -eo pipefail
#########################################################################
# File Name: example_All_pipeline.sh
# Created on: 2022-12-19 22:22:33
# Author: JFF
# Last Modified: 2024-02-20 15:03:11
# Description: 
# Usage: 
# Input: 
# Output: 
#########################################################################
python Part-All-eQTac_pipeline.py \
	-p test_data/test.positive.bed \
	-ex test_data/test.exclude.bed \
	-pre test_data/test.pre.bed \
	--geno test_data/test.geno \
	--snp test_data/test.geno.snplist \
	-fa test_data/test.hg19.chr17.fa \
	-exp test_data/test.exp_residual \
	-n 100 \
	-o output_eQTac \
	-t 3 -l 10 -k 6 -c 10 -g 2 -e 0.01 
