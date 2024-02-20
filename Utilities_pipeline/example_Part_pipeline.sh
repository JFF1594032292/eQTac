set -eo pipefail
#########################################################################
# File Name: example_Part_pipeline.sh
# Created on: 2022-12-26 13:47:05
# Author: JFF
# Last Modified: 2024-02-20 14:21:52
# Description: 
# Usage: 
# Input: 
# Output: 
#########################################################################
python Part-1-Train_model.py \
	-p test_data/test.positive.bed \
	-ex test_data/test.exclude.bed \
	-o output_eQTac_part \
	-t 3 -l 10 -k 6 -c 10 -g 2 -e 0.01 -f 2.5

python Part-2-Generate_PRE_fa.py \
	-pre test_data/test.pre.bed \
	--geno test_data/test.geno \
	--snp test_data/test.geno.snplist \
	-fa test_data/test.hg19.chr17.fa \
	-o output_eQTac_part

python Part-3-Predict_PRE_score.py \
	-m output_eQTac_part/test.positive.pos.svmmodel.3_10_6_0.01.model.txt \
	-l output_eQTac_part/test.geno.snplist.bed--test.pre.bed.pre_snplist.ld_info \
	-mfa output_eQTac_part/test.geno.snplist.bed--test.pre.bed.pre_snplist.ld_info.snplist.bed.mutate.fa \
	-geno test_data/test.geno \
	-snp output_eQTac_part/test.geno.snplist.bed--test.pre.bed.pre_snplist \
	-T 1 \
	-o output_eQTac_part

python Part-4-Calculate_eQTac_correlation.py \
	-pre output_eQTac_part/test.geno.vcf.gz.PRE_score \
	-exp test_data/test.exp_residual \
	-n 20 \
	-o output_eQTac_part
