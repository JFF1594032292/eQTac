# eQTac
---
EQTac is a method to predict the potential regulatory elements (PREs) and their target genes, based on the eQTL datasets, The only additional data was ATAC-seq or ChIP-seq peak data. 
## Dependence
### Python packages
```
numpy >= 1.22.4
pandas >= 1.4.3
pybedtools >= 0.8.2
pysam >= 0.16.0.1
rpy2 >= 3.4.2
scipy >= 1.8.1
```
### Other software
```
plink >= v1.90b6.24
bedtools >= v2.30.0
R >= 3.6.1
    r-gkmSVM >= 0.8.0
```
## Installation

## Input data
1. Data used in model training:
    (1) **Positive sets in bed format.** It's usually the peak data from ATAC-seq or ChIP-seq, we recomended to trim peaks to the core region (e.g. summits $\pm$ 100bp). See test_data/test.positive.bed
    (2) **Excluded sets in bed format.** It's usually the peak data from ATAC-seq or ChIP-seq, but with more relaxed thresholds (e.g. p=0.2). These region will be removed from genrated negative regions, in order to remove potential positive sequences from negative sets. See test_data/test.exclude.bed.
    (3) **Fasta file with .fai index.** Usually the human genome sequnce file in fasta format. See test.hg19.chr17.fa.
2. Data used in eQTac calculation.
    (1) 
## Usage pattern
We provided three level patterns: (1) pipeline level. (2) part level. (3) function level.
### Pipeline level pattern
For the function level pattern, we provide a script: Part-All-eQTac_pipeline.py.
It can be used as follow:
```
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
```
### Part level pattern
For the function level pattern, we provide four scripts:
```
python Part-1-Train_model.py \
	-p test_data/test.positive.bed \
	-ex test_data/test.exclude.bed \
	-o output_eQTac_part \
	-t 3 -l 10 -k 6 -c 10 -g 2 -e 0.01

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
	-n 50 \
	-o output_eQTac_part
```
### Function level pattern
For the function level pattern, we provide a series of functions:
```
from eQTac.get_nullseq import get_nullseq
from eQTac.filter_bkg import filter_bkg
from eQTac.generate_snp_dict import generate_snp_dict
from eQTac.generate_PRE import generate_PRE
from eQTac.generate_mut_fa import generate_mut_fa
from eQTac.geno2score import geno2score
from eQTac.eQTac_correlation import eQTac_correlation
from eQTac.eQTac_permutation import eQTac_permutation
from eQTac.control_FDR import control_FDR
```
These functions can be used to construct the whole pipeline. The Part-All-eQTac_pipeline.py file used this

## Notes