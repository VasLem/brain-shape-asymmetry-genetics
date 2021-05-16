#!/bin/bash
# excluding SNPs with minor allele frequency <1%, Hardy–Weinberg equilibrium test P value <1e−7 and imputation INFO score <0.7
# keeping only biallelic SNPs
set -e
export PATH=$PATH:.


plink1 --noweb --bfile $1 --maf 0.01 --hwe 0.0000001 --allow-no-sex --proxy-impute-threshold 0.7   --make-bed --out "results/genomics/processed_`basename $1`" --pheno $2

gcta64 --bfile "results/genomics/processed_`basename $1`"  --thread-num $3 --make-grm --out processed_grm
