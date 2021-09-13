#!/bin/bash
# excluding SNPs with minor allele frequency <1%, Hardy–Weinberg equilibrium test P value <1e−7 and imputation INFO score <0.7
# keeping only biallelic SNPs
set -e
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=$PATH:$SCRIPT_DIR


plink1 --noweb --bfile $1 --maf 0.01 --hwe 0.0000001 --allow-no-sex --proxy-impute-threshold 0.7   --make-bed --out "results/genomics/processed_`basename $1`" --pheno $2

gcta64 --bfile "results/genomics/processed_`basename $1`"  --thread-num $3 --make-grm --out "results/genomics/processed_grm_`basename $1`"
gcta64 --grm "results/genomics/processed_grm_`basename $1`" --grm-cutoff 0.05 --make-grm --out "results/genomics/pruned_cutoff_05_grm_`basename $1`"
gcta64 --grm "results/genomics/pruned_cutoff_05_grm_`basename $1`" 


