#!/bin/bash
for i in {01..31}; do
echo Handling Partition $i
./ldsc/munge_sumstats.py --sumstats ../results/genomeDemo/STAGE00DATA/CCAPart$i.csv \
 --n-min 19644 --snp "rsID" --p P-value --signed-sumstats ChiScore,1 --no-alleles  --delim , --out v1/MUNGE_OUTPUT$i

done