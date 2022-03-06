#!/bin/bash
set -e

echo Filtering european population
cd ../../SAMPLE_DATA/1000GenomeProject
mkdir -p euro
for i in {1..22}; do
echo Handling chromosome $i
bcftools-1.15/build/bin/bcftools view --force-samples -Ov -S  eu_ids chr$i.1kg.phase3.v5a.vcf.gz | java -jar bref3.08Feb22.fa4.jar > euro/chr$i.1kg.phase3.v5a.b37.euro.bref3
done