#!/bin/bash
for i in {01..31}; do
i=$(printf "%02d" $i)
echo Handling Partition $i
gunzip < ../results/genomeDemo/STAGE00DATA/CCAPart$i.csv.gz > data.csv
./ldsc/munge_sumstats.py --sumstats data.csv \
 --n-min 19644 --snp "rsID" --p P-value --a1 "A1" --a2 "A2" --signed-sumstats ChiScore,1  --delim , --out v1/MUNGE_OUTPUT$i
./ldsc/ldsc.py --h2 v1/MUNGE_OUTPUT$i.sumstats.gz --ref-ld-chr ../SAMPLE_DATA/eur_w_ld_chr/ --w-ld-chr ../SAMPLE_DATA/eur_w_ld_chr/ --out v2/MUNGE_OUTPUT$i
done

for i in {1..31}; do
i=$(printf "%02d" $i)
if [[ "$i" = "01" ]]; then
cat v2/MUNGE_OUTPUT$i.log | sed -nr  's/^(.*):\s+([0-9\.]+).*$/\1\t\2/p'| head -n4  | sed "1s/^/\t$i\n/">  ret.csv
else
paste ret.csv <(cat v2/MUNGE_OUTPUT$i.log | sed -nr  's/^.*:\s+([0-9\.]+).*$/\1/p'| head -n4 | sed "1s/^/$i\n/") -d '\t'  > temp && mv temp ret.csv
fi
done
mv ret.csv ../results/genomeDemo/STAGE00DATA/ldsc_heritability_results.csv