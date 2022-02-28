#!/bin/bash
set -e
cd ../python
if [[ $1 == 1 ]]; then
    dataset=STAGE00DATA
elif [[ $1 == 2 ]]; then
    dataset=BATCH2_2021_DATA
else
    echo "Need to supply dataset index, 1 (for STAGE00DATA) or 2 (for BATCH2_2021_DATA)"
    exit
fi

mkdir -p ../results/ldsc/$dataset/munged/
mkdir -p ../results/ldsc/$dataset/h2/

for i in {1..31}; do
    i=$(printf "%02.f" $i)
    echo Handling Partition $i
    gunzip <../results/genomeDemo/$dataset/CCAPart$i.csv.gz >data.csv
    ./ldsc/munge_sumstats.py --sumstats data.csv \
        --n-min 1000 --snp "rsID" --p P-value --a1 "A1" --a2 "A2" --signed-sumstats ChiScore,1 --delim , --out ../results/ldsc/$dataset/munged/par$i
    ./ldsc/ldsc.py --h2 ../results/ldsc/$dataset/munged/par$i.sumstats.gz --ref-ld-chr ../SAMPLE_DATA/eur_w_ld_chr/ --w-ld-chr ../SAMPLE_DATA/eur_w_ld_chr/ --out ../results/ldsc/$dataset/h2/par$i
done

echo "Collecting stats"
for i in {1..31}; do
    i=$(printf "%02.f" $i)
    if [[ "$i" = "01" ]]; then
        cat ../results/ldsc/$dataset/h2/par$i.log | sed -nr 's/^(.*):\s+([0-9\.]+).*$/\1\t\2/p' | head -n4 | sed "1s/^/\t$i\n/" >ret.csv
    else
        paste ret.csv <(cat ../results/ldsc/$dataset/h2/par$i.log | sed -nr 's/^.*:\s+([0-9\.]+).*$/\1/p' | head -n4 | sed "1s/^/$i\n/") -d '\t' >temp && mv temp ret.csv
    fi
done
mv ret.csv ../results/ldsc/$dataset/h2_results.csv
cd -
