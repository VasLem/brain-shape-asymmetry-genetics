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

mkdir -p ../results/ldsc/$dataset/rg/
for i in {1..30}; do
    par_i=$(printf "%02.f" $i)
    echo Handling Partition $par_i
    files="../results/ldsc/$dataset/munged/par$par_i.sumstats.gz"
    for ((j = $i + 1; j <= 31; j++)); do
        par_j=$(printf "%02.f" $j)
        files="$files,../results/ldsc/$dataset/munged/par$par_j.sumstats.gz"
    done
    ./ldsc/ldsc.py --rg $files --ref-ld-chr ../SAMPLE_DATA/eur_w_ld_chr/ --w-ld-chr ../SAMPLE_DATA/eur_w_ld_chr/ --out ../results/ldsc/$dataset/rg/par$par_i
done
cd -

ret="../results/ldsc/$dataset/correlation.csv"
rm $ret
for i in {1..30}; do
    par_i=$(printf "%02.f" $i)
    if [[ $i = 1 ]]; then
        first_line=2
    else
        first_line=3
    fi
    cat ../results/ldsc/$dataset/rg/par$par_i.log | sed -n '/^Summary of Genetic/,/^Analysis finished/{p;/^Analysis finished/q}' | tail -n +$first_line | head -n -2 >>$ret
done