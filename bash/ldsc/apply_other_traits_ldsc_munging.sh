#!/bin/bash
set -e
cd ../../python
ROOT_DIR=".."
for trait in face brain_shape; do
    if [[ $trait -eq 'face' ]]; then
        N_PARTS=63
        MEDIAN=74.8
    else
        N_PARTS=285
        MEDIAN=285
    fi
    MUNGED_DIR=$ROOT_DIR/results/ldsc/$trait/munged
    H2_DIR=$ROOT_DIR/results/ldsc/$trait/h2
    mkdir -p $MUNGED_DIR
    mkdir -p $H2_DIR
    i=1;
    while [ $i -le $N_PARTS ]; do
        par_i=$(printf "%02.f" $i)
        echo Handling Partition $i
        gunzip <$ROOT_DIR/results/other_traits_gwas/$trait/CCAPart$par_i.csv.gz >data.csv
        ./ldsc/munge_sumstats.py --sumstats data.csv \
            --n-min 1000 --snp "rsID" --p P-value --a1 "A1" --a2 "A2" \
            --signed-sumstats ChiScore,$MEDIAN --delim , --out $MUNGED_DIR/par$par_i
        ./ldsc/ldsc.py --h2 $MUNGED_DIR/par$par_i.sumstats.gz \
            --ref-ld-chr $ROOT_DIR/SAMPLE_DATA/eur_w_ld_chr/ \
            --w-ld-chr $ROOT_DIR/SAMPLE_DATA/eur_w_ld_chr/ --out $H2_DIR/par$par_i
        i=$(( $i + 1 ))
    done
done
