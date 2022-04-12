#!/bin/bash
set -e
source ../get_input_args.sh $1 $2 $3 $4
cd ../../python


ROOT_DIR=".."

MUNGED_DIR=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/munged
RG_DIR=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/rg
mkdir -p $RG_DIR
for i in {1..30}; do
    par_i=$(printf "%02.f" $i)
    echo Handling Partition $par_i
    files="$MUNGED_DIR/par$par_i.sumstats.gz"
    for ((j = $i + 1; j <= 31; j++)); do
        par_j=$(printf "%02.f" $j)
        files="$files,$MUNGED_DIR/par$par_j.sumstats.gz"
    done
    ./ldsc/ldsc.py --rg $files --ref-ld-chr $ROOT_DIR/SAMPLE_DATA/eur_w_ld_chr/ --w-ld-chr $ROOT_DIR/SAMPLE_DATA/eur_w_ld_chr/ --out $RG_DIR/par$par_i
done

ret="$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/rg/auto_correlation.csv"
rm -f $ret
for i in {1..30}; do
    par_i=$(printf "%02.f" $i)
    if [[ $i = 1 ]]; then
        first_line=2
    else
        first_line=3
    fi
    cat $RG_DIR/par$par_i.log | sed -n '/^Summary of Genetic/,/^Analysis finished/{p;/^Analysis finished/q}' | tail -n +$first_line | head -n -2 >>$ret
done