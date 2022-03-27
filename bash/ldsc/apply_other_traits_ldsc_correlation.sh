#!/bin/bash
set -e
source ../get_input_args.sh $1 $2 $3 $4
cd ../../python

ROOT_DIR=".."
MUNGED_DIR=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/munged
traits=$(ls $ROOT_DIR/SAMPLE_DATA/OTHER_TRAITS_GWAS)
RG_DIR=$ROOT_DIR/results/$MODALITY/ldsc/other_traits/$DATASET_ID/rg
for i in {1..31}; do
    par_i=$(printf "%02.f" $i)
    echo Handling Partition $par_i
    as_file="$MUNGED_DIR/par$par_i.sumstats.gz"
    tr_files=
    for trait in $traits; do
        mkdir -p $RG_DIR
        MUNGED_TRAIT_DIR=$ROOT_DIR/SAMPLE_DATA/OTHER_TRAITS_GWAS/$trait
        tr_files=$tr_files,$MUNGED_TRAIT_DIR/munged.sumstats.gz
    done
    ./ldsc/ldsc.py --rg $as_file$tr_files --ref-ld-chr $ROOT_DIR/SAMPLE_DATA/eur_w_ld_chr/ --w-ld-chr $ROOT_DIR/SAMPLE_DATA/eur_w_ld_chr/ --out $RG_DIR/par$par_i
done

mkdir -p $ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/rg/
ret="$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/rg/other_traits_correlation.csv"
rm -f $ret
for i in {1..31}; do
    par_i=$(printf "%02.f" $i)
    if [[ $i = 1 ]]; then
        first_line=2
    else
        first_line=3
    fi
    cat $RG_DIR/par$par_i.log | sed -n '/^Summary of Genetic/,/^Analysis finished/{p;/^Analysis finished/q}' | tail -n +$first_line | head -n -2 >>$ret
done

