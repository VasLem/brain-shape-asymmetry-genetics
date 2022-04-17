#!/bin/bash
set -e
ROOT_DIR=".."
source ../get_input_args.sh 1 $2 $3 $4
MUNGED_DIR_DIS=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/munged

source ../get_input_args.sh 2 $2 $3 $4
MUNGED_DIR_REP=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/munged


cd ../../python


RG_DIR=$ROOT_DIR/results/$MODALITY/ldsc/disVsRep/rg
mkdir -p $RG_DIR
for i in {1..31}; do
    par_i=$(printf "%02.f" $i)
    echo Handling Partition $par_i
    
    dis_file="$MUNGED_DIR_DIS/par$par_i.sumstats.gz"
    rep_file="$MUNGED_DIR_REP/par$par_i.sumstats.gz"
    ./ldsc/ldsc.py --rg $dis_file,$rep_file --ref-ld-chr $ROOT_DIR/SAMPLE_DATA/eur_w_ld_chr/ --w-ld-chr $ROOT_DIR/SAMPLE_DATA/eur_w_ld_chr/ --out $RG_DIR/par$par_i
done

ret="$ROOT_DIR/results/$MODALITY/ldsc/disVsRep/rg/discovery_vs_replication_correlation.csv"
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
