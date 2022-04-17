#!/bin/bash
set -e
source ../get_input_args.sh $1 $2 $3 $4
cd ../../python
ROOT_DIR=".."
MUNGED_DIR=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/munged
H2_DIR=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/h2
mkdir -p $MUNGED_DIR
mkdir -p $H2_DIR

for i in {1..31}; do
    i=$(printf "%02.f" $i)
    echo Handling Partition $i
    gunzip <$ROOT_DIR/results/$MODALITY/meta_analysis/$DATASET_ID/CCAPart$i.csv.gz >$MUNGED_DIR/tmp.csv
    ./ldsc/munge_sumstats.py --sumstats $MUNGED_DIR/tmp.csv \
        --n-min 1000 --snp "rsID" --p P-value --a1 "A1" --a2 "A2" \
        --signed-sumstats ChiScore,1 --delim , --merge-alleles $ROOT_DIR/SAMPLE_DATA/w_hm3.snplist --out $MUNGED_DIR/par$i
    rm $MUNGED_DIR/tmp.csv
done
