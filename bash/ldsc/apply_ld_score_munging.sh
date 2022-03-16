#!/bin/bash
set -e
source ../get_input_args.sh $1 $2 $3
cd ../../python
ROOT_DIR=".."
MUNGED_DIR=$ROOT_DIR/results/ldsc/$DATASET_ID/munged
H2_DIR=$ROOT_DIR/results/ldsc/$DATASET_ID/h2
mkdir -p $MUNGED_DIR
mkdir -p $H2_DIR

for i in {1..31}; do
    i=$(printf "%02.f" $i)
    echo Handling Partition $i
    gunzip <$ROOT_DIR/results/meta_analysis/$DATASET_ID/CCAPart$i.csv.gz >data.csv
    ./ldsc/munge_sumstats.py --sumstats data.csv \
        --n-min 1000 --snp "rsID" --p P-value --a1 "A1" --a2 "A2" \
        --signed-sumstats ChiScore,1 --delim , --out $MUNGED_DIR/par$i
done
