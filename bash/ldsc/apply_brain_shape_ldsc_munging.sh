#!/bin/bash
set -e
cd ../../python
ROOT_DIR=".."
MODALITY=asymmetry
N_PARTS=285
MEDIAN=285
MUNGED_DIR=$ROOT_DIR/results/brain_shape/ldsc/munged
H2_DIR=$ROOT_DIR/results/brain_shape/ldsc/h2
mkdir -p $MUNGED_DIR
mkdir -p $H2_DIR
i=1;
while [ $i -le $N_PARTS ]; do
    par_i=$(printf "%02.f" $i)
    echo Handling Partition $i
    gunzip <$ROOT_DIR/SAMPLE_DATA/BRAIN_SHAPE_PARTITIONS/par$par_i.csv.gz >data.csv
    ./ldsc/munge_sumstats.py --sumstats data.csv \
        --n-min 1000 --snp "rsID" --p P-value --a1 "A1" --a2 "A2" \
        --signed-sumstats ChiScore,$MEDIAN --delim , --out $MUNGED_DIR/par$par_i
    i=$(( $i + 1 ))
done
