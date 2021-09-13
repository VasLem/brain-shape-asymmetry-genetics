#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=$PATH:$SCRIPT_DIR
echo hi
dir=results/genomics/$1
mkdir -p $dir
N=5
i=0
echo $N
for c in {1..2976}
do
echo $c
((i=i%N)); ((i++==0)) && wait
gcta64 --reml --grm "$dir" --pheno $2 --mpheno $c --out "$dir/pheno$c"  --thread-num $3 &
done