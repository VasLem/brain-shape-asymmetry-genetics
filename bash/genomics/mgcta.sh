#!/bin/bash
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
export PATH=$PATH:$SCRIPT_DIR
dir=dirname $1
N=5
i=0
numOfParams=head -1 $1 | sed 's/[^,]//g' | wc -c
numOfParams-=2
for c in {1..$numOfParams}
do
((i=i%N)); ((i++==0)) && wait
gcta64 --reml --grm "$dir" --pheno $2 --mpheno $c --out "$dir/pheno$c"  --thread-num $3 &
done