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
if [ "$#" -gt 1 ]; then
partition=$2
else
partition=1
fi
partition=$(printf "%02.f" $partition)

echo "Processing dataset " $dataset " partition" $partition
path=../results/FUMA/$dataset/par$partition/GenomicRiskLoci.txt
out=../results/GREAT/$dataset/genes_par$partition.txt
mkdir -p ../results/GREAT/$dataset
awk 'NR>1{printf ("chr%s %s %s %s\n",$4, $7, $8, $3)}' $path > $out