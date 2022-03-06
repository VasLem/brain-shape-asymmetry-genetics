set -e 
source ../get_input_args $1 $2 $3
cd ../../python

if [ "$#" -gt 3 ]; then
partition=$4
else
partition=1
fi
partition=$(printf "%02.f" $partition)

echo "Processing dataset " $DATASET_ID " partition" $partition
path=../results/FUMA/$DATASET_ID/par$partition/GenomicRiskLoci.txt
out=../results/GREAT/$DATASET_ID/genes_par$partition.txt
mkdir -p ../results/GREAT/$DATASET_ID
awk 'NR>1{printf ("chr%s %s %s %s\n",$4, $7, $8, $3)}' $path > $out