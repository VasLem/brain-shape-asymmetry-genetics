set -e 

source ../get_input_args.sh $1 $2 $3 $4
cd ../../python

if [ "$#" -gt 3 ]; then
partition=$5
else
partition=1
fi
partition=$(printf "%02.f" $partition)

echo "Processing partition" $partition

echo "Processing GREAT"
path=../results/$MODALITY/GREAT/$DATASET_ID/par$partition.txt
great=`cat $path | awk 'NR>1{print $1}'`

echo "Processing FUMA"
path=../results/$MODALITY/FUMA/$DATASET_ID/par$partition/genes.txt
fuma=$(cat $path | awk 'NR>1{print $2}')
mkdir -p ../results/$MODALITY/genes/$DATASET_ID
out=../results/$MODALITY/genes/$DATASET_ID/par$partition.txt
echo "Extracting unique genes list to $out"
printf "$great\n$fuma"|sort -n|uniq > $out