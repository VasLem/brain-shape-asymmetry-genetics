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

echo "Processing partition" $partition

echo "Processing GREAT"
path=../results/GREAT/$dataset/par$partition.txt
great=`cat $path | awk 'NR>1{print $1}'`

echo "Processing FUMA"
path=../results/FUMA/$dataset/par$partition/genes.txt
fuma=$(cat $path | awk 'NR>1{print $2}')
mkdir -p ../results/genes/$dataset
out=../results/genes/$dataset/par$partition.txt
echo "Extracting unique genes list to $out"
printf "$great\n$fuma"|sort -n|uniq > $out