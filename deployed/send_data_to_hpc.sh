#!/bin/bash
set -e
DATASET=$1
SSH_CONNECTION=$2
SSH_FOLDER=$3

if [[ -z "$DATASET" ]]; then
DATASET=1
fi
if [ $DATASET -ne 1 -a $DATASET -ne 2 ]; then
    echo "Invalid value provided for dataset, only 1 or 2 accepted"
fi
if [[ -z "$SSH_CONNECTION" ]]; then
SSH_CONNECTION=hpc
fi
if [[ -z "$SSH_FOLDER" ]]; then
SSH_FOLDER=vlThesis
fi

USER_FOLDER=`ssh hpc "set" | grep -i HOME | sed -rn 's/HOME=(.*)/\1/p'`
DATA_FOLDER=`sed  's/user/data/g' <<< $USER_FOLDER`

cd ..
echo Sending data
cat deployed/sample_data_paths$DATASET
rsync -av -zvr -LK --progress --files-from="deployed/sample_data_paths$DATASET"  .  $SSH_CONNECTION:$DATA_FOLDER/$SSH_FOLDER
echo Sending code
cd -
rsync -zvrLK --progress *  $SSH_CONNECTION:$USER_FOLDER/$SSH_FOLDER

