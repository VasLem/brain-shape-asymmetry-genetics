if [[ $1 == 1 ]]; then
    dataset=STAGE00DATA
elif [[ $1 == 2 ]]; then
    dataset=BATCH2_2021_DATA
else
    echo "Need to supply dataset index, 1 (for STAGE00DATA) or 2 (for BATCH2_2021_DATA)"
    exit
fi

impute_id=$2
if [[ -z $impute_id ]]; then
impute_id='median_imputed'
fi

subsampled_id=$3
if [[ -z $subsampled_id ]]; then
subsampled_id='not_subsampled'
fi

echo "DATASET:$dataset"
echo "IMPUTE_ID:$impute_id"
echo "SUBSAMPLED_ID:$subsampled_id"
DATASET_ID=$dataset/$impute_id$/$subsampled_id