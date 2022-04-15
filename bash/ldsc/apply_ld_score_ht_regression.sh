#!/bin/bash
source ../get_input_args.sh $1 $2 $3 $4
cd ../../SAMPLE_DATA
ROOT_DIR=".."
MUNGED_DIR=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/munged
H2_CTS_DIR=$ROOT_DIR/results/$MODALITY/ldsc/$DATASET_ID/h2-cts

mkdir -p $MUNGED_DIR
mkdir -p $H2_CTS_DIR

for cts_name in Cahoy Multi_tissue_chromatin GTEx_brain; do  #Multi_tissue_gene_expr
    H2_CTS_SET_DIR=$H2_CTS_DIR/$cts_name
    mkdir -p $H2_CTS_SET_DIR
    baseline=$ROOT_DIR/SAMPLE_DATA/1000G_EUR_Phase3_baseline
    weights=$ROOT_DIR/SAMPLE_DATA/weights_hm3_no_hla
    seg_ldscores=$ROOT_DIR/SAMPLE_DATA/${cts_name}_1000Gv3_ldscores

    if [[ ! -d $baseline ]]; then
        wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
        tar -xvzf 1000G_Phase3_baseline_ldscores.tgz
        mv 1000G_EUR_Phase3_baseline $baseline
        rm 1000G_Phase3_baseline_ldscores.tgz
    fi
    if [[ ! -d $weights ]]; then
        wget https://data.broadinstitute.org/alkesgroup/LDSCORE/weights_hm3_no_hla.tgz
        tar -xvzf weights_hm3_no_hla.tgz
        mv weights_hm3_no_hla $weights
        rm weights_hm3_no_hla.tgz
    fi

    if [[ ! -d $seg_ldscores ]]; then
        wget https://data.broadinstitute.org/alkesgroup/LDSCORE/LDSC_SEG_ldscores/${cts_name}_1000Gv3_ldscores.tgz
        tar -xvzf ${cts_name}_1000Gv3_ldscores.tgz
        rm ${cts_name}_1000Gv3_ldscores.tgz
    fi
    echo "Performing analysis"
    N=10
    for i in {1..31}; do
        ((j=j%N)); ((j++==0)) && wait
        i=$(printf "%02.f" $i)
        echo Handling Partition $i
        ../python/ldsc/ldsc.py --h2-cts $MUNGED_DIR/par$i.sumstats.gz \
            --ref-ld-chr $baseline/baseline. \
            --out $H2_CTS_SET_DIR/par$i\
            --ref-ld-chr-cts $cts_name.ldcts \
            --w-ld-chr $weights/weights. &
    done

    echo "Collecting stats"
    for i in {1..31}; do
        i=$(printf "%02.f" $i)
        if [[ "$i" = "01" ]]; then
            cat $H2_CTS_SET_DIR/par$i.log | sed -nr 's/^(.*):\s+([0-9\.]+).*$/\1\t\2/p' | head -n4 | sed "1s/^/\t$i\n/" >ret.csv
        else
            paste ret.csv <(cat $H2_CTS_SET_DIR/par$i.log | sed -nr 's/^.*:\s+([0-9\.]+).*$/\1/p' | head -n4 | sed "1s/^/$i\n/") -d '\t' >temp && mv temp ret.csv
        fi
    done
    mv ret.csv $H2_CTS_SET_DIR/h2_results.csv

    ret=$H2_CTS_SET_DIR/h2_results.csv
    rm -f $ret
    for i in {1..31}; do
        par_i=$(printf "%02.f" $i)
        if [[ $i = 1 ]]; then
            head -n 1 $H2_CTS_SET_DIR/par$par_i.cell_type_results.txt | awk -F '\t' -v OFS='\t' '{ $(NF+1) = "partition"; print }' >> $ret
        fi
        tail -n +2 $H2_CTS_SET_DIR/par$par_i.cell_type_results.txt | awk -v VARIABLE=$i -F '\t' -v OFS='\t' '{ $(NF+1) = VARIABLE; print }' >>$ret
    done
done

cd -
