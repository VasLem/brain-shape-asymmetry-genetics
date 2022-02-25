#!/bin/bash
set -e
{
module purge &&
module load matlab/R2019a &&
LOCAL=0
} || {
echo Assuming Local Deployment
MCC_PATH=/SOFTWARE/Matlab/R2019a/bin/
LOCAL=1
}

mkdir -p deployed
cd matlab
${MCC_PATH}mcc -mv -R -singleCompThread -R -nodisplay  demoBrainAsymmetryGenome.m -d ../deployed/ \
 -a ./AIDFUNCTIONS -a ./FUNCTIONS -a ./BrainAsymmetrySignificanceAnalysis \
 -a ../SNPLIB
cd ..
if [[ $LOCAL ]]; then
rsync -av -zvr -LK --progress deployed/demoBrainAsymmetryGenome hpc:/user/leuven/338/vsc33862/vlThesis
fi

