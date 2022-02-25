# STEPS
## Go to my directory
```
cd /usr/local/micapollo01/MIC/DATA/STUDENTS/vlemon0/code/deployed
```
## Send files required.
The syntax is ./send_data_to_hpc.sh DATASET (1 for discovery or 2 for replication dataset) SSH_URL _DIR. 
The following is an example for my case, I am user vsc33862, that saves the dataset 1 to
a folder named vlThesis (keep this folder as is, so that the pbs script does not need to change), 
using hpc as SSH_URL (I have made an SSH hook).
The allocation to VSC_DATA or VSC_HOME is done automatically, all the large files are saved in VSC_DATA.
> :warning: **Files of both datasets cannot exist at the same time in the data folder due to the 75GB limitation!**
```
./send_data_to_hpc.sh 1 hpc vlThesis
```
## Switch to remote, change the following commands as needed.
```
ssh hpc
cd vlThesis
```
## To run an array job for all the chromosomes:
```
./run_genome.sh
```
## To run only chromosome 20, a small one:
```
./run_genome_chr20.sh
```

