#!/bin/bash
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
fi
module load worker/1.6.10-intel-2018a
qsub  -t 1-22  run_genome.pbs $1
