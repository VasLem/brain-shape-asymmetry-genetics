#!/bin/bash
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
fi
qsub -t 20 run_genome.pbs $1
