#!/bin/bash
if [ $# -eq 0 ]
  then
    echo "No arguments supplied"
fi
qsub run_hierarchical.pbs $1
