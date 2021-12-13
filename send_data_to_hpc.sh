#!/bin/bash

rsync -av -zvr -LK --progress --files-from=sample_data_paths  .  hpc:/data/leuven/338/vsc33862/Thesis
