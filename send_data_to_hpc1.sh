#!/bin/bash

rsync -av -zvr -LK --progress --files-from=sample_data_paths1  .  hpc:/data/leuven/338/vsc33862/Thesis
rsync -zvrLK --progress --files-from=code_data_paths  .  hpc:/user/leuven/338/vsc33862/imagen
