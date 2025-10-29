#!/bin/bash

python3 /scratch/project_2009813/JOHANNA/MAFFT/999scripts/batch_conservation_v2.py \
  --threshold 0.6 --min-length 10 --merge-distance 10 \
  --output-dir conservation_results_60pct_10aa_merge10