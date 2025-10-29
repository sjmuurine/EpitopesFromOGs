#!/bin/bash

cd /scratch/project_2009813/JOHANNA/MAFFT/Data/OrthgroupApproach/

for og in OG0000089 OG0000067 OG0000090 OG0000034 OG0000019 OG0000083 OG0000072 OG0000080; do
  echo '=== '$og' ==='
  grep -A 5 'Region' conservation_results_60pct_8aa_merge10/${og}_conservation.txt
  echo ''
done > orthogroup_regions.txt

echo "Output saved to orthogroup_regions.txt"
