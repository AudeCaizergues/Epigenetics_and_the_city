#!/bin/bash
### PURPOSE OF THE SCRIPT : from _GpG.txt files -> filter on coverage 
### (10X min) for RDA analyses


cd /PATH/TO/06_results/output_extractor/
ls *_CpG.txt > liste_CpG_files

mkdir /PATH/TO/06_results/rda

for line in $(cat liste_CpG_files)
do
awk '$5>=10' "$line" > /PATH/TO/06_results/rda/"$line"_filtered10cov.txt
done


