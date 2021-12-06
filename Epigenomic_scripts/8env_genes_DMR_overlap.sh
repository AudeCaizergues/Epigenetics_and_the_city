#!/bin/bash
### PURPOSE OF THE SCRIPT: find with which genes the DMR overlap

### Copy DMR files
#mkdir /PATH/TO/06_results/genes_DMR/
cp /PATH/TO/06_results/methylkit/methylkit_bar/dmr_tiles_bar_allow_missing_10pourcent_env.tsv /PATH/TO/06_results/genes_DMR/
cp /PATH/TO/06_results/methylkit/methylkit_mtp/dmr_tiles_mtp_allow_missing_10pourcent_env.tsv /PATH/TO/06_results/genes_DMR/
cp /PATH/TO/06_results/methylkit/methylkit_var/dmr_tiles_var_allow_missing_10pourcent_env.tsv /PATH/TO/06_results/genes_DMR/
cd /PATH/TO//06_results/genes_DMR/

awk 'NR > 1 { print }' < dmr_tiles_bar_allow_missing_10pourcent_env.tsv > DMR_bar_prov.txt
awk '{$1=""}1' DMR_bar_prov.txt | awk '{$1=$1}1' > DMR_bar_2prov.txt
awk -v OFS="\t" '$1=$1' DMR_bar_2prov.txt > DMR_bar.txt
awk 'NR > 1 { print }' < dmr_tiles_mtp_allow_missing_10pourcent_env.tsv > DMR_mtp_prov.txt
awk '{$1=""}1' DMR_mtp_prov.txt | awk '{$1=$1}1' > DMR_mtp_2prov.txt
awk -v OFS="\t" '$1=$1' DMR_mtp_2prov.txt > DMR_mtp.txt
awk 'NR > 1 { print }' < dmr_tiles_var_allow_missing_10pourcent_env.tsv > DMR_var_prov.txt
awk '{$1=""}1' DMR_var_prov.txt | awk '{$1=$1}1' > DMR_var_2prov.txt
awk -v OFS="\t" '$1=$1' DMR_var_2prov.txt > DMR_var.txt
rm *prov.txt

### find genes overlapping with DMR
bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_var.txt > genes_5kb_DMR_var.txt
bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_bar.txt > genes_5kb_DMR_bar.txt
bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_mtp.txt > genes_5kb_DMR_mtp.txt

awk '{print $4}' genes_5kb_DMR_bar.txt | sort | uniq > liste_genes_5kb_DMR_bar.txt
awk '{print $4}' genes_5kb_DMR_mtp.txt | sort | uniq > liste_genes_5kb_DMR_mtp.txt
awk '{print $4}' genes_5kb_DMR_var.txt | sort | uniq > liste_genes_5kb_DMR_var.txt
cat liste_genes_5kb_DMR_bar.txt liste_genes_5kb_DMR_mtp.txt liste_genes_5kb_DMR_var.txt | sort | uniq > liste_genes_5kb_DMR_toutesvilles.txt

### DMR intersections between cities

awk 'NR > 1 { print }' < matrix_DMR_10percent_intersection2ou3villes.txt > DMR_2ou3villes_prov.txt
awk '{print $1"\t"$2"\t"$3}' DMR_2ou3villes_prov.txt > DMR_2ou3villes.txt

awk 'NR > 1 { print }' < matrix_DMR_10percent_intersection2villes.tsv > DMR_2villes_prov.txt
awk '{print $1"\t"$2"\t"$3}' DMR_2villes_prov.txt > DMR_2villes.txt

awk 'NR > 1 { print }' < matrix_DMR_10percent_intersection3villes.tsv > DMR_3villes_prov.txt
awk '{print $1"\t"$2"\t"$3}' DMR_3villes_2prov.txt > DMR_3villes.txt

bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_2ou3villes.txt > genes_5kb_DMR_2ou3villes.txt
bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_2villes.txt > genes_5kb_DMR_2villes.txt
bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_3villes.txt > genes_5kb_DMR_3villes.txt

awk '{print $4}' genes_5kb_DMR_2ou3villes.txt | sort | uniq > liste_genes_5kb_DMR_2ou3villes.txt
awk '{print $4}' genes_5kb_DMR_2villes.txt | sort | uniq > liste_genes_5kb_DMR_2villes.txt
awk '{print $4}' genes_5kb_DMR_3villes.txt | sort | uniq > liste_genes_5kb_DMR_3villes.txt




