#!/bin/bash
### PURPOSE OF THE SCRIPT: find with which genes the sex DMR overlap



### Copier les fichiers DMR
#mkdir /PATH/TO/06_results/genes_DMR/
cd /PATH/TO/06_results/sex_DMR/genes_sex

awk 'NR > 1 { print }' < dmr_tiles_bar_allow_missing_10pourcent_sex.tsv > DMR_bar_prov.txt
awk '{$1=""}1' DMR_bar_prov.txt | awk '{$1=$1}1' > DMR_bar_2prov.txt
awk -v OFS="\t" '$1=$1' DMR_bar_2prov.txt > DMR_bar_sex.txt
awk 'NR > 1 { print }' < dmr_tiles_mtp_allow_missing_10pourcent_sex.tsv > DMR_mtp_prov.txt
awk '{$1=""}1' DMR_mtp_prov.txt | awk '{$1=$1}1' > DMR_mtp_2prov.txt
awk -v OFS="\t" '$1=$1' DMR_mtp_2prov.txt > DMR_mtp_sex.txt
awk 'NR > 1 { print }' < dmr_tiles_var_allow_missing_10pourcent_sex.tsv > DMR_var_prov.txt
awk '{$1=""}1' DMR_var_prov.txt | awk '{$1=$1}1' > DMR_var_2prov.txt
awk -v OFS="\t" '$1=$1' DMR_var_2prov.txt > DMR_var_sex.txt
rm *prov.txt

### find genes overlapping with DMR
bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_var_sex.txt > genes_5kb_DMR_var_sex.txt
bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_bar_sex.txt > genes_5kb_DMR_bar_sex.txt
bedtools window -bed -w 5000 -a genes_pmajor_ok.bed -b DMR_mtp_sex.txt > genes_5kb_DMR_mtp_sex.txt

awk '{print $4}' genes_5kb_DMR_bar_sex.txt | sort | uniq > liste_genes_5kb_DMR_bar_sex.txt
awk '{print $4}' genes_5kb_DMR_mtp_sex.txt | sort | uniq > liste_genes_5kb_DMR_mtp_sex.txt
awk '{print $4}' genes_5kb_DMR_var_sex.txt | sort | uniq > liste_genes_5kb_DMR_var_sex.txt


cat liste_genes_5kb_DMR_bar_sex.txt liste_genes_5kb_DMR_mtp_sex.txt liste_genes_5kb_DMR_var_sex.txt | sort | uniq > liste_genes_5kb_DMR_sex_toutesvilles.txt


### Building a list of genes for circos plot
awk  '{print $1 "        " $2 "        " $4}' genes_10kb_DMR_bar.txt | sort | uniq > liste_genes_10kb_DMR_bar_pour_circos.txt
awk  '{print $1 "        " $2 "        " $4}' genes_10kb_DMR_mtp.txt | sort | uniq > liste_genes_10kb_DMR_mtp_pour_circos.txt
awk  '{print $1 "        " $2 "        " $4}' genes_10kb_DMR_var.txt | sort | uniq > liste_genes_10kb_DMR_var_pour_circos.txt

cat liste_genes_10kb_DMR_bar_pour_circos.txt liste_genes_10kb_DMR_mtp_pour_circos.txt liste_genes_10kb_DMR_var_pour_circos.txt | sort | uniq > liste_genes_10kb_DMR_toutesvilles_pour_circos.txt
