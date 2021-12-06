#!/bin/bash
### PURPOSE OF THE SCRIPT: extract CpG
#prerequisites
module load bamtools/2.4.1

#Global variables
PATHOUTPUT="/PATH/TO/06_results/output_extractor/"
GENOMEFOLDER="/PATH/TO/REFERENCE_GENOME_FOLDER/04_reference"
PATHTOFILE="/PATH/TO/MAPPING_FOLDER/05_mapped/bismark_mapping"
NCPU=10

cd PATH/TO/MAPPING_FOLDER/05_mapped/bismark_mapping

ls *.bam > liste_temp
sed -r 's/_aR1_bismark_bt2_pe.bam//g' liste_temp > liste

conda activate

for line in $(cat liste)
do

#Methylation calls
bismark_methylation_extractor --gzip -p --bedGraph --scaffolds --cytosine_report --genome_folder "$GENOMEFOLDER" --multicore "$NCPU" -o "$PATHOUTPUT" "$PATHTOFILE"/"$line"_aR1_bismark_bt2_pe.bam

done 2>&1 | tee /PATH/TO/LOG_FILES/98_log_files/bismark_meth_extractor.log
