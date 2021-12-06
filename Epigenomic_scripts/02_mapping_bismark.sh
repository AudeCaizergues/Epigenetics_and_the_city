#!/bin/bash
### PURPOSE OF THE SCRIPT: map RRBS reads on reference genome

#### Prepare reference ####
#variables
GENOME_FOLDER="/PATH/TO/REFERENCE_GENOME_FOLDER/04_reference"
#prepare genome
bismark_genome_preparation $GENOME_FOLDER  2>&1 | tee /PATH/TO/LOG_FILES/98_log_files/refindex_bismark.log


##### Mapping ####
#global variables
GENOME_FOLDER="/PATH/TO/REFERENCE_GENOME_FOLDER/04_reference"
DATAFOLDER="/PATH/TO/TRIMMED_DATA_FOLDER/03_trimmed"
DATAOUTPUT="/PATH/TO/MAPPING_FOLDER/05_mapped/bismark_mapping"
N="-N 1"			
L="-L 20"			
p="-p 4"			

conda activate

for base in $(cat list_of_individuals)
do
zcat "$DATAFOLDER"/"$base"_R1_001.fastq.gz >"$DATAFOLDER"/"$base"_aR1.fq
zcat "$DATAFOLDER"/"$base"_R2_001.fastq.gz >"$DATAFOLDER"/"$base"_aR2.fq
bismark $N $L $p -q $GENOME_FOLDER -1 "$DATAFOLDER"/"$base"_aR1.fq -2 "$DATAFOLDER"/"$base"_aR2.fq
rm "$DATAFOLDER"/"$base"_aR*.fq
done 2>&1 | tee /PATH/TO/LOG_FILES/98_log_files/mapping_bismark.log

conda deactivate 

cd PATH/TO/MAPPING_FOLDER/05_mapped/bismark_mapping

for line in $(cat list_of_individuals)
do
samtools sort "$line"_aR1_bismark_bt2_pe.bam > "$line"_sorted.bam
samtools index "$line"_sorted.bam 
done


