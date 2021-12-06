#!/bin/bash
### PURPOSE OF THE SCRIPT: trim data with fastp

INPUT="/PATH/TO/DATA_FOLDER/02_data"
OUTPUT="/PATH/TO/TRIMMED_DATA_FOLDER/03_trimmed"
LENGTH=40
QUAL=15

for line in $(cat list_of_individuals)
do
./fastp --thread 12 --compression 1 -i "$INPUT"/"$line"_R1_001.fastq.gz -I "$INPUT"/"$line"_R2_001.fastq.gz \
	--length_required="$LENGTH" --qualified_quality_phred="$QUAL" \
	--correction \
	--trim_tail1=1 --trim_tail2=1 \
	-o "$OUTPUT"/"$line"_R1_001.fastq.gz -O "$OUTPUT"/"$line"_R2_001.fastq.gz
done 2>&1 | tee /PATH/TO/LOG_FILES/98_log_files/fastp_trim.log
