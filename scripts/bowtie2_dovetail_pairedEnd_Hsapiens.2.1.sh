#!/bin/bash
# script for mapping paired end ChIP-seq/ATAC-seq reads with bowtie2 when there is one pair of fastq files per sample
# variables include 
# $1 index- path to premade bowtie2 index and index suffix: egz. "mm39"
# $2 sample
# $3 outFolder
#
# output: sorted bam file with bai index in the output folder & stats file in the "trimmedFastq" folder
# for compatibility with UCSC browser only "standard" chromosomes are kept; reads aligned to mitochondrial DNA
# discarded as they are highly enriched in ATAC-seq due to preferential Tn5 integration into mtDNA
#!/bin/bash
#
# USAGE EXAMPLE:
# cd /home/micgdu/Analysis/Epidermis/geneExpression_Fan_Nayat/trimmedFastq/ATAC_seq_Nayat/
# samples=($(ls -f *.fq.gz))
# index="/home/micgdu/GenomicData/genomicIndices/Mmusculus/bowtie2/mm39"
# outFolder="/home/micgdu/Analysis/Epidermis/geneExpression_Fan_Nayat/Bowtie2alingments_Nayat"
# 
# for ((i=0;i<=${#samples[@]};i++)); do
#   /home/micgdu/workflows/RNAseq/scripts/bowtie2_SE_singleFile_3may.sh $index "${samples[i]}" $outFolder &
# done
#
# $1 index- path to premade bowtie2 index and index suffix: egz. "mm39"
# $2 sample
# $3 outFolder
# paired end files format:
# basal_ctr_Rep1_MKDL250010461-1A_235HVLLT3_L5_1_val_1.fq.gz
# basal_ctr_Rep1_MKDL250010461-1A_235HVLLT3_L5_2_val_2.fq.gz 

R1="$2"

if [[ "$R1" == *"_R1_001_val_1.fq.gz" ]]; then
  R2="${R1/_R1_001_val_1.fq.gz/_R2_001_val_2.fq.gz}"
elif [[ "$R1" == *"_1_val_1.fq.gz" ]]; then
  R2="${R1/_1_val_1.fq.gz/_2_val_2.fq.gz}"
else
  echo "ERROR: Don't recognize R1 naming: $R1" >&2
  exit 1
fi

bowtie2 --minins 20 --maxins 1200 --dovetail -q --phred33 -p 12 --no-unal -x  $1 -1 $R1 -2 $R2 -S $2.sam --met-file $2.metrics.txt 2> $2.txt
samtools view  -@ 12 -bS $2.sam > $2.bam
samtools sort -@ 12 $2.bam -o $2.sorted_stChr.bam #name with "stChr" for comaptybility
samtools index  -@ 12 -b $2.sorted_stChr.bam $2.sorted_stChr.bai
mv $2.sorted_stChr.bam $2.sorted_stChr.bai $2.txt $3
rm $2.sam
rm $2.bam