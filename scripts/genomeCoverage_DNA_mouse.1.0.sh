#!/bin/bash
# script for generates coverage raw & read number (rpm) normalized coverage in the form of begraph file
# non standard chromosomes and chrM are removed - browser compatibility; chrM => avoid excess in ATAC-seq
# first retrieves number of mapped reads from bam alignments "samtools view -c"
# then calculates scalingFactor: 1/(reads number x 1e6)
# additionally script sorts normalized bedGrpah files and generates bigWig files with bedGraphToBigWig programme

#samtools view -b $1 chrM > $1_chrM.bam
samtools index $1
samtools view -@ 10 -b $1 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > $1_stChr.bam
samtools sort -@ 10 $1_stChr.bam -o $1_SortStChr.bam
samtools index -@ 10 $1_SortStChr.bam

bedtools genomecov -ibam $1 -bga > $1.bedGraph
gzip $1.bedGraph

readsNumber=$(samtools view -@ 10 -c $1_SortStChr.bam)
scalingFactor=$(awk "BEGIN {print 1e6/$readsNumber}")
/home/micgdu/kentutils/fetchChromSizes $2 > $2.chrom.sizes
bedtools genomecov -ibam $1_SortStChr.bam -bga -scale $scalingFactor > $1_norm.bedGraph
awk '$1 ~ /chr/' $1_norm.bedGraph > $1_normLim.bedGraph
LC_COLLATE=C sort -k1,1 -k2,2n $1_normLim.bedGraph > $1_Snorm.bedGraph
/home/micgdu/kentutils/bedGraphToBigWig $1_Snorm.bedGraph $2.chrom.sizes $1_Snorm.bw
gzip $1_Snorm.bedGraph

rm $1_stChr.bam
rm $1_SortStChr.bam
rm $1_norm.bedGraph
rm $1_normLim.bedGraph

