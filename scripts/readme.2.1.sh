#changes fastq2tracks 1.5 => 2.1

1.a) changing major problem => picard_deduplication.1.0
java -Xmx128G -jar /home/micgdu/software/picard.jar MarkDuplicates => java -Xmx32G -jar /home/micgdu/software/picard.jar MarkDuplicates

1.b) adjusting batch script => name of subscript
java -Xmx128G -jar /home/micgdu/software/picard.jar MarkDuplicates => java -Xmx32G -jar /home/micgdu/software/picard.jar MarkDuplicates

2.a) changing bowtie2_dovetail_pairedEnd_MMusculus.1.5; remoing filtering with chromsomes => possibly causes faulty picard filtering
original; threads 15 =>12 ; adding multithredian to samtools

"bowtie2 --minins 20 --maxins 1200 --dovetail -q --phred33 -p 15 --no-unal -x  $1 -1 $R1 -2 $R2 -S $2.sam --met-file 2> $2.txt
samtools view -bS $2.sam > $2.bam
samtools sort $2.bam -o $2.sorted.bam
samtools index -b $2.sorted.bam $2.sorted.bai
samtools view -b $2.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY chrM > $2.sorted_stChr.bam
samtools index -b $2.sorted_stChr.bam $2.sorted_stChr.bai
mv $2.sorted_stChr.bam $2.sorted_stChr.bai $2.txt $3
rm $2.sam
rm $2.bam
rm $2.sorted.bam
rm $2.sorted.bai"

new:
"bowtie2 --minins 20 --maxins 1200 --dovetail -q --phred33 -p 12 --no-unal -x  $1 -1 $R1 -2 $R2 -S $2.sam --met-file "$2.metrics.txt 2> $2.txt
samtools view  -@ 12 -bS $2.sam > $2.bam
samtools sort -@ 12 $2.bam -o $2.sorted_stChr.bam #name with "stChr" for comaptybility
samtools index  -@ 12 -b $2.sorted.bam $2.sorted_stChr.bai
mv $2.sorted_stChr.bam $2.sorted_stChr.bai $2.txt $3
rm $2.sam
rm $2.bam"

2.b) bowtie2_dovetail_pairedEnd_Hsapiens.2.1 => now different name the same script; before filtering through chromsomes was species specific

2.c & 2.d) adjusting names of batch scripts respetively

