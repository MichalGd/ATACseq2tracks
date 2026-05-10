#!/bin/bash
# reads deduplication with picard
# requires bam file with read group tags; samtools addreplacerg

java -Xmx32G -jar /home/micgdu/software/picard.jar MarkDuplicates \
    -INPUT $1/$3 \
    -OUTPUT $2/$3_dedup.bam \
    -METRICS_FILE $2/$3_dedup_rep.txt \
    -OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
    -REMOVE_DUPLICATES true \
    -ASSUME_SORT_ORDER coordinate \
    -CREATE_INDEX true \
    -TMP_DIR /tmp