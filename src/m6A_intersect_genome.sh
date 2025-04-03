#!/bin/bash

# generate bed file.[chr,start,end,name,score,strand]
tail -n +2 ../data/example_m6A.csv | awk '{print $1,$2,$2+1,$5,$10,$3}'  OFS='\t'  > ../data/example_m6A.bed
bedtools intersect -a ../data/example_m6A.bed -b ../genecode/gencode.v38.annotation.gtf -s -wb > ../data/m6A_intersect_genome.gtf

