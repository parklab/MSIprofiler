#!/bin/bash

echo "Download fasta files for chrs (hg19)"

mkdir -p chrs_fa

for i in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y
do
	echo $i
	wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr${i}.fa.gz
	gunzip chr${i}.fa.gz
	mv chr${i}.fa ./chrs_fa
done

echo "Fasta files downloaded successfully"

