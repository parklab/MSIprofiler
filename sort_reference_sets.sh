#!/bin/bash



for i in {{1..22},X,Y}; do
sort -k 2 reference_set_${i}.txt > reference_set_${i}_sorted.txt
#sort -k 2 reference_set_${i}_up_to_5rep_in_flank.txt > reference_set_${i}_up_to_5rep_in_flank_sorted.txt
done
