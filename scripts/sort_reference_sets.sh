#!/bin/bash



for i in {{1..22},X,Y}; do
sort -k 2 reference_set_${i}.txt > reference_set_${i}_sorted.txt
done
