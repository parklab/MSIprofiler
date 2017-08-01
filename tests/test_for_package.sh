#!/bin/bash
#cd /n/data1/hms/dbmi/park/DATA/ICGC/WGS/coad-us/COAD-US::c83d38fc-e011-4f75-a100-96513611f3e9
#
#
## get test tumor
#samtools view PCAWG.376a4b9b-94b9-4fd7-b3f0-96712f54428a.bam "22:16050273-16250000" -Sb > test_tumor.bam
#samtools index test_tumor.bam
#
## get test normal
#samtools view ./6eb072a7-5c3b-4adb-b3d3-78303bb8f6a6/PCAWG.9c2006cb-79e6-49c1-993d-dad0ee88edf9.bam "22:16050273-16250000"  -Sb > test_normal.bam
#samtools index test_normal.bam
##17125365
#head -n 1000 /n/data1/hms/dbmi/park/icortes/MSIprofiler/reference_set_22_sorted.txt  > reference_set_22_sorted_test.txt

# phased test
time python /n/data1/hms/dbmi/park/icortes/MSIprofiler/MSIprofiler_2.py --tumor_bam test_tumor.bam  --normal_bam test_normal.bam  --bed germline_calls_22_sel1k.bed --chromosomes 22 --fasta /n/data1/hms/dbmi/park/icortes/MSIprofiler/chrs_fa/  --output_prefix test --mode phased --nprocs 8  --reference_set reference_set_22_sorted_test.txt  --min_coverage 5 --min_MS_length 6 --flank_size 5 --rus 1 2 3 4 5 6

# unphased test
time python /n/data1/hms/dbmi/park/icortes/MSIprofiler/MSIprofiler_2.py --tumor_bam test_tumor.bam  --normal_bam test_normal.bam  --bed germline_calls_22_sel1k.bed --chromosomes 22 --fasta /n/data1/hms/dbmi/park/icortes/MSIprofiler/chrs_fa/  --output_prefix test --mode unphased --nprocs 8  --reference_set reference_set_22_sorted_test.txt  --min_coverage 5 --min_MS_length 6 --flank_size 5 --rus 1 2 3 4 5 6

