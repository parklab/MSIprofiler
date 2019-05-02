# MSIprofiler [![Build Status](https://travis-ci.com/parklab/MSIprofiler.svg?token=EkzyvwdZ2jcY78ErmS88&branch=master)](https://travis-ci.com/parklab/MSIprofiler) [![codecov](https://codecov.io/gh/parklab/MSIprofiler/branch/master/graph/badge.svg?token=PMpspJNu0Z)](https://codecov.io/gh/parklab/MSIprofiler)

MSIprofiler is a tool to detect microsatellite instability (MSI) from sequencing data.
Its functionalities include:

  - Detection of microsatellites (MS) in the reference genome.
  - Haplotype-specific MSI detection. In this case heterozygous SNPs are used for phasing the microsatellites.
  - Detection of MSI accross the entire genome (not phased).

The underlying principles of MSIprofiler and the details about its validation
using experimental data from TCGA can be found in the following two publications:

**The landscape of microsatellite instability in colorectal and endometrial cancer genomes**
Tae-Min Kim et al. Cell. 2013 (http://www.sciencedirect.com/science/article/pii/S0092867413012919)

and 

**A molecular portrait of microsatellite instability across multiple cancers**
Isidro Cortes-Ciriano et al. Nat. Commun. 2017 (https://www.nature.com/articles/ncomms15180)

Requests for use of the Software for or on behalf of for-profit entities or for any commercial purposes, please contact:
Office of Technology Development
Harvard University
Smith Campus Center, Suite 727E
1350 Massachusetts Avenue
Cambridge, MA 02138 USA
Telephone: (617) 495-3067
E-mail: otd@harvard.edu

# Pre-reqs

MSIprofiler is written entirely in python and has been developed using the following python version:
- **python**: 2.7.6 (default, Aug  3 2015, 17:43:52)  [GCC 4.4.7 20120313 (Red Hat 4.4.7-11)] 

We use pip and virtual enviornments to take care of our dependency management
- [**pip**](https://pip.pypa.io/en/stable/installing/)
- [**virtualenvwrapper**](https://virtualenvwrapper.readthedocs.io/en/latest/install.html#installation) (optional, but recommended)


# Installation

`$ git clone https://github.com/parklab/MSIprofiler.git && cd MSIprofiler`

`$ mkvirtualenv MSIprofiler-env` (optional, but recommended)

`$ pip install -r requirements.txt`

# Running Tests
`$ python tests.py`

## Detection of microsatellites in the human genome
  - MSIprofiler uses a high-confidence reference set of microsatellites derived from the human genome.
 Mapping reads to highly repetitive regions is challenging and the presence of concatenated microsatellites can hamper a correct detection of microsatellite lengths. Therefore, MSIprofiler only considers MS repeats whose flanking regions do not contain MS repeats of more than 5 bases (e.g. AAAAAA). This permits to correctly align the flaking regions and reduce the false positive rate.

To generate the reference set of MS repeats, first download the fasta sequences for the chromosomes by running the file:

```sh
$ ./scripts/download_chromosomes_fa.sh
```

- Once the fasta sequences are downloaded, run the following two scripts in the root directory of MSIprofiler: 
```sh
$ python scripts/get_reference_set_from_fasta.py
$ ./scripts/sort_reference_sets.sh
```

These scripts will generate one file per chromosome containing the reference sets of MS repeats (coordinate-sorted).

Make sure to add the current directory to your PYTHONPATH:
export PYTHONPATH=$PYTHONPATH:$PWD

<!-- ### Detection of microsatellites

![examples_detection](https://user-images.githubusercontent.com/5588266/28093196-42fedc6a-6664-11e7-9e8a-04d555a88e7e.png)
-->

## MSIprofiler parameters

Once the reference sets are ready, MSIprofiler can be used.

Information on the parameters can be accessed by typing:
```sh
$ python msi_profiler.py --help
```

Example of usage:
```sh
python path_to_msi_profiler/msi_profiler.py  --tumor_bam tumor.bam  --normal_bam normal.bam --bed hets_SNPs_chr22.bed --chromosomes 22  --fasta path_to_MSIprofiler/chrs_fa/ --output_prefix example_chr22  --mode unphased --nprocs 8  --reference_set path_to_reference_sets_folder --min_coverage 8 --min_MS_length 6 --flank_size 10 --rus 1 2 3 4 5 6  --tolerated_mismatches 0
```

## Haplotype-specific detection of MSI

MSIprofiler can detect haplotype-specific MSI by phasing microsatellites and heterozygous SNPs detected in the germline (i.e. normal/control sample).
The steps followed by MSIprofiler for the detection of haplotype-specific MSI are:

- Detect whether the reads covering each allele of the input germline heterozygous SNP also contain MS repeats that are present in the reference set. This step is applied across all input SNPs to both the normal and the tumor/case samples using the input bam files.
- Compare the distribution of MS lengths (i.e. read length distributions) in the normal and tumor/case samples using the Kolmogorov-Smirnov test.

To calculate haplotype-specific MSI, the parameter "mode" needs to be set to 'phased'.

The bed files containing the heterozygous SNPs in 0-based coordinates need to have the following columns: chr, start, end, ref and alt.
For instance, an entry would look like:
7	20607	20608	A	G

### Example of usage
```sh
python msi_profiler.py --tumor_bam test_tumor.bam  --normal_bam test_normal.bam  --bed germline_het_SNPs.bed --chromosomes 21 22 X --fasta ./chrs_fa/  --output_prefix example_unphased --mode unphased --nprocs 2  --reference_set path_to_reference_sets_folder  --min_coverage 8 --min_MS_length 6 --flank_size 5 --rus 1 2 3 4 5 6 
```

### Example of output
22	16116649	AGAAGAAG	3	8	0.3	G	16116721	8,8,8,8,8,8	8,8,8,8,8,8,8	1.0
The columns correspond to:
1. chromosome
2. start of the MS repeat (0-based)
3. microsatellite repeat
4. microsatellite repeat motif length
5. microsatellite length in the reference genome
6. GC contents in the flanking regions in the reference genome (2000 bases)
7. SNP allele
8. SNP start (0-based)
9. length of the MS repeats detected in the normal/control sample
10. length of the MS repeats detected in the tumor/case sample
11. Kolmogorov-Smirnov P value

## Detection of MSI (unphased) by comparing read-length distributions across both alleles

MSIprofiler detecs MSI directly from the sequencing data without considering phasing information in a similar manner as
the experimental assays customarily used for MSI detection (e.g. capillary sequencing-based fragment length assay). 
To this aim, the lengths of a given repeat present in the reference set are measured using the reads from both the normal and tumor/case samples. Next, the Kolmogorov-Smirnov test is used to test for a significant difference in the distributions. 
These steps are applied across all MS repeats in the reference set sufficiently covered by the sequencing data (default is 10 reads).

In this case, the parameter "mode" needs to be set to 'unphased'.

### Example of usage
```sh
python msi_profiler.py --tumor_bam test_tumor.bam  --normal_bam test_normal.bam  --bed None --chromosomes 22 --fasta ./chrs_fa/  --output_prefix example_unphased --mode unphased --nprocs 2  --reference_set path_to_reference_sets_folder  --min_coverage 8 --min_MS_length 6 --flank_size 5 --rus 1 2 3 4 5 6 
```

### Example of output
22	16138703	16138712	ACAAGACAAG	5	10	0.4	10,10,10,10,10,10	10,10,10,10,10,10,10,10,10,10,10,10	1.0	high
The columns correspond to:
1. chromosome
2. start of the MS repeat (0-based)
3. end of the MS repeat (0-based)
4. MS repeat
5. repeat motif length
6. length of the repeat in the reference genome
7. length of the MS repeat in the normal/control sample
8. length of the MS repeat in the tumor/case sample
9. Kolmogorov-Smirnov P value
10. confidence (see below)

### Confidence of the unphased calls

The assumptions of this methodology are that:
- both alleles are at the same copy number (usually 1:1)
- both alleles are represented evenly in the sequencing data; i.e. there is no allelec bias during the library preparation or sequencing steps.

These assumptions are however not always met in tumor samples due to aneuploidy (e.g., focal amplifications or deletions, whole-chromosome gains/losses, etc..). Moreover, intratumor heterogeneity represents a potential confounding factor if an MS repeat is altered in only a subset of the cancer cells (i.e. subclonal mutation). We note that intratumor heterogeneity, as well as tumor purity, are the main reasons why it is challenging to define a pan-cancer statistical model to genotype microsatellites in tumor samples in a similar way as it is done for normal cases.

<!--In the case of single-cell sequencing data where the polymerase phi29 is generally used for amplification (e.g. multiple displacement amplification or MDA), large regions of the genome are not amplified evenly (allelic imbalance or even allelic dropout). -->

These issues are relevant for MSI detection, as the imbalance between alleles represents a source of false positives for heterozygous microsatellites. To account for this, we assign a confidence level to each of the unphased calls based on whether the 
MS repeat under consideration is homozygous or heterozygous in the germline. 

If an MS repeat is homozygous in the germline and it is mutated in any of the copies present in the tumor/case sample,
the read length distribution will change and the mutation will be detected using the pipeline presented above. 
This is due to the fact that the read length distribution for a homozygous microsatellite in the normal sample is unimodal (or close to unimodal depending on the stutter noise).

By contrast, in the case of heterozygous MS repeats the read length distribution is bimodal, with each mode corresponding
to one alelle. If one allele is lost in the tumor/case sample (e.g. one-chromosome loss or focal deletion), the read length distributions between the normal and tumor samples would differ significantly, hence leading to a false positive. This situation can also happen if one copy has been amplified hundreds of times (e.g. in a double minute chromosome). In such a case, the high imbalance in copy number between the two alleles would be refleceted in the seqeuncing data, as the faction of reads coming from the unamplified allele would be underrepresented. Hence, the read length distribution would differ significantly in this case even if there is no real mutation.

<!--Ideally, only MSI calls calculated for MS repeats located in regions at copy number of 2 without loss of heterozygosity should be considered. -->
Given that copy number information is not always available,
we assign high-confidence to the unphased calls made on MS repeats that are homozygous in the germline,
and low confidence to those calls made on heterozygous MS repeats.
We consider that an MS repeat is homozygous in the normal if at least 70% of the reads support the same MS length. 

<!-- Overall, we recommend using phased calls whenever possible.-->

# A comment on the number of mismatches in the flanking regions and the length of these

Based on our experience, we recommend to consider flanking regions of at least 10 bases and allow for no mismatches in these to get conservative calls. We recommend to inspect candidate mutations manually when using more lenient parameter values.

The following example illustrates why allowing mismatches in the flanks can lead to false positive calls (and also illustrates why refining the reference sets is necessary). 

Consider the microsatellite repeat:
chr9	89149744	89149756	AAAAAAAAAAAAA intergenic; length: 13

The reads found in one of the alleles of a given patient look like (base qualities and additional information have been removed for the sake of clarity) :

<!-- HSQ700642:192:C13FVACXX:2:2214:4777:60283   99  9   89149699    60  100M    =   89149723    124  -->                               
ATTGCACAATACATGACCTAATGGAAATGTGAGAATA **TTTTAGTG** *AAAAAAAAAAAAA* **TAAAAAGA** AGCAGCAAAGATCCAACCAAATGAGATCCATATG

The length of this repeat supported by the data is 13. The right-hand flanking region (i.e. 8 downstream bases in this example) is **TAAAAAGA**.

The reads supporting the other allele look like:
<!-- SQ700642:208:D1D6WACXX:2:1116:11062:89332  83  9   89149707    60  100M    =   89149522    -285   -->                
ATACATGACCTAATGGAAATGTGAGAATA **TTTTAGTG**  *AAAAAAA* **TAAAAATA** AAAAGAAGCAGCAAAGATCCAACCAAATGAGATCCATATGGGATGGGT   

The estimated length for this repeat in these reads is 7, due to the SNP located in the middle of the repeat. In this case, the right-hand flanking region is **TAAAAATA**.

One of the SNP alleles (i.e. T) interrupts the polyA motif. The 8 bases on the right, only differ in one base with the reference flank region:  reference: TAAAAA**G**A _vs_ TAAAAA**T**A.

If we allowed for one mismatch in the flanking region, the estimated lengths for this microsatellite in this patient would be 7 (SNP allee T) and 13 (SNP allele A). 

Thus, allowing for a mismatch can lead to an incorrect estimation of the true length of the microsatellite.


# Notes
Currently, MSIprofiler considers by default the hg19 assembly of the human genome.
However, MSIprofiler can be used with other assemblies if the bam files and the fasta sequences from which the reference sets are derived are concordant.

<!--
# How to cite
The details of MSIprofiler have been published in:
XXX
-->

# Contact
If you have any questions or suggestions please contact us: 
- Isidro Cortes Ciriano: isidrolauscher at gmail.com  or isidro\_cortesciriano at hms.harvard.edu
- Peter J Park: peter\_park at hms.harvard.edu

# Acknowledgements
We would like to acknowledge the support of the Chan Zuckerberg Initiative for their Collaborative Computational Tools for the Human Cell Atlas program.



