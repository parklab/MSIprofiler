# MSIprofiler

MSIprofiler is a tool to detect microsatellite instability (MSI) from sequencing data.
Its functionalities include:

  - Detection of microsatellites (MS) in the reference genome
  - Haplotype-specific MSI detection using heterozygous SNPs for phasing the microsatellites
  - Detection of MSI accross the entire genome (not phased)

The underlying principles of MSIprofiler and the details about its validation
using experimental data from TCGA can be found in the following two publications:

**The landscape of microsatellite instability in colorectal and endometrial cancer genomes**
Tae-Min Kim et al. Cell. 2013 (http://www.sciencedirect.com/science/article/pii/S0092867413012919)

and 

**A molecular portrait of microsatellite instability across multiple cancers**
Isidro Cortes-Ciriano et al. Nat. Commun. 2017 (https://www.nature.com/articles/ncomms15180)

MSIprofiler is free for academic use only. 
For non-academic use, please email Dr. Tatiana Demidova-Rice at Harvard University Office of Technology Development (tatiana\_demidova-rice@harvard.edu)

# Requirements

MSIprofiler is written entirely in python and has been developed using the following python and library versions:

- python: 2.7.6 (default, Aug  3 2015, 17:43:52)  [GCC 4.4.7 20120313 (Red Hat 4.4.7-11)] 
- numpy: 1.10.1
- scipy: 0.17.0
- csv: 1.0
- multiprocessing: 0.70a1
- argparse: 1.1
- pysam: 0.10.0
- bisect: 2.1


# Usage

First clone the repository to your computer:

```sh
$ git clone https://github.com/parklab/MSIprofiler.git
```

## Detection of microsatellites in the human genome
  - MSIprofiler uses a refined reference set of microsatellites derived from the human genome. Given that mapping reads to highly repetitive regions is challenging and the presence of concatenated microsatellites can hamper a correct detection of microsatellite lengths, MSIprofiler only considers a refined set of MS repeats whose flanking regions do not contain MS repeats. This permits to correctly align the flaking regions. 
MSIprofiler discards soft-clipped bases.
  To generate the reference set of MS repeats, first download the fasta sequences for the chromosomes by running the file:

	  ```sh
	  $ ./download_chromosomes_fa.sh
	  ```

- Once the fasta sequences are downloaded, run in the root of the directory of MSIprofiler the following script: 

	  ```sh
	  $ python get_reference_set_from_fasta.py
	  ```
This will generate one file per chromosome containing the refined reference sets of MS repeats.

Make sure to add the current directory to your PYTHONPATH:
export PYTHONPATH=$PYTHONPATH:$PWD

Note: the reference set coordinated are 1-based.

The refinement consists of discarding MS repeats whose flanking regions also comprise MS repeats of at most 5 bases.  


## MSIprofiler parameters

Once the reference sets are ready, MSIprofiler can be used.

Information on the parameters can be accessed by typing:
```sh
$ python MSIprofiler.py --help
```

- mode
- genomic_region

Parallel 



## Haplotype-specific detection of MSI

MSIprofiler can detect haplotype-specific msi by phasing microsatellites and heterozygous SNPs detected in the germline.

bed file in 0-based coordinates

The steps followed by MSIprofiler for the detection of haplotype-specific MSI are:

- Detect whether the reads covering each germline heterozygous SNP also contain MS repeats that are present in the reference set. 
This step is applied to both the normal and tumor/case samples using the input bam files.
- Compare the distribution of MS lengths (i.e. read length distributions) in the normal and tumor/case samples using the KS test

To calculate haplotype-specific MSI, the parameter "mode" needs to be set to 'both' or 'phased'.


## Detection of MSI (unphased) by comparing read-length distributions across both alleles
- msi directly from the sequencing data without considering phasing information in a similar manner as
the experimental assays (e.g. capillary sequencing-based fragment length assay). 
To this aim, the length of the repeats is measured from the reads in both the normal and tumor/case samples.
Next, the Kolmogorov-Smirnov test is used to test for a significant difference in the distributions.


*  Confidence of the calls

The assumption of this methodology are that:
- both alleles are at the same copy number (usually 1:1)
- both alleles are represented evenly; i.e. both alleles are amplified and sequences evenly.

These assumptions are however not met in tumor samples due to aneuploidy (e.g., focal amplifications or deletions, whole-chromosome gains/losses, etc..).
Moreover, intratumor heterogeneity represents a potential confounding factor if an MS repeat is altered in only a subset of populations of cancer cells (i.e. subclonal mutation).

In the case of single-cell sequencing data where the polymerase phi29 was used for amplification (e.g. multiple displacement amplification or MDA), 
large regions of the genome are not amplified evenly (allelic imbalance). 

The imbalance between alleles represents a source of false positives for heterozygous microsatellites. 
To account for this, we assign a confidence level for each of the unphased calles based on whether the 
microsatellite repeat under consideration is homozygous or heterozygous in the germline. 

If a MS repeat is homozygous in the germline and is mutated in any of the copies present in the tumor/case sample,
the mutation will be detected using the pipeline presented above even if some copies are lost. 
This is due to the fact that the loss of one or multiple copies does not affect the read length distribution, 
as this is unimodal (Fig XX).

By contrast, in the case of heterozygous SNPs, the distribution of read lengths is altered if one of the alleles is lost (e.g. one-chromosome loss or focal deletion).
In such a case, a significant difference between the read length distributions of the germline and tumor/case
samples would represent a false positive (Fig XXX). 

We thus assign high-confidence to the calls made on MS repeats that are homozygous in the germline,
and low confidence to those calls made on heterozygous MS repeats.
To determine whether an MS repeat is homo- or heterozygous in a given normal sample (i.e. in the germline)
we test whether...XXX



copy number alterations in cancer (relationship to confidence)
imposibilidad de statistical model due to heterogeneity (tal vez se puede controlar purity)


# Best practices for single-cell sequencing data
only phased
this is due to the allelic dropout (in MDA)

<!--
# Note on the number of tolerated mismatches and size of the flanking regions
The following example illustrates the importance of a careful choice of the size of the flanking region and the number of tolerated mismatches in these.
Consider the microsatellite repeat:
chr9	89149744	89149756	AAAAAAAAAAAAA	mono intergenic 13

reads found in a patient (base qualities and additional information have been removed for the sake of clarity).


HSQ700642:192:C13FVACXX:2:2214:4777:60283   99  9   89149699    60  100M    =   89149723    124                                                                              
ATTGCACAATACATGACCTAATGGAAATGTGAGAATA *TTTTAGTG*  *AAAAAAAAAAAAA* *TAAAAAGA* AGCAGCAAAGATCCAACCAAATGAGATCCATATG
Detected MS length: 13

SQ700642:208:D1D6WACXX:2:1116:11062:89332  83  9   89149707    60  100M    =   89149522    -285                                                                             

ATACATGACCTAATGGAAATGTGAGAATA *TTTTAGTG*  *AAAAAAA* *TAAAAATA* AAAAGAAGCAGCAAAGATCCAACCAAATGAGATCCATATGGGATGGGT   
Detected MS length: 7
due to a SNP in the middle of the read.
Given that we allowed one mismatch in the flanking region
-->


# Notes
Currently, MSIprofiler only considers the hg19 assembly of the human genome.
However, can be used with other assemblies if the bam files and the reference set are concordant.

# How to cite
The details of MSIprofiler have been published in:
XXX




