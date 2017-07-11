# MSIprofiler

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

MSIprofiler is free for academic use **only**. 
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
$ cd MSIprofiler
```

## Detection of microsatellites in the human genome
  - MSIprofiler uses a high-confidence reference set of microsatellites derived from the human genome.
 Mapping reads to highly repetitive regions is challenging and the presence of concatenated microsatellites can hamper a correct detection of microsatellite lengths. Therefore, MSIprofiler only considers MS repeats whose flanking regions do not contain MS repeats of more than 5 bases (e.g. AAAAAA). This permits to correctly align the flaking regions and reduce the false positive rate.
In addition, MSIprofiler discards soft-clipped bases.

To generate the reference set of MS repeats, first download the fasta sequences for the chromosomes by running the file:

```sh
$ ./download_chromosomes_fa.sh
```

- Once the fasta sequences are downloaded, run the following two scripts in the root directory of MSIprofiler: 
```sh
$ python get_reference_set_from_fasta.py
$ ./sort_reference_sets.sh
```

These scripts will generate one file per chromosome containing the reference sets of MS repeats (coordinate-sorted).

Make sure to add the current directory to your PYTHONPATH:
export PYTHONPATH=$PYTHONPATH:$PWD

Note: the coordinates of the reference MS sets are 1-based.


## MSIprofiler parameters

Once the reference sets are ready, MSIprofiler can be used.

Information on the parameters can be accessed by typing:
```sh
$ python MSIprofiler.py --help

usage: MSIprofiler.py [-h] --tumor_bam TUMOR_BAM --normal_bam NORMAL_BAM --bed BED
                  --fasta FASTA --reference_set REFERENCE_SET --output_prefix
                  OUTPUT_PREFIX --mode MODE --genomic_region GENOMIC_REGION
                  --nprocs NPROCS [-ru RUS] [--min_MS_length MIN_MS_LENGTH]
                  [--max_MS_length MAX_MS_LENGTH]
                  [--mapping_quality MAPPING_QUALITY]
                  [--flank_size FLANK_SIZE] [--min_coverage MIN_COVERAGE]
                  [--tolerated_mismatches TOLERATED_MISMATCHES]

optional arguments:
  -h, --help            show this help message and exit
  --tumor_bam TUMOR_BAM
                        Tumor or case (e.g. single cell) bam file name
  --normal_bam NORMAL_BAM
                        Normal or control bam file name
  --bed BED             Input bed file containing heterozygous SNPs. Those of genotype 0/1 are preferred. The input bed files need to be in 0-based coordinates.
  --fasta FASTA         Input fasta reference file name
  --reference_set REFERENCE_SET
                        Input reference set of microsatellites
  --output_prefix OUTPUT_PREFIX
                        Prefix for the output files
  --mode MODE           The value of this parameter sets whether MSIprofiler will detect MSI focusing only on microsatellites phased with germline SNPs (phased), all microsatellites contained in the reference sets that can be detected in the input bam files (unphased), or both (both).
  --nprocs NPROCS       Number of processes to be launched
  -ru RUS               MS repeat units to be considered (e.g. mono, di, tri, tetra, ..). Specify these as integers (e.g. 1, 2, 3, 4, ..).
  --min_MS_length MIN_MS_LENGTH
                        Minimum length of microsatellites to be considered.
                        Minimum available is 6; default is 10.
  --max_MS_length MAX_MS_LENGTH
                        Maximum length of microsatellites to be considered.
                        Maximum available is 60; default is 60.
  --mapping_quality MAPPING_QUALITY
                        Minimum mapping quality. Default is 40.
  --flank_size FLANK_SIZE
                        Minimum length of the flanking regions. Default is 10
  --min_coverage MIN_COVERAGE
                        Minimum coverage at each MS locus -both in the case
                        and control bams-. Default is 10
  --tolerated_mismatches TOLERATED_MISMATCHES
                        Maximum number of tolerated mismatches in the flanking
                        regions. Default is 0
```

## Haplotype-specific detection of MSI

MSIprofiler can detect haplotype-specific MSI by phasing microsatellites and heterozygous SNPs detected in the germline (i.e. normal/control sample).
The steps followed by MSIprofiler for the detection of haplotype-specific MSI are:

- Detect whether the reads covering each allele of the input germline heterozygous SNP also contain MS repeats that are present in the reference set. This step is applied across all input SNPs to both the normal and the tumor/case samples using the input bam files.
- Compare the distribution of MS lengths (i.e. read length distributions) in the normal and tumor/case samples using the Kolmogorov-Smirnov test.

To calculate haplotype-specific MSI, the parameter "mode" needs to be set to 'both' or 'phased'.


## Detection of MSI (unphased) by comparing read-length distributions across both alleles

MSIprofiler detecs MSI directly from the sequencing data without considering phasing information in a similar manner as
the experimental assays customarily used for MSI detection (e.g. capillary sequencing-based fragment length assay). 
To this aim, the length of a given repeat present in the reference set is measured using the reads from both the normal and tumor/case samples. Next, the Kolmogorov-Smirnov test is used to test for a significant difference in the distributions. These steps are applied across all MS repeats in the reference set sufficiently covered by the sequencing data (default is 10 reads).

To calculate haplotype-specific MSI, the parameter "mode" needs to be set to 'both' or 'unphased'.


*Confidence of the calls

The assumptions of this methodology are that:
- both alleles are at the same copy number (usually 1:1)
- both alleles are represented evenly in the sequencing data; i.e. there is no allelec bias during the library preparation or sequencing steps.

These assumptions are however not always met in tumor samples due to aneuploidy (e.g., focal amplifications or deletions, whole-chromosome gains/losses, etc..). Moreover, intratumor heterogeneity represents a potential confounding factor if an MS repeat is altered in only a subset of the cancer cells (i.e. subclonal mutation). We note that intratumor heterogeneity, as well as tumor purity, are the main reasons why it is challenging to define a pan-cancer statistical model to genotype microsatellites in tumor samples in a similar way as it is done for normal cases.

In the case of single-cell sequencing data where the polymerase phi29 is generally used for amplification (e.g. multiple displacement amplification or MDA), large regions of the genome are not amplified evenly (allelic imbalance or even allelic dropout). 

These issues are relevant for MSI detection, as the imbalance between alleles represents a source of false positives for heterozygous microsatellites. To account for this, we assign a confidence level to each of the unphased calls based on whether the 
MS repeat under consideration is homozygous or heterozygous in the germline. 

If an MS repeat is homozygous in the germline and it is mutated in any of the copies present in the tumor/case sample,
the read length distribution will change and the mutation will be detected using the pipeline presented above. 
This is due to the fact that the read length distribution for a homozygous microsatellite in the normal sample is unimodal (or close to unimodal depending on the stutter noise).

By contrast, in the case of heterozygous MS repeats the read length distribution is bimodal, with each mode corresponding
to one alelle. If one allele is lost in the tumor/case sample (e.g. one-chromosome loss or focal deletion), the read length distributions between the normal and tumor samples would differ significantly, hence leading to a false positive. This situation can also happen if one copy has been amplified hundreds of times (e.g. in a double minute chromosome). In such a case, the high imbalance in copy number between the two alleles would be refleceted in the seqeuncing data, as the faction of reads coming from the unamplified allele would be underrepresented. Hence, the read length distribution would differ significantly in this case even if there is no real mutation.

Ideally, only MSI calls calculated for MS repeats located in regions at copy number of 2 without loss of heterozygosity should be considered. Given that copy number information is not always available,
we assign high-confidence to the unphased calls made on MS repeats that are homozygous in the germline,
and low confidence to those calls made on heterozygous MS repeats.
We consider that an MS repeat is homozygous in the normal if at least 70% of the reads support the same MS length. 


# Best practices for single-cell sequencing data

We recommend using phased calls whenever possible,

# A comment on the number of mismatches in the flanking regions and the length of these

Based on our experience, we recommend to consider flanking regions of at least 10 bases and allow for no mismatches in these. 
The following example illustrates why allowing mismatches in the flanks can lead to false positive calls (and also illustrates why refining the reference sets is necessary). 

Consider the microsatellite repeat:
chr9	89149744	89149756	AAAAAAAAAAAAA intergenic; length: 13

The reads found in one of the alleles of a given patient looked like (base qualities and additional information have been removed for the sake of clarity) :

HSQ700642:192:C13FVACXX:2:2214:4777:60283   99  9   89149699    60  100M    =   89149723    124                                  
ATTGCACAATACATGACCTAATGGAAATGTGAGAATA <span style="color:blue">some *TTTTAGTG* text</span>  *AAAAAAAAAAAAA* *TAAAAAGA* AGCAGCAAAGATCCAACCAAATGAGATCCATATG
The length of this repeat supported by the data is: 13

whereas those phased with the other allele looked like:
SQ700642:208:D1D6WACXX:2:1116:11062:89332  83  9   89149707    60  100M    =   89149522    -285                   
ATACATGACCTAATGGAAATGTGAGAATA *TTTTAGTG*  *AAAAAAA* *TAAAAATA* AAAAGAAGCAGCAAAGATCCAACCAAATGAGATCCATATGGGATGGGT   
The length of this repeat supported by the data is: 7

due to a SNP in the middle of the read.
If we allowed for one mismatch in the flanking region, the length of this microsatellite in this patient would be 7 and 13. 
The right-hand flanking region is *TAAAAAGA*. 
One of the SNP alleles (i.e. T; shown in green) interrupts the polyA motif. The 8 bases on the right, only differ in one base with the reference flank region. 
Thus, allowing for a mismatch can lead to an incorrect estimation of the true length of the microsatellite.


# Notes
Currently, MSIprofiler considers by default the hg19 assembly of the human genome.
However, MSIprofiler can be used with other assemblies if the bam files and the fasta sequences from which the reference sets are derived are concordant.

<!--
# How to cite
The details of MSIprofiler have been published in:
XXX
-->



