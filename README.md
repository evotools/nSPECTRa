
# nSPECTRa
## Mutation spectra analyses
*nSPECTRa* is a nextflow workflow that implements several methods to compute the mutation spectra for a given set of variants and multiple genome alignments which includes the ancestral genomes.

## Requirements
The workflow comes with an anaconda environment which delivers most of the dependencies.
However, you need to install:
 1. [anaconda](https://www.anaconda.com/products/individual)
 2. [nextflow](https://www.nextflow.io/)
 3. if run with conda, install [gcc]() > 5, tested with version 7.3.0
 4. [beagle v5](https://faculty.washington.edu/browning/beagle/beagle.html#download) or newer jar file. The workflow can download this automatically if not specified.
 5. [relate](https://myersgroup.github.io/relate/) software suite
Most of the remaining dependencies are downloaded by *nSPECTRa* at runtime.

## Inputs
To run the workflow effectively, you'll need:
 1. A vcf file with the genotypes for the samples to analyse, tbi-indexed
 2. An HAL alignment archive inclusive of the ancestral genome and the same genome used to call the variants in the VCF
 3. If running relate, the population lists (one text file with the IDs of the samples in each population; e.g. Pop1.txt includes the list of samples in Pop1) or a poplist prepared as described in the relate documentation
 4. Identifier of the reference genome in the HAL archive
 5. Identifier of the ancestral genome in the HAL archive
 6. An effective population size estimate to use at imputation time

To generate the HAL alignment archive, we suggest to use [cactus](https://github.com/ComparativeGenomicsToolkit/cactus) which automatically infers the ancestral sequence.

## Workflow overview
The workflow constitutes of five separate components, repesented in figure below:

![Flowchart](https://github.com/evotools/nSPECTRa/blob/bc8d089c75f8ca7625e8decb64058c1b6f230c5b/imgs/WorkflowComponents.png)

In short, there are six separate units:
1. VCF preprocessing: the vcf file gets imputed using either Beagle or SHAPEIT4, and then annotated using the VEP.
2. Ancestral genome generation: the workflow generate a version of the reference genome carrying the ancestral allele estimated by [cactus](https://github.com/ComparativeGenomicsToolkit/cactus)
3. Exclude constrained elements: the workflow extract the constrained elements using halPhyloP from the [HAL suite](https://github.com/ComparativeGenomicsToolkit/hal) alignments.
4. K-mer mutation spectra: compute the mutation spectra for each individual at different K-mers using [mutyper](https://github.com/harrispopgen/mutyper/)
5. Sequential dinucleotide mutation: compute the sequential dinucleotide mutation (SDM) spectra for the different individuals
6. Compute the mutation spectra evolution: compute the trend of the different mutation classes in the different populations using [relate](https://myersgroup.github.io/relate/)

The different components can be switches off depending on the analytical requirements and public resources available.

## Algorithm choice
The user can choose between two different imputation algorithms and three algorithms to compute the mutation spectra.

### Imputation
The user can choose to impute either with [shapeit]() v4 or [beagle]() v5 or newer. 
Both are supported, but only shapeit is included in the anaconda environment, which is the reason why it is the default choice.
In case the user prefers to use beagle, this needs to be downloaded manually and the path provided with the flag `--beagle`. 
If `--beagle` is not provided, the workflow will download the version 5.2 automatically.  
The user can choose what algorithm to use, but keep in mind that they are mutually exclusive, and shapeit will be prioritized over beagle.

### Mutation spectra
The software currently supports three software to compute the mutation spectra:
 1. mutyper: the default choice, it can compute the spectra for different K-mer size;
 2. relate: used to define the mutation spectra changes over time and an accurate estimate of the effective population size;
 3. sdm: method to define sequential dinucleotide polymorphisms, MNPs and adjacent SNPs from [Prendergast et al., 2019](https://academic.oup.com/gbe/article/11/3/759/5299487).

Mutiple choices are possible, and the user is free to use all of them by setting `--algorithm sdm,mutyper,relate`. 

## Input pre-filtering
We recommend to pre-filter the vcf file to obtain samples with a reasonable coverage (>8-10 mean DP) and with variants with low call rate removed from the dataset (CCR > 90%) and minor allele count MAC >= 2.