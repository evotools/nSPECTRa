
# nSPECTRa
## Mutation spectra analyses
*nSPECTRa* is a nextflow workflow that implements several methods to compute the mutation spectra for a given set of variants and multiple genome alignments which includes the ancestral genomes.

## Requirements
The workflow comes with an anaconda environment which delivers most of the dependencies.
However, you need to install:
 1. One of the following package/container managers:
   - [anaconda](https://www.anaconda.com/products/individual)
   - [singularity](https://sylabs.io/) or 
   - [docker](https://www.docker.com/)
 2. [nextflow](https://www.nextflow.io/)
 3. if run with conda, install [gcc]() > 5, tested with version 7.3.0
 4. [beagle v5](https://faculty.washington.edu/browning/beagle/beagle.html#download) or newer jar file. The workflow can download this automatically if not specified.
 5. [relate](https://myersgroup.github.io/relate/) software suite

## Quicker environment installation with mamba
The anaconda environment can take up to several hours to install due to the large number of dependencies it has to retrieve.
You can speed up the installation using [mamba](https://mamba.readthedocs.io/en/latest/), that is an alternative anaconda package manager focusing on speed of installation and dependencies resolution.
To create an environment suitable for *nSPECTRa*, first create a new anaconda environment with [mamba](https://mamba.readthedocs.io/en/latest/) installed:
```
conda create -n mamba -c conda-forge mamba
mamba create -n nspectr -f ./nSpectr/environment.yml
conda activate nspectr
```

Then, simply pass it to nextflow as follow:
```
nextflow run main.nf [PARAMETERS HERE] -profile conda -with-conda '/PATH/TO/ANACONDA/myanaconda/nspectr'
```

## Inputs
To run the workflow effectively, you'll need:
 1. A vcf file with the genotypes for the samples to analyse, tbi-indexed
 2. An HAL alignment archive inclusive of the ancestral genome and the same genome used to call the variants in the VCF
 3. If running relate, the population lists (one text file with the IDs of the samples in each population; e.g. Pop1.txt includes the list of samples in Pop1)
 4. Identifier of the reference genome in the HAL archive
 5. Identifier of the ancestral genome in the HAL archive
 6. An effective population size to use for the imputation

To generate the HAL alignment archive, we suggest to use [cactus](https://github.com/ComparativeGenomicsToolkit/cactus) which automatically infers the ancestral sequence.

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
 2. relate: still under test;
 3. dinuc: method to define adjacent subsequent mutations from [Prendergast et al., 2019](https://academic.oup.com/gbe/article/11/3/759/5299487).

Mutiple choices are possible, and the user is free to use all of them by setting `--algorithm dinuc,mutyper,relate`. 

## Input pre-filtering
We recommend to pre-filter the vcf file to obtain samples with a reasonable coverage (>8-10 mean DP) and with variants with low call rate removed from the dataset (CCR > 90%).