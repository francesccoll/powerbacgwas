# PowerBacGWAS: Power calculations for Bacterial GWAS

PowerBacGWAS is a computational pipeline to conduct power calculations for bacterial GWAS. It uses existing collections of bacterial genomes to establish the sample sizes required to detect statistical significant associations for a given genotype frequency and effect size (or phenotype heritability). It supports a range of genomic variation including SNPs, indels, and variation in gene content (pan-genome). Here, we make the code available, and provide installation and [usage](https://github.com/francesccoll/powerbacgwas/wiki#usage) instructions. PowerBacGWAS can be applied to any bacterial population Here we applied it to three different bacterial species: _Enterococcus faecium_, _Klebsiella pneumoniae_, and _Mycobacterium tuberculosis_.

# Docker/Nextflow Installation

The easiest and recommended way to install and run PowerBacGWAS is via its Docker/Nextflow implementation.

You will need to:

1. Install [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) 
2. Install [Nextflow](https://www.nextflow.io/)
3. Download PowerBacGWAS Nextflow files from GitHub:

```console
git clone https://github.com/francesccoll/powerbacgwas/
cd powerbacgwas/nextflow
nextflow run main.nf --help
```

See the [PowerBacGWAS](https://github.com/francesccoll/powerbacgwas/wiki#usage) wiki page for examples of Nextflow commands.

# Manual Installation
PowerBacGWAS consists of a set of Python and R scripts that would work provided that all required dependencies below (both python modules and software) are installed in your local machine.

## Required dependencies

### Software
* [Python3](https://www.python.org/) version >= 3.6.9
* [R](https://www.r-project.org/) version >= 3.6.3
* [PastML](https://github.com/evolbioinfo/pastml) version >= 1.9.20
* [plink](https://zzz.bwh.harvard.edu/plink/download.shtml) version >= PLINK v1.90b6.17 64-bit (28 Apr 2020)
* [GCTA](https://cnsgenomics.com/software/gcta/#Download) version >= version 1.93.2
* [pyseer](https://github.com/mgalardini/pyseer) version >= pyseer 1.3.7-dev
* [bcftools](https://github.com/samtools/bcftools/releases/) version >= 1.9
* [bgzip](https://github.com/samtools/htslib/releases/) version >= 1.9
* [tabix](https://github.com/samtools/tabix) version >= 1.9

### Python Modules
* [cyvcf2](https://github.com/brentp/cyvcf2) version >= 0.20.8
* scipy version >= 1.5.4
* numpy version >= 1.19.4
* pandas version >= 1.1.4
* PyVCF version >= 0.6.8

### R libraries
* optparse >= 1.6.6
* ggplot2 version >= 3.3.1

### From Source

Download the latest release from this github repository or clone it.

```console
git clone https://github.com/francesccoll/powerbacgwas/
```
```console
cd powerbacgwas/
```

As the pipeline uses scripts from [PastML](https://github.com/evolbioinfo/pastml) and [PySeer](https://github.com/mgalardini/pyseer), clone their GitHub directories into the downloaded _powerbacgwas_ folder:
```console
git clone https://github.com/evolbioinfo/pastml
```
```console
git clone https://github.com/mgalardini/pyseer
```

# Usage and Tutorials

Please read the [PowerBacGWAS](https://github.com/francesccoll/powerbacgwas/wiki#usage) wiki page for full usage instructions and tutorials.

# License

PowerBacGWAS is a free software, licensed under [GNU General Public License v3.0](https://github.com/francesccoll/powerbacgwas/blob/main/LICENSE)

# Feedback/Issues

Use the [issues page](https://github.com/francesccoll/powerbacgwas/issues) to report on installation and usage issues.

# Citation
_Not available yet_




