# PowerBacGWAS: Power calculations for Bacterial GWAS


# Installation
PowerBacGWAS consists of a set of python scripts that should work provided that all required dependencies below (both python modules and software) are installed in your local machine.

## Required dependencies

### Software
* [Python3](https://www.python.org/) version >= 3.6.9
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

