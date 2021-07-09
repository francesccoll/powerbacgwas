# Use a LTS Ubuntu version as parent image
FROM ubuntu:20.04
# NOTE: without interactive dialogue needed to install r-base
ARG DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# metadata
LABEL base.image="ubuntu:20.04"
LABEL container.version="0.1"
LABEL software="powerbacgwas"
LABEL software.version="1.0"
LABEL description="PowerBacGWAS is a computational pipeline to conduct power calculations for bacterial GWAS"
LABEL website="https://github.com/francesccoll/powerbacgwas"
LABEL lisence="https://github.com/francesccoll/powerbacgwas/blob/main/LICENSE"
LABEL maintainer="Francesc Coll"
LABEL maintainer.email="francesc.coll@lshtm.ac.uk"


# Install general dependencies
# NOTE: 'libbz2-dev liblzma-dev ncurses-dev zlib1g-dev libcurl4-openssl-dev libssl-dev' needed to install HTSlib and bcftools
# NOTE: 'gfortran libopenblas-dev liblapack-dev' needed to install scipy (https://github.com/scipy/scipy/issues/9005)
RUN apt-get update && apt-get install --no-install-recommends -y \
    curl \
    build-essential \
    automake \
    pkg-config \
    zlib1g-dev \
    unzip \
    autoconf \
    check \
    libtool \
    libsubunit-dev \
    git \
    wget \
    python3.8 \
    python3.8-dev \
    python3-setuptools \
    python3-pip \
    r-base \
    libbz2-dev liblzma-dev ncurses-dev zlib1g-dev libcurl4-openssl-dev libssl-dev libcrypto++6 \
    libfreetype-dev libfreetype6-dev libpng-dev libpng-tools gfortran gfortran-9 libgfortran-9-dev libopenblas-dev liblapack-dev


# Installing htslib & bcftools (includes bgzip version >= 1.9 & tabix version >= 1.9)

ARG btversion=1.9
RUN wget https://github.com/samtools/htslib/releases/download/${btversion}/htslib-${btversion}.tar.bz2 && \
  tar -vxjf htslib-${btversion}.tar.bz2 && \
  cd htslib-${btversion} && \
  make install && \
  wget https://github.com/samtools/bcftools/releases/download/${btversion}/bcftools-${btversion}.tar.bz2 && \
  tar -vxjf bcftools-${btversion}.tar.bz2 && \
  cd bcftools-${btversion} && \
  make install


# Install R and R libraries

RUN R -e "install.packages('optparse',dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
  R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')"


# Installing PastML version == 1.9.20
RUN pip3 install pastml==1.9.20

# RUN pip3 install pastml==1.9.33


# Making PastML calculate_changes.py script available

RUN mkdir /opt/software/ && \
    cd /opt/software/ && \
    git clone https://github.com/evolbioinfo/pastml && \
    chmod +x /opt/software/pastml/pastml/utilities/calculate_changes.py && \
    ln -s /opt/software/pastml/pastml/utilities/calculate_changes.py /usr/local/bin/calculate_changes.py


# Installing plink version >= PLINK v1.90b6.17 64-bit (version PLINK v1.90b6.22 64-bit (16 Apr 2021) installed)
# adapted from done in: https://hub.docker.com/r/gelog/plink/dockerfile

RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20210416.zip && \
    mkdir /opt/plink_tmp/ && \
    unzip plink_linux_x86_64_20210416.zip -d /opt/plink_tmp/ && \
    rm plink_linux_x86_64_20210416.zip && \
    cp /opt/plink_tmp/plink /usr/local/bin && \
    rm -rf /opt/plink_tmp/


# Installing GCTA version >= version 1.93.2

RUN wget https://cnsgenomics.com/software/gcta/bin/gcta_1.93.2beta.zip && \
    mkdir /opt/gcta_tmp/ && \
    unzip gcta_1.93.2beta.zip -d /opt/gcta_tmp/ && \
    rm gcta_1.93.2beta.zip && \
    cp /opt/gcta_tmp/gcta_1.93.2beta/gcta64 /usr/local/bin && \
    rm -rf /opt/gcta_tmp/


# Installing pyseer Python module dependencies (as specified in the requirements.txt file)
# NOTE: freetype needed to install matplotlib
# NOTE: gfortran is needed to install glmnet_python
# NOTE: do not change versions from those in Pyseer requirements.txt as this may lead to compilation errors
# NOTE: cython needed to install pandas

# Installing pyseer version >= pyseer 1.3.7-dev
# NOTE: the command 'python3 -m pip install pyseer' did not work
# NOTE: do not use https://github.com/johnlees/pyseer directory as it is not the most updated
# NOTE: 'ln -s /usr/bin/python3 /usr/bin/python' used so that pyseer can run with header '#!/usr/bin/env python' using python3

RUN cd /opt/software/ && \
    git clone https://github.com/mgalardini/pyseer && \
    cd /opt/software/pyseer && \
    pip3 install -r requirements.txt && CYTHONIZE=1 pip3 install -e . && \
    chmod +x /opt/software/pyseer/pyseer-runner.py && \
    ln -s -f /opt/software/pyseer/pyseer-runner.py /usr/local/bin/pyseer && \
    chmod +x /opt/software/pyseer/scripts/phylogeny_distance.py && \
    ln -s /opt/software/pyseer/scripts/phylogeny_distance.py /usr/local/bin/phylogeny_distance.py && \
    ln -s /usr/bin/python3 /usr/bin/python


# Install Python dependencies for PowerBacGWAS scripts
# NOTE: the versions of scipy, numpy and pandas from Pyseer requirements.txt are different from the ones needed here
# NOTE: numpy==1.19.4 and scipy==1.6.0 are needed to run scripts simulate_binary_phenotype_*.py (ImportError: numpy.core.multiarray failed to import) 
# NOTE: scikit-learn==0.24.2 had to be re-installed to use newly installed numpy version
# NOTE: BioPython 1.78 installed to use annotate_nodes_newick.py
#  pytest \
#  pytest-cov \
#  pandas==1.1.4 \

RUN pip3 uninstall -y numpy scipy scikit-learn

RUN pip3 install \
  PyVCF==0.6.8 \
  numpy==1.19.4 \
  scipy==1.6.0 \
  scikit-learn==0.24.2 \
  biopython==1.78

# Install cyvcf2 from source (installation using pip3 failed)
# RUN pip3 install --no-binary cyvcf2 --no-cache-dir cyvcf2
# NOTE: the latest version of htslib had to be used, instead of the one on cyvcf2 GitHub directory, for cyvcf2 to be installed
# NOTE: the libssl-dev had to be installed to avoid compilation error
# NOTE: cyvcf2 must be installed after installing numpy

# RUN apt-get install libcrypto++6

RUN git clone --recursive https://github.com/brentp/cyvcf2 && \
  cd /opt/cyvcf2/ && \
  cd htslib && \
  autoheader && \
  autoconf && \
  ./configure --enable-libcurl && \
  make && \
  cd /opt/cyvcf2/ && \
  pip3 install -r requirements.txt && \
  CYTHONIZE=1 pip3 install -e .


# And finally, making all PowerBacGWAS scripts available and executable on image

WORKDIR /usr/local/bin
COPY ancestral_state_reconstruction.py prepare_gwas_runs.py	 prepare_vcf_file.py sample_casual_variants_from_roary.py simulate_phenotype_using_gcta.py ancestral_state_reconstruction_roary.py prepare_gwas_runs_roary.py process_gwas_runs.py sample_casual_variants_from_vcf.py subsample_plink_fam.py calculate_changes_per_region.py prepare_gwas_runs_subsampling.py roary_to_pastml_matrix.py simulate_binary_phenotype_roary.py vcf_to_pastml_matrix.py plink_to_pyseer_phenotype_file.py prepare_gwas_runs_subsampling_roary.py roary_to_plink_files.py simulate_binary_phenotype_vcf.py vcf_to_plink_files.py plot_gwas_runs.R plot_gwas_runs_subsampling.R annotate_nodes_newick.py subsample_vcf_samples.py subsample_pangenome_samples.py .
RUN chmod +x ancestral_state_reconstruction.py prepare_gwas_runs.py	 prepare_vcf_file.py sample_casual_variants_from_roary.py simulate_phenotype_using_gcta.py ancestral_state_reconstruction_roary.py prepare_gwas_runs_roary.py process_gwas_runs.py sample_casual_variants_from_vcf.py subsample_plink_fam.py calculate_changes_per_region.py prepare_gwas_runs_subsampling.py roary_to_pastml_matrix.py simulate_binary_phenotype_roary.py vcf_to_pastml_matrix.py plink_to_pyseer_phenotype_file.py prepare_gwas_runs_subsampling_roary.py roary_to_plink_files.py simulate_binary_phenotype_vcf.py vcf_to_plink_files.py plot_gwas_runs.R plot_gwas_runs_subsampling.R annotate_nodes_newick.py subsample_vcf_samples.py subsample_pangenome_samples.py

