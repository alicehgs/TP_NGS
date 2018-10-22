#!/bin/bash

# If on Biosphere (IFB), use '/mnt/data/tools' as the installation directory
# INSTALL_DIR=/mnt/data/tools
INSTALL_DIR=$(pwd)

########################################################################################################################
# Bedtools 2
#   Version: 2.25.0
#   Licence: GNU Public License (Version 2).
#   Author: Quinlan AR and Hall IM
#   Repository: https://github.com/arq5x/bedtools2
#   Citation: Quinlan AR and Hall IM, 2010. BEDTools: a flexible suite of utilities for comparing genomic features.
#             Bioinformatics. 26, 6, pp. 841–842.
#   DOI: https://doi.org/10.1093/bioinformatics/btq033
########################################################################################################################

# Install libz necessary for compiling Bedtools
sudo apt install libz-dev

# Download (clone git repository)
cd ${INSTALL_DIR}
wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
tar -zxvf bedtools-2.25.0.tar.gz

# Create the executable
cd bedtools2
make

# Export the path of the executable (such that bwa can be lunched from anywhere)
echo 'export PATH='$(pwd)'/bin/:${PATH}' >> ~/.bashrc
source ~/.bashrc

########################################################################################################################
# VcfTools
#   Version: 0.1.17
#   Licence: GNU Lesser General Public License version 3.0
#   Author: Adam Auton, Petr Danecek and Anthony Marcketta
#   Repository: https://github.com/vcftools/vcftools
#   Citation: Petr Danecek et al, The Variant Call Format and VCFtools, Bioinformatics (2011)
#   DOI: https://doi.org/10.1093/bioinformatics/btr330
########################################################################################################################

# Install autoconf and pkg-config necessary for configuring VcfTools
sudo apt install autoconf
sudo apt install pkg-config

# Download (clone git repository)
cd ${INSTALL_DIR}
git clone https://github.com/vcftools/vcftools.git
cd vcftools

# Create the executable
./autogen.sh
./configure --prefix=$(pwd)
make
make install

# Export the path of the executable (such that bwa can be lunched from anywhere)
echo 'export PATH='$(pwd)'/bin:${PATH}' >> ~/.bashrc
echo 'export PERL5LIB='$(pwd)'/src/perl/' >> ~/.bashrc
source ~/.bashrc

########################################################################################################################
# Python custom scripts
########################################################################################################################

pip3 install BioPython --user
pip3 install matplotlib --user

cd ${INSTALL_DIR}

for PYTHON_SCRIPT in ./*.py
do
chmod a+x ${PYTHON_SCRIPT}
done

echo 'export PATH='$(pwd)':${PATH}' >> ~/.bashrc
source ~/.bashrc