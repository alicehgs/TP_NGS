#!/bin/bash
installationDirectory="/mnt/data/tools"
mkdir -p ${installationDirectory}

########################################################################################################################
# FastQC
#   Version: 0.11.7
#   Licence: BSD, MIT
#   Author: Simon Andrews
#   URL: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#   Citation: Simon Andrews. FastQC: a quality control tool for high throughput sequence data. (2010).
########################################################################################################################

# Download and extract
cd ${installationDirectory}
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.7.zip
unzip fastqc_v0.11.7.zip

# Export the path of the executable (such that fastqc can be lunched from anywhere)
cd FastQC
chmod 755 fastqc
echo 'export PATH='$(pwd)':${PATH}' >> ~/.bashrc


########################################################################################################################
# BWA-MEM
#   Version: 0.7.17-r1194-dirty
#   Licence: GPLv3
#   Author: Heng Li (lh3@me.com)
#   Repository: https://github.com/lh3/bwa
#   Citation: Heng Li. Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM (2013).
#   DOI: https://arxiv.org/abs/1303.3997v2
########################################################################################################################

# Install libz necessary for compiling bwa
sudo apt-get install libz-dev

# Download (clone git repository)
cd ${installationDirectory}
git clone https://github.com/lh3/bwa.git

# Create the executable
cd bwa
make

# Export the path of the executable (such that bwa can be lunched from anywhere)
echo 'export PATH='$(pwd)':${PATH}' >> ~/.bashrc


########################################################################################################################
# SAMtools
#   Version: 1.9
#   Licence: BSD, MIT
#   Author: Heng Li
#   URL: http://www.htslib.org/
#   Repository: https://github.com/samtools/samtools
#   Citation: Heng Li, et al. The Sequence alignment/map (SAM) format and SAMtools. Bioinformatics, 25, 2078-9 (2009).
#   DOI:  https://doi.org/10.1093/bioinformatics/btp352
########################################################################################################################

# Download and extract
cd ${installationDirectory}
wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2
tar xvjf samtools-1.9.tar.bz2

# Create the executable
cd samtools-1.9
./configure --prefix=$(pwd)
make
make install

# Export the path of the executable (such that samtools can be lunched from anywhere)
echo 'export PATH=${PATH}:'$(pwd)'/bin' >> ~/.bashrc


########################################################################################################################
# Integrative Genomics Viewer (IGV)
#   Version: 2.4.14
#   Licence: MIT
#   Author: James T. Robinson
#   URL: http://software.broadinstitute.org/software/igv/
#   Repository: https://github.com/igvteam/igv/
#   Citation: James T. Robinson, et al. Integrative Genomics Viewer. Nature Biotechnology 29, 24–26 (2011).
#   DOI: https://doi.org/10.1038/nbt.1754
########################################################################################################################

# Download and extract
cd ${installationDirectory}
wget http://data.broadinstitute.org/igv/projects/downloads/2.4/IGV_2.4.14.zip
unzip IGV_2.4.14.zip

# Export the path of the executable (such that igv.sh can be lunched from anywhere)
cd IGV_2.4.14
mv igv.sh igv
chmod 755 igv
echo 'export PATH='$(pwd)':${PATH}' >> ~/.bashrc


########################################################################################################################
# Genome Analysis ToolKit (GATK)
#   Version: 3.8-1-0
#   Licence: BSD 3-Clause (https://github.com/broadinstitute/gatk/blob/master/LICENSE.TXT)
#   Author: Broad Institute, Inc (https://github.com/broadinstitute/gatk/blob/master/AUTHORS.TXT)
#   URL: https://software.broadinstitute.org/gatk/
#   Repository: https://github.com/broadinstitute/gatk
#   Citation: Aaron McKenna, et al.  The Genome Analysis Toolkit: a MapReduce framework for analyzing next-
#             generation DNA sequencing data. Genome research (2010).
#   DOI: https://dx.doi.org/10.1101%2Fgr.107524.110
########################################################################################################################


# Download and extract
cd ${installationDirectory}
wget https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive\&version=3.8-1-0-gf15c1c3ef -O gatk-3.8-1.tar.bz2
tar -xjf gatk-3.8-1.tar.bz2

# Export the path of the executable (such that gatk can be launched from anywhere)
echo 'export GATK='${installationDirectory}'/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar' >> ~/.bashrc

########################################################################################################################
# Picard Suite
#	Version: 2.0
#	URL: https://broadinstitute.github.io/picard/
#	Repository: https://github.com/broadinstitute/picard.git
########################################################################################################################

# Install java Jdk necessary to install Picard
sudo apt-get install openjdk-8-jdk

# Download and install
cd ${installationDirectory} 
git clone https://github.com/broadinstitute/picard.git
cd picard/
./gradlew shadowJar

# Export name of the directory
echo 'export PICARD='$(pwd)'/build/libs/picard.jar' >> ~/.bashrc

source ~/.bashrc
