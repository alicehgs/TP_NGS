#!/bin/bash

cd /mnt/data/
mkdir 1000genomes
cd 1000genomes/

#### ------------------------------- PART I
#### 1.1-Polymorphism dataset

# Download the .vcf file containing the polymorphism dataset (aligned to the chromosome 20 of GRch38 human reference)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr20_GRCh38.genotypes.20170504.vcf.gz
gunzip ALL.chr20_GRCh38.genotypes.20170504.vcf.gz


#### 1.2-Filtering out non-coding polymorphism
# We want to select only deleterious polymorphisms. We assume that any polymorphism located inside a conding sequence is either deleterious or neutral, and therefore select all the exonic polymorphisms for our analysis.

# Download the .gtf file of chromosome 20 of the GRch38 human reference.
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.94.chr.gtf.gz

# We want to extract the start and stop position of all coding sequence from the above .gtf file. 
# gtf_to_bed.py extracts these information and outputs a .bed file. 
gtf_to_bed.py -g Homo_sapiens.GRCh38.94.chr.gtf


# Some exons can overlap. Therefore, we want to merge them. 
# First, the .bed file must be sorted.
bedtools sort -i Homo_sapiens.GRCh38.94.chr.bed > Homo_sapiens.GRCh38.94.chr.bed.sorted
# Then, we merge the overlaping exons.
bedtools merge -c 4 -o distinct -i Homo_sapiens.GRCh38.94.chr.bed.sorted > Homo_sapiens.GRCh38.94.chr.bed.sorted.merged
# Finally, we filter the .vcf file using the merged and sorted .bed file
bedtools intersect -header -wb -a ALL.chr20_GRCh38.genotypes.20170504.vcf -b Homo_sapiens.GRCh38.94.chr.bed.sorted.merged > Chr20_intersect.vcf


#### 1.3-Filtering by population or super-population
# We should study the polymorphism dataset by population or super-population to ensure that the condition of panmixie is fulfilled. 
# To filter the data by population or super-population, we first use extract_pop.py. 
# extract_pop.py takes in argument a metadata file that links each individual to its population and the name of the population we want to isolate. It returns a .txt file containing only the ID of the individuals that belong to the population of interest.
# Then, the .txt file can be used to filter the .vcf file, using vcftools.

# Download the metadata file
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel

# As an example:
POP="GBR"
extract_pop.py -p integrated_call_samples_v3.20130502.ALL.panel -k ${POP} # output: GBR.txt
vcftools --vcf Chr20_intersect.vcf --keep GBR.txt --recode --out Chr20_GBR.vcf


#### 1.3-Analysis
# See vcf_analysis.py


# As there are lots of different populations and super-population, we can make a loop to automate the filtering by population and the analysis. 

for MUT_TYPE in "NonSyn" "Syn" "Stop"
# for CLASS in Chr20_intersect.*.vcf
do 
    # Super population
    for SUPER_POP in "EUR" "AFR" "EAS" "AMR" "SAS"
    do
    extract_pop.py -p integrated_call_samples_v3.20130502.ALL.panel -k ${SUPER_POP}
    vcftools --vcf Chr20_intersect.${MUT_TYPE}.vcf --keep ${SUPER_POP}.txt --recode --out Chr20.${MUT_TYPE}_${SUPER_POP}
    vcf_analysis.py -v Chr20.${MUT_TYPE}_${SUPER_POP}.recode.vcf
    done

    # Populations
    for POP in "GBR" "FIN" "CHS" "PUR" "CDX" "JPT" "CLM" "IBS" "PEL" "PJL" "KHV" "LWK" "ACB" "GWD" "ESN" "BEB" "MSL" "MXL" "STU" "ITU" "CEU" "YRI" "CHB" "ASW" "TSI" "GIH"
    do
    extract_pop.py -p integrated_call_samples_v3.20130502.ALL.panel -k ${POP}
    vcftools --vcf Chr20_intersect.${MUT_TYPE}.vcf --keep ${POP}.txt --recode --out Chr20.${MUT_TYPE}_${POP}
    vcf_analysis.py -v Chr20.${MUT_TYPE}_${POP}.recode.vcf
    done

done





#### ------------------------------- PART II
#### 2.1-Polymorphism classification
# See vcf_coding_polymorphism.py
# vcf_coding_polymorphism split a .vcf into 3 .vcf files, each of them containing a certain type of polymorphism (Synonymous, Non-Synonymous or Stop). 
# To classify the polymorphisms... 

# Download the transcript sequence .FASTA file.
wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/cds/Homo_sapiens.GRCh38.cds.all.fa.gz
gunzip Homo_sapiens.GRCh38.cds.all.fa.gz


vcf_coding_polymorphism.py -f Homo_sapiens.GRCh38.cds.all.fa -v Chr20_intersect.vcf -g Homo_sapiens.GRCh38.94.chr.gtf 


#### 2.2-Meta-analysis
vcf_meta_analysis.py -n Chr20_intersect.NonSyn.vcf -o Chr20_intersect.Stop.vcf -s Chr20_intersect.Syn.vcf -p integrated_call_samples_v3.20130502.ALL.panel





