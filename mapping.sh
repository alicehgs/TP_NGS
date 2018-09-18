#!/bin/bash
mkdir -p /mnt/data/variant_calling
cd /mnt/data/variant_calling

########################################################################################################################
# Requirements:
#	Java (version 8)
#	FastQC (version 0.11.7)
#	BWA-MEM (version 0.7.17-r1194-dirty)
#	SAMtools (version 1.9)
#	IGV (version 2.4.14)
#	GATK (version 3.3)
########################################################################################################################

java -version
fastqc -version
bwa #not installed
samtools
java -jar ${GATK} --help
java -jar ${PICARD}

##########################################################
## Download, extract and index the reference chromosome ##
##########################################################

# Download the reference Human chromosome (chromosome 20) from Ensembl
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed reference sequence (.fa.gz)
wget ftp://ftp.ensembl.org/pub/release-93/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz -O Homo_sapiens.Chr20.fa.gz

# Extract the reference chromosome
# Command: gunzip
# Input: compressed reference sequence (.fa.gz)
# Ouput: reference sequence (remove .gz file)
gunzip Homo_sapiens.Chr20.fa.gz

# Index the reference chromosome
# Command: bwa index
# Input: reference (.fa)
# Ouput: indexed reference (.fa.amb, .fa.ann, .fa.bwt, fa.pac, .fa.sa)
bwa index Homo_sapiens.Chr20.fa

######################################################
## Mapping of a family trio to the reference genome ##
######################################################

# The sequences are from an East Asian (Kinh Vietnamese) family forming a trio : daughter/mother/father
# Data available at http://www.internationalgenome.org/data-portal/sample/HG02024
# Daughter:
#       StudyId: SRP004063
#       SampleName: HG02024
#       Library: Pond-206419
#       ExperimentID: SRX001595
#       RunId: SRR822251
#       PlatformUnit: C1E0PACXX121221.6.tagged_373
#       InstrumentModel: Illumina HiSeq 2000
#       InsertSize: 160
# Mother:
#       StudyId: SRP004063
#       SampleName: HG02025
# Father:
#       SampleName: HG02026

#############################
## Mapping of the daughter ##
#############################

# Download paired sequencing reads for the daughter (SampleName: HG02024, RunId: SRR822251)
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed sequencing reads (.fastq.gz)
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_1.filt.fastq.gz

# Second pair
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/HG02024/sequence_read/SRR822251_2.filt.fastq.gz


# Map the paired sequencing reads against the reference Human chromosome 20
# Command: bwa mem
# Options: -M (Mark shorter split hits as secondary for GATK compatibility)
#          -t [number of CPU] (multi-threading)
# Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
# Ouput: alignment (.sam)
bwa mem -t 4 -M Homo_sapiens.Chr20.fa SRR822251_1.filt.fastq.gz SRR822251_2.filt.fastq.gz > SRR822251.sam

# (Optional)
# Compute summary statistics of the alignment
# Command: samtools flagstats
# Input: alignment (.sam)
# Ouput: text file (human and computer readable)
samtools flagstat SRR822251.sam > SRR822251.sam.flagstats # readable with cat

# Compress the alignment and filter unaligned reads
# Command: samtools view
# Options: -@ [number of CPU] (multi-threading) # here 4
#	   -S (input format is auto-detected)
# 	   -b (output BAM)
#	   -h (include header in output)
#          -f [flag] (include reads with all  of the FLAGs in INT present)
# 	        flag=3 for read paired & read mapped in proper pair, see
#	        https://broadinstitute.github.io/picard/explain-flags.html
# Input: alignment (.sam)
# Ouput: compressed alignment (.bam)
samtools view -Sbh -@ 4 -f 3 SRR822251.sam > SRR822251.bam


# Sort the alignment
# Command: samtools sort
# Input: compressed alignment (.bam)
# Ouput: sorted and compressed alignment (.bam)
samtools sort -@ 4 SRR822251.bam > SRR822251.sorted.bam

# Add Read group (cf https://gatkforums.broadinstitute.org/gatk/discussion/6472/read-groups)
# Command: gatk AddOrReplaceReadGroups
# Input: alignment (.bam) and read group (Read group identifier, DNA preparation library identifier, Platform, Platform Unit, Sample) # I=FILE en .bam
# Ouput: annotated alignment (.bam) # O=FILE en.bam
java -jar ${PICARD} AddOrReplaceReadGroups I=SRR822251.sorted.bam \
                                         O=daughter.bam \
                                         RGID=SRR822251 RGLB=Pond-206419 RGPL=illumina \
                                         RGPU=C1E0PACXX121221.6.tagged_373 RGSM=HG02024 RGPI=160


# (Optional)
# Compute statistics of the alignment
# Command: samtools-stats
# Input: alignment (.bam)
# Ouput: text file (human and computer readable)
samtools stats daughter.bam > daughter.bam.stats # bug ici 

# (Optional)
# Plot statistics of the alignment
# Command: plot-bamstats
# Input: statistics text file (output of samtools-stats)
# Ouput: plots (.png)
plot-bamstats -p ./plot daughter.bam.stats

# Index the alignment
# Command: samtools index
# Input: alignment (.bam)
# Ouput: indexed alignment (.bam.bai)
samtools index daughter.bam



###########################
## Mapping of the mother ##
###########################

# Variables definition
FTP_SEQ_FOLDER=ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3 # Ftp folder from 1000Genomes project
RUN_ID= # Read group identifier
SAMPLE_NAME= # Sample
INSTRUMENT_PLATFORM= # Platform/technology used to produce the read
LIBRARY_NAME= # DNA preparation library identifier
RUN_NAME= # Platform Unit
INSERT_SIZE= # Insert size

# Looking for the mother's SAMPLE_NAME
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.pedwget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_g1k.ped
head -n 1 20130606_g1k.ped > daughter.ped
grep "HG02024" 20130606_g1k.ped >> daughter.ped # in daughter.ped we can find the parents' SAMPLE_NAME. For the mother:
SAMPLE_NAME=HG02025

# Then, we download the file containing sequencing information.
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20130502.phase3.analysis.sequence.index

# We focus only on one specific sequencing run: 
RUN_ID=SRR361100
head -n 1 20130502.phase3.analysis.sequence.index > run_${RUN_ID}_mother.txt
grep "${SAMPLE_NAME}" 20130502.phase3.analysis.sequence.index | grep "${RUN_ID}" 20130502.phase3.analysis.sequence.index >> run_${RUN_ID}_mother.txt

# In run_SRR361100_mother.txt, we can find the following information:
INSTRUMENT_PLATFORM=ILLUMINA
LIBRARY_NAME=Catch-88584
RUN_NAME=BI.PE.110902_SL-HBC_0182_AFCD046MACXX.2.tagged_851.srf
INSERT_SIZE=96

 

# Download paired sequencing reads for the mother
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: compressed sequencing reads (.fastq.gz)
wget ${FTP_SEQ_FOLDER}/data/${SAMPLE_NAME}/sequence_read/${RUN_ID}_1.filt.fastq.gz
wget ${FTP_SEQ_FOLDER}/data/${SAMPLE_NAME}/sequence_read/${RUN_ID}_2.filt.fastq.gz

# Map, filter, and sort the paired sequencing reads of the mother against the reference genome
# Command: bwa mem && samtools view && samtools sort
# Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
# Ouput: sorted alignment (.bam)
bwa mem -t 4 -M Homo_sapiens.Chr20.fa ${RUN_ID}_1.filt.fastq.gz ${RUN_ID}_2.filt.fastq.gz | samtools view -Sbh -@ 4 -f 3 | samtools sort -@ 4 > ${RUN_ID}.sorted.bam 
# En vrai, utiliser beaucoup de pipes à la suite n'est pas idéal : c'est plus difficile d'identifier l'origine d'un problème. 

# Add Read group
# Command: gatk AddOrReplaceReadGroups
# Input: alignment (.bam) and read group
# Ouput: alignment (.bam)
java -jar ${PICARD} AddOrReplaceReadGroups I=${RUN_ID}.sorted.bam O=mother.bam \
                                         RGID=${RUN_ID} RGLB=${LIBRARY_NAME} RGPL=${INSTRUMENT_PLATFORM} \
                                         RGPU=${RUN_NAME} RGSM=${SAMPLE_NAME} RGPI=${INSERT_SIZE}

# Index the alignment
# Command: samtools index
# Input: alignment (.bam)
# Ouput: indexed alignment (.bam.bai)
samtools index mother.bam


###########################
## Mapping of the father ##
###########################

# Variables definition
SAMPLE_NAME=HG02026 # Sample

# Download index file containing sequencing runs information
# Command: wget
# Input: url (http:// or ftp://)
# Ouput: text file (.index)
wget ${FTP_SEQ_FOLDER}/20130502.phase3.analysis.sequence.index -O 20130502.phase3.index
# En fait, on l'a déjà !

# Filter paired exome sequencing runs related to father (HG02026)
# Command: grep && grep -v (-v: select non-matching lines
# Input: tab-separated values file (.index)
# Ouput: filtered comma-separated values file (.index)
grep ${SAMPLE_NAME} 20130502.phase3.index | grep "exome" | grep 'PAIRED' | grep -v 'Catch-88526' | grep -v 'Solexa' | grep -v 'from blood' | grep -v '_1.filt.fastq.gz' | grep -v '_2.filt.fastq.gz' | sed 's/\t/,/g' > father.index
# 'Catch-88526' est de très mauvaise qualité, c'est pourquoi on l'exclut


# File containing the list of alignments (each line is a .bam file)
# This file is necessary to merge multiple alignments into a single alignment.
# Command: touch
# Input: file name
# Ouput: empty file (.bamlist)
touch father.bamlist # creation of an empty file which name is "empty file"

# for each sequencing run (the first 10), align to the reference, sort, add read group and index

head -6 father.index | while IFS="," read FASTQ_FILE MD5 RUN_ID STUDY_ID STUDY_NAME CENTER_NAME SUBMISSION_ID SUBMISSION_DATE SAMPLE_ID SAMPLE_NAME POPULATION EXPERIMENT_ID INSTRUMENT_PLATFORM INSTRUMENT_MODEL LIBRARY_NAME RUN_NAME RUN_BLOCK_NAME INSERT_SIZE LIBRARY_LAYOUT PAIRED_FASTQ WITHDRAWN WITHDRAWN_DATE COMMENT READ_COUNT BASE_COUNT ANALYSIS_GROUP
do

    # Variables definition
    FASTQ_FILE_1=${FASTQ_FILE/.filt.fastq.gz/_1.filt.fastq.gz} # Path of the fasta file in the FTP folder
    FASTQ_FILE_2=${FASTQ_FILE/.filt.fastq.gz/_2.filt.fastq.gz} # Path of the fasta file in the FTP folder (pairing file)

    # Download paired sequencing reads for the father
    # Command: wget
    # Input: url (http:// or ftp://)
    # Ouput: compressed sequencing reads (.fastq.gz)
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_1} -O ${RUN_ID}_1.filt.fastq.gz
    wget ${FTP_SEQ_FOLDER}/${FASTQ_FILE_2} -O ${RUN_ID}_2.filt.fastq.gz

    # Map, filter, and sort the paired reads of the sequencing run against the reference genome
    # Command: bwa mem && samtools view && samtools sort
    # Input: indexed reference (.fa), and compressed sequencing reads (.fastq.gz)
    # Ouput: sorted alignment (.bam)
bwa mem -t 4 -M Homo_sapiens.Chr20.fa ${RUN_ID}_1.filt.fastq.gz ${RUN_ID}_2.filt.fastq.gz | samtools view -Sbh -@ 4 -f 3 | samtools sort -@ 4 > ${RUN_ID}.sorted.bam

    # Add Read group
    # Command: gatk AddOrReplaceReadGroups
    # Input: alignment (.bam) and read group
    # Ouput: alignment (.bam)
    java -jar ${PICARD} AddOrReplaceReadGroups I=${RUN_ID}.sorted.bam O=${RUN_ID}.sorted.RG.bam \
                                         RGID=${RUN_ID} RGLB=${LIBRARY_NAME} RGPL=${INSTRUMENT_PLATFORM} \
                                         RGPU=${RUN_NAME} RGSM=${SAMPLE_NAME} RGPI=${INSERT_SIZE}

    # Append the file name (.bam) to the list of alignments that will be merged
    echo ${RUN_ID}.sorted.RG.bam >> father.bamlist
done

# Merge the list of alignments into a single file
# Command: samtools merge
# Input: file containing the list of alignments (each line is a .bam file)
# Ouput: alignment (.bam)
samtools merge father.bam -b father.bamlist # father.bam = output file / father.bamlist = input.file

# Index the alignment
# Command: samtools index
# Input: alignment (.sam or .bam)
# Ouput: indexed alignment (.sam.bai or .bam.bai)
samtools index father.bam
