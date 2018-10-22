#!/bin/bash

# ------- Polymorphism dataset 
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
gunzip ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
gunzip ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

#wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/* #download all the files in the folder



# ------- Filtering out non-coding polymorphism
wget ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.chr.gtf.gz
gunzip Homo_sapiens.GRCh38.94.chr.gtf.gz


gtf_to_bed.py -g Homo_sapiens.GRCh38.94.chr.gtf
# returns 
#   Homo_sapiens.GRCh38.94.chr.gtf_to_bed_errors.txt
#   Homo_sapiens.GRCh38.94.chr.bed
# now, data must be filtered

# merges overlapping BED/GFF/VCF entries into a single interval.
bedtools sort -i Homo_sapiens.GRCh38.94.chr.bed > Homo_sapiens.GRCh38.94.chr.bed.sorted

bedtools merge -c 4 -o distinct -i Homo_sapiens.GRCh38.94.chr.bed.sorted > Homo_sapiens.GRCh38.94.chr.bed.sorted.merge

# keeps all the polymorphisms in the exons
# -wb: write the original entry in B for each overlap
bedtools intersect -header -wb -a ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf -b Homo_sapiens.GRCh38.94.chr.bed.sorted > ALL.chr21.filtered.vcf