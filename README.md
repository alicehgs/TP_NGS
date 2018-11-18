Project: Analyzing association patterns of deleterious polymorphisms in human populations
Date: September-Octobre 2018
Description: The project consists in reproducing the analysis of Sohail et al., Science (2017) to study association patterns of deleterious polymorphisms in diverse human populations using sequencing data from the 1000 Genomes project (more information about the 1000 Genomes project at http://www.internationalgenome.org/).   

Reference: 	Sohail et al., Science (2017)
		Negative selection in humans and fruit flies involves synergistic epistasis.
		Science 356.6337 (2017): 539-542.
		https://doi.org/10.1126/science.aah523874
Sohail et al. concludes that deleterious mutations are always associated with synergistic epistasis, i.e that an individual cannot accumulate multiple deleterious mutations. 


File list:
	installation.sh			
	installation_2.sh
	
	mapping.sh
	variant-calling.sh
	trio-analysis.sh
	
	analysis.sh
	--- gtf_to_bed.py
	--- vcf_coding_polymorphism.py
	--- vcf_meta_analysis.py

	meta_analysis.R


All the scripts can be downloaded at https://github.com/alicehgs/TP_NGS.git


 ---------------------------------------------------------------------------------------
|                                   1. Installation 					|
 ---------------------------------------------------------------------------------------
Files :	installation.sh	
	installation_2.sh

Run installation.sh and installation_2.sh to install all the tools that are needed for the project. 


 ---------------------------------------------------------------------------------------
|                                  2. Variant calling 					|
 ---------------------------------------------------------------------------------------
Files : mapping.sh
	variant-calling.sh
	trio-analysis.sh

This part of the project aimed at constructing a pipeline for variant calling, starting from raw sequencing data.

mapping.sh exemplifies how to map raw sequencing data on a reference. 


 -----------------  bwa index   -----------     
| Reference (.fa) | ---------> | Indexed   | ---
 -----------------             | reference |	|			         samtools   -----------------  samtools   -----------           gatk            -----------  samtools   -----------
				-----------     | bwa mem   ------------------	   view    | Compressed      |   sort    | Sorted    | AddOrReplaceReadGroups  | Annotated |   index   | Indexed   |
						---------> | Alignment (.sam) | ---------> | and filtered    | --------> | alignment | ----------------------> | alignment | --------> | alignment |   
						|	    ------------------		   | alignment (.bam)|    	 | (.bam)    |                         | (.bam)    |           | (.bam.bai)|
						|			 		    -----------------             -----------			        -----------             -----------
 ---------------------------                    | 														     |			     
| Sequencing reads (.fastq) | ------------------														     | samtools stats
 ---------------------------   																	     |
																				     ∨ 
																			     -----------------
																			    | Quality control |
																			     -----------------

variant-calling.sh corrects the mapping and perform variant calling.
The correction of the alignment is performed using the following Gatk tools:
	1. Correction for PCR artifacts: Gatk MarkDuplicates
	2. Local realignment around indels: Gatk RalignerTargetCretor + Gatk IndelRealigner
	3. Base quality recalibration: Gatk BaseRecalibrator + Gatk PrintReads
More information about the Gatk tools we found can be found at https://software.broadinstitute.org/gatk/documentation/tooldocs/current/.
Before running the variant calling step, we additionally removed all the reads mapped in intronic regions since we started from exome sequencing data. 

Run trio-analysis.sh to asses the efficiency of the variant calling pipeline using known family ties to identify unlikely variants (gatk PhaseByTransmission and VariantEval). The pipeline should be adjusted until the number of unlikely variants is negligible. 

--> This pipeline can be used to infer variants in thousands of individuals. The variant calling results are recorded in a .vcf file. 



 ---------------------------------------------------------------------------------------
|          3. Analysis of association patterns between deleterious mutations 		|
 ---------------------------------------------------------------------------------------
Files : analysis.sh
	--- gtf_to_bed.py
	--- vcf_coding_polymorphism.py
	--- vcf_meta_analysis.py
	
	meta_analysis.R

This part of the project aims at identifying the interaction pattern of deleterious mutations (defined as nonsense mutations in our analysis) in different human populations. 
We can distinguish three types of interactions between deleterious mutations:
	1. No interaction (or no epistasis)
	2. Attraction (or antagonistic epistasis)
	3. Repulsion (or synergistic epistasis)

We can define the subset of deleterious mutations included in the analysis using a cut-off for minor allele count. Sohail et al. Science (2017) only included singletons in their analysis of human populations.

Run analysis.sh to perform the analysis of Sohail et al., Science (2017) on a polymorphism dataset from the 1000 Genomes project. analysis.sh calls 3 .py scripts: gtf_to_bed.py, vcf_coding_polymorphism.py and vcf_meta_analysis.py. For more information about each .py script, see the comments provided in analysis.sh. 
meta_analysis.R is the R script used to represent graphically the .tsv files (previously converted into .csv files) output by vcf_meta_analysis.py.  


--> Results:
	- By contrast with the results obtained by Sohail et al., Science (2017), we found no interaction between deleterious mutations in most populations (about 70%). 
	- We inferred in about 30% of the populations either synergistic or antagonistic epistasis. However, the type of epistasis was found to be very sensitive to the choice of the cut-off: synergistic epistasis was inferred mostly for small cut-offs whereas antagonistic epistasis was inferred for high cut-offs. 
	- We also observed to some extent a geographic segregation of the populations for which we inferred either epistasis or antagonistic epistasis. 


