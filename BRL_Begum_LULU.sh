#######################################
#     DNA-based arthropod diversity assessment in Amazonian iron mine lands show ecological succession towards undisturbed reference sites
# 		date: 10/09/19 
#		Rainy Season / Leray primer
# pipeline: Begum (https://github.com/shyamsg/Begum)
#######################################

## CHECK FASTQ FILES WITH FASTQC. Transfer data files to the correspondant pool. 
gunzip *gz
module load java
module load fastqc fastqc/v0.11.8
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- fastqc TOG_Q3MU_Leray_1_BrasilR_S4_L001_R1_001.fastq -threads 2
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- fastqc TOG_Q3MU_Leray_1_BrasilR_S4_L001_R2_001.fastq -threads 2
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- fastqc TOG_Q3MU_Leray_2_BrasilR_S5_L001_R1_001.fastq -threads 2
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- fastqc TOG_Q3MU_Leray_2_BrasilR_S5_L001_R2_001.fastq -threads 2
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- fastqc TOG_Q3MU_Leray_3_BrasilR_S6_L001_R1_001.fastq -threads 2
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- fastqc TOG_Q3MU_Leray_3_BrasilR_S6_L001_R2_001.fastq -threads 2
#Check the output files

## ADAPTER REMOVAL 
#Load AdapterRemoval
module load AdapterRemoval/v2.2.2   

#see help file
AdapterRemoval -h

# Remove adapters (adjust parameters accordingly)
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- AdapterRemoval --file1 TOG_Q3MU_Leray_1_BrasilR_S4_L001_R1_001.fastq.gz --file2 TOG_Q3MU_Leray_1_BrasilR_S4_L001_R2_001.fastq.gz --mm 0.05 --minlength 100 --shift 5 --basename ~/pool1_merged  --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse --threads 2
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- AdapterRemoval --file1 TOG_Q3MU_Leray_2_BrasilR_S5_L001_R1_001.fastq.gz --file2 TOG_Q3MU_Leray_2_BrasilR_S5_L001_R2_001.fastq.gz --mm 0.05 --minlength 100 --shift 5 --basename ~/pool2_merged  --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse --threads 2
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- AdapterRemoval --file1 TOG_Q3MU_Leray_3_BrasilR_S6_L001_R1_001.fastq.gz --file2 TOG_Q3MU_Leray_3_BrasilR_S6_L001_R2_001.fastq.gz --mm 0.05 --minlength 100 --shift 5 --basename ~/pool3_merged  --trimns --trimqualities --qualitybase 33 --minquality 28 --minalignmentlength 50 --collapse --threads 2

#Concatenate the merged collapsed files 
cat pool1_merged.collapsed pool1_merged.collapsed.truncated > pool1_merged.fastq
cat pool2_merged.collapsed pool2_merged.collapsed.truncated > pool2_merged.fastq
cat pool3_merged.collapsed pool3_merged.collapsed.truncated > pool3_merged.fastq


################## BEGUM ############################                                              
module load Begum

#Adding "tag" to the tags
awk '{print $1"\tTag"$2"\tTag"$3"\tPool"$4}' BRL_sampleInfo.txt  > BRL_sampleInfo.txt 
od -c leray_primers.txt | less

module load python/v2 

# Make all files ready
# In the working folder the following 4 files should exist:
			# BRL_poolInfo.txt 
			# leray_primers.txt
			# BRL_sampleInfo.txt #First column is sample, second TagF, third TagR, fourth, Pool.
			# BRL_tags.txt 

# Sort files with Begum                                                

# Navigate to project folder and there:
module load python/v2 

# With 2 primer mismatches
mkdir begum_2mismatches
xsbatch -c 1 --mem-per-cpu 8000 -- Begum sort -p leray_primers.txt -t BRL_tags.txt  -s BRL_sampleInfo.txt  -l BRL_poolInfo.txt -pm 2 -d begum_2mismatches -o begum_lerayRainy &> begum.log &

cd begum_2mismatches
cat begum_lerayRainy_pool1.summaryCounts 	#To check the number of reads
cat begum_lerayRainy_pool2.summaryCounts
cat begum_lerayRainy_pool3.summaryCounts 

mkdir filtered
Begum filter -i begum_lerayRainy -s ~/BRL_sampleInfo.txt -p 0.66 -m 10 -l 280 -d filtered -o filtered_0.66_10_280 

cd filtered
grep x filtered_0.66_10_280.fna
grep pos filtered_0.66_10_280.fna

#Check if your sequences are all the same length (roughly)
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen += length($0)}END{print seqlen}' filtered_0.66_10_280.fna > lengthseq_280.txt
	
	
## OTU CLUSTERING
#Convert working-file to required format.

#See help file
#From DAme: https://github.com/shyamsg/DAMe/blob/master/bin/convertToUSearch.py
python ~/convertToUSearch.py -h

#Convert
python ~/convertToUSearch.py -i filtered_0.66_10_280.fna -lmin 290 -lmax 330  

# Using sumaclust and sumatra for assess clustering Parameters
#From DAme: https://github.com/shyamsg/DAMe/blob/master/bin/assessClusteringParameters.py
module load sumaclust/v1.0.20  
module load sumatra/v1.0.20

xsbatch -c 1 --mem-per-cpu 8000 -- python ~/assessClusteringParameters.py -i FilteredReads.forsumaclust.fna -mint 0.90 -minR 0.95 -step 0.005 -t 4 -o LerayBR_clusterParams.pdf 

# Create OTU clusters with sumaclust
sumaclust -e FilteredReads.forsumaclust.fna -F OTUs_97_sumaclust.fna

## ASSING TAXONOMY
#Create a blast file
module load blast+/v2.6.0 
xsbatch -c 2 --time 180 --mem-per-cpu 11000 -- blastn -query table_BrazilRainyLeray_97_min280.txt.blast.txt -out  blast_BrazilR_LerayB_97_min280.output.txt -db nt -num_threads 2

#the  blast_BrazilR_LerayB_97_min280.output.txt file is the file to use with Megan. Open it using Megan to get the rough taxonomies.
	
#	open Megan / import from blast/ collapse tree/ choose ranking  (each level at the time) and safe each file: export - csv - readname_to_taxonname
#	each of these files can then be added to the main excel sheet 
#	Export as .csv for each result


##########################
#FOR  LULU:
#First produce a blastdatabase with the OTUs

module load blast+/v2.6.0 

#First produce a blastdatabase with the OTUs
makeblastdb -in table_BrazilRainyLeray_97.txt.blast.txt -parse_seqids -dbtype nucl


# Produce a match list, and blast the OTUs against the database
blastn -db table_BrazilRainyLeray_97.txt.blast.txt -outfmt '6 qseqid sseqid pident' -out match_list_BR_LerayB.txt -qcov_hsp_perc 80 -perc_identity 84 -query table_BrazilRainyLeray_97.txt.blast.txt

###############################

#Run the curation with LULU (done in R).

#Remove the sequences from the file and save it as filename_forLULU.txt.

R studio

library(devtools)

#install_github("tobiasgf/lulu")
library(lulu)
library(dplyr)

#Import dataset the files
#Have to remove anything else than the otus, sample and numbers from the otutable.
otutabL <- read.csv("table_BrazilRainyLeray_97_forLULU.txt",sep='\t',header=TRUE,as.is=TRUE, row.names = 1)
matchlist <- read.csv("match_list_BR_LerayB.txt", sep='\t',header=FALSE,as.is=TRUE)

curated_result <- lulu(otutabL, matchlist)

#No. of OTUs discarded: 
curated_result$discarded_count

#No. of OTUs retained: 
curated_result$curated_count

#Export new OTU table to working directory:
cleaned_OTUtab <- curated_result$curated_table
write.table(cleaned_OTUtab, file = "cleaned_OTUtab_BR_L.txt", sep = "\t")	

discarded_OTU <- curated_result$discarded_otus
write.table(discarded_OTU, file = "discarded_OTU_BR_L.txt", sep = "\t")	

#Create a final table, where you have OTUs after LULU + taxonomy assigned using MEGAN (remember to delete OTUs from the MEGAN files).
