#####################################################################
# Gene-centric analysis for noncoding rare variants of ncRNA
# genes using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/28/2022
#####################################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
chr <- as.numeric(commandArgs(TRUE)[1])
basedir <- commandArgs(TRUE)[2]

## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata"))) 
## Null model
obj_nullmodel <- get(load(paste0(basedir, "/staar_null_model/obj_nullmodel.Rdata")))

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "SNV"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## Annotation_dir
Annotation_dir <- "annotation/info/FunctionalAnnotation"
## Annotation channel
Annotation_name_catalog <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/Annotation_name_catalog.Rdata"))) 
# Or equivalently
# Annotation_name_catalog <- read.csv("/path_to_the_file/Annotation_name_catalog.csv")
## Use_annotation_weights
Use_annotation_weights <- TRUE
## Annotation name
Annotation_name <- c("CADD","LINSIGHT","FATHMM.XF","aPC.EpigeneticActive","aPC.EpigeneticRepressed","aPC.EpigeneticTranscription",
                     "aPC.Conservation","aPC.LocalDiversity","aPC.Mappability","aPC.TF","aPC.Protein")

# Masks to exclude
# (ncRNA masks are included in this noncoding mask file)
#mask_file <- "/home/STAAR/potential_masks/potential_noncoding_masks_chr_5_11_20_any_pval_eq_1.csv"

## output path
output_path <- paste0(basedir, "/Gene_Centric_ncRNA_Analysis/")
print(paste("Will create directory:", output_path))
dir.create(output_path)
## output file name
output_file_name <- "TOPMed_F5_LDL_ncRNA"
## input array id from batch file (Harvard FAS RC cluster)
#arrayid <- as.numeric(commandArgs(TRUE)[1])
# use chromosome as only input 

###########################################################
#           Main Function 
###########################################################
## gene number in job
#gene_num_in_array <- 100 
#group.num.allchr <- ceiling(table(ncRNA_gene[,1])/gene_num_in_array)
#sum(group.num.allchr)
#
#chr <- which.max(arrayid <= cumsum(group.num.allchr))
#group.num <- group.num.allchr[chr]
#
#if (chr == 1){
#  groupid <- arrayid
#}else{
#  groupid <- arrayid - cumsum(group.num.allchr)[chr-1]
#}

ncRNA_gene_chr <- ncRNA_gene[ncRNA_gene[,1]==chr,]
sub_seq_num <- dim(ncRNA_gene_chr)[1]

#if(groupid < group.num)
#{
#  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
#}else
#{
#  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
#}
#
### exclude large ncRNA masks
#if(arrayid==117)
#{
#  sub_seq_id <- setdiff(sub_seq_id,53)
#}
#
#if(arrayid==218)
#{
#  sub_seq_id <- setdiff(sub_seq_id,19)
#}
#
#if(arrayid==220)
#{
#  sub_seq_id <- setdiff(sub_seq_id,c(208,274))
#}
#
#if(arrayid==221)
#{
#  sub_seq_id <- setdiff(sub_seq_id,311)
#}
#
#if(arrayid==156)
#{
#  sub_seq_id <- setdiff(sub_seq_id,41)
#}
#
#if(arrayid==219)
#{
#  sub_seq_id <- setdiff(sub_seq_id,103)
#}

# new way to exclude masks--by explicit gene name within a file
# don't do this
#mask_df <- read.csv(mask_file)
## remove masks from gene info table.
## this mask df has mask genes from all chrs, but should still be able to use it here like this.
#print("Removing mask genes")
#ncRNA_gene_chr = ncRNA_gene_chr[which(!ncRNA_gene_chr$ncRNA %in% mask_df$hgnc_symbol), ]
## and re-set the row index value to loop over..
#sub_seq_num <- dim(ncRNA_gene_chr)[1]

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

# still not used
#genes <- genes_info

results_ncRNA <- c()
#for(kk in sub_seq_id)
# loop over each gene in the current chromosome instead
for(kk in 1:sub_seq_num)
{
  #print(kk)
  gene_name <- ncRNA_gene_chr[kk,2]
  print(gene_name)
  results <- c()
  results <- try(ncRNA(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                       rare_maf_cutoff=0.01,rv_num_cutoff=2,
                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name))
  
  results_ncRNA <- rbind(results_ncRNA,results)
}

save(results_ncRNA,file=paste0(output_path,output_file_name,"_",chr,".Rdata"))

seqClose(genofile)

