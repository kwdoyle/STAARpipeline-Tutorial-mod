#####################################################################
# Gene-centric analysis for noncoding rare variants of protein-coding 
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
afthresh <- as.numeric(commandArgs(TRUE)[3])

print(paste("Using MAF cutoff of", afthresh))

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
#mask_file <- "/home/STAAR/potential_masks/potential_noncoding_masks_chr_5_11_20_any_pval_eq_1.csv"

## output path
output_path <- paste0(basedir, "/Gene_Centric_Noncoding_Analysis/")
print(paste("Will create directory:", output_path))
dir.create(output_path)
## output file name
output_file_name <- "TOPMed_F5_LDL_Noncoding"
## input array id from batch file (Harvard FAS RC cluster)
#arrayid <- as.numeric(commandArgs(TRUE)[1])
# use chromosome as only input

###########################################################
#           Main Function 
###########################################################
## gene number in job
#gene_num_in_array <- 50 
#group.num.allchr <- ceiling(table(genes_info[,2])/gene_num_in_array)
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

genes_info_chr <- genes_info[genes_info[,2]==chr,]
sub_seq_num <- dim(genes_info_chr)[1]

#if(groupid < group.num)
#{
#  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):(groupid*gene_num_in_array)
#}else
#{
#  sub_seq_id <- ((groupid - 1)*gene_num_in_array + 1):sub_seq_num
#}
#
### exclude large noncoding masks
#jobid_exclude <- c(21,39,44,45,46,53,55,83,88,103,114,127,135,150,154,155,163,164,166,180,189,195,200,233,280,285,295,313,318,319,324,327,363,44,45,54)
#sub_seq_id_exclude <- c(1009,1929,182,214,270,626,741,894,83,51,611,385,771,493,671,702,238,297,388,352,13,303,600,170,554,207,724,755,1048,319,324,44,411,195,236,677)
#
#for(i in 1:length(jobid_exclude))
#{
#  if(arrayid==jobid_exclude[i])
#  {
#    sub_seq_id <- setdiff(sub_seq_id,sub_seq_id_exclude[i])
#  }
#}

# new way to exclude masks--by explicit gene name within a file
# don't do this
#mask_df <- read.csv(mask_file)
## remove masks from gene info table.
## this mask df has mask genes from all chrs, but should still be able to use it here like this.
#print("Removing mask genes")
#genes_info_chr = genes_info_chr[which(!genes_info_chr$hgnc_symbol %in% mask_df$hgnc_symbol), ]
## and re-set the row index value to loop over..
#sub_seq_num <- dim(genes_info_chr)[1]

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

# still not used in this script either..
#genes <- genes_info

results_noncoding <- c()
#for(kk in sub_seq_id)
# loop over each gene in the current chromosome
for(kk in 1:sub_seq_num)
{
  #print(kk)
  gene_name <- genes_info_chr[kk,1]
  print(gene_name)
  results <- Gene_Centric_Noncoding(chr=chr,gene_name=gene_name,genofile=genofile,obj_nullmodel=obj_nullmodel,
                                    rare_maf_cutoff=afthresh,rv_num_cutoff=2,
                                    QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
                                    Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
                                    Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name)
  
  results_noncoding <- append(results_noncoding,results)
}

save(results_noncoding,file=paste0(output_path,output_file_name,"_",chr,".Rdata"))

seqClose(genofile)

