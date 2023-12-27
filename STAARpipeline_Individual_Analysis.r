###########################################################
# Individual analysis using STAARpipeline
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 12/28/2022
###########################################################
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


## Number of jobs for each chromosome
jobs_num <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/jobs_num.Rdata")))
## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
## Null model
obj_nullmodel <- get(load(paste0(basedir, "/staar_null_model/obj_nullmodel.Rdata")))
paste0("Null model is: ", basedir, "/staar_null_model/obj_nullmodel.Rdata")

## QC_label
QC_label <- "annotation/filter"
## variant_type
variant_type <- "variant"
## geno_missing_imputation
geno_missing_imputation <- "mean"

## output path
#output_path <- "/home/STAAR/TOPMed_Full_Cohort/Individual_Variant_Analysis/"
output_path <- paste0(basedir, "/Individual_Variant_Analysis/")
print(paste("Will create directory:", output_path))
dir.create(output_path)
## output file name
output_file_name <- "TOPMed_F5_LDL_Individual_Analysis"
## input array id from batch file (Harvard FAS RC cluster)
# I think I can also just specify the chr to use instead of an array id
#arrayid <- as.numeric(commandArgs(TRUE)[1])

###########################################################
#           Main Function 
###########################################################
#chr <- which.max(arrayid <= cumsum(jobs_num$individual_analysis_num))
#group.num <- jobs_num$individual_analysis_num[chr]

#if (chr == 1){
#  groupid <- arrayid
#}else{
#  groupid <- arrayid - cumsum(jobs_num$individual_analysis_num)[chr-1]
#}

## aGDS file
agds.path <- agds_dir[chr]
genofile <- seqOpen(agds.path)

# this first part, groupid-1, is usually 0 anyway, so don't need it?
#start_loc <- (groupid-1)*10e6 + jobs_num$start_loc[chr]
start_loc <- jobs_num$start_loc[chr]
# always select the true end location.
# not sure why this min between the start + 10e6 and the true end is used..
# MAYBE they devised some way to run this in segments?
# but then the start locs would need to change for each "run" as well..
# I don't think that's done within Association_Analysis_PreStep..
# No, the start is always just the start that's in the gds. no modifications are done.
#end_loc <- start_loc + 10e6 - 1
#end_loc <- min(end_loc,jobs_num$end_loc[chr])
end_loc <- jobs_num$end_loc[chr]

a <- Sys.time()
results_individual_analysis <- c()
if(start_loc <= end_loc)
{
  results_individual_analysis <- Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,
                                                     genofile=genofile,obj_nullmodel=obj_nullmodel,mac_cutoff=20,
                                                     QC_label=QC_label,variant_type=variant_type,
                                                     geno_missing_imputation=geno_missing_imputation)
}
b <- Sys.time()
b - a

# note: switching this to save using the chr number might mess up other parts of the pipeline.
# will need to modify them later
#save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",arrayid,".Rdata"))
save(results_individual_analysis,file=paste0(output_path,output_file_name,"_",chr,".Rdata"))

seqClose(genofile)

