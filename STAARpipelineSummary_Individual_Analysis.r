###########################################################
# Summarization and visualization of individual analysis
# results using STAARpipelineSummary
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 03/10/2022
###########################################################
rm(list=ls())
gc()

## load required packages
library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)
library(STAARpipelineSummary)

# for now, to operate on only some chromosomes, use my modified function
source("/home/STAAR/STAARpipeline-Tutorial/individ_analysis_summary_custom.R")
environment(Individual_Analysis_Results_Summary) <- asNamespace('STAARpipelineSummary')   
assignInNamespace("Individual_Analysis_Results_Summary", Individual_Analysis_Results_Summary, ns = "STAARpipelineSummary")

###########################################################
#           User Input
###########################################################
basedir <- commandArgs(TRUE)[1]

## Number of jobs for each chromosome
jobs_num <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/jobs_num.Rdata")))
## aGDS directory
agds_dir <- get(load(paste0(basedir, "/AssociationAnalysisPrestep/agds_dir.Rdata")))
## Known loci
#known_loci <- get(load("/path_to_the_file/TOPMed_F5_LDL_known_loci_genome_LD_pruning.Rdata"))
# For now, set this to NULL
known_loci <- NULL
## Null model
obj_nullmodel <- get(load(paste0(basedir, "/staar_null_model/obj_nullmodel.Rdata")))

## results path
input_path <- paste0(basedir, "/Individual_Variant_Analysis/")
output_path <- input_path
## results name
individual_results_name <- "TOPMed_F5_LDL_Individual_Analysis"

## QC_label
QC_label <- "annotation/filter"
## geno_missing_imputation
geno_missing_imputation <- "mean"
## method_cond
method_cond <- "optimal"
## alpha level
alpha <- 5E-08

###########################################################
#           Main Function 
###########################################################
Individual_Analysis_Results_Summary(agds_dir=agds_dir,jobs_num=jobs_num,input_path=input_path,output_path=output_path,
                                    individual_results_name=individual_results_name,
                                    obj_nullmodel=obj_nullmodel,known_loci=known_loci,
                                    method_cond=method_cond,
                                    QC_label=QC_label,geno_missing_imputation=geno_missing_imputation,
                                    alpha=alpha,manhattan_plot=TRUE,QQ_plot=TRUE)

