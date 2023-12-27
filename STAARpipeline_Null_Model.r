###########################################################
# fit STAAR null model
# Xihao Li, Zilin Li
# Initiate date: 11/04/2021
# Current date: 01/06/2023
###########################################################
rm(list=ls())
gc()

library(gdsfmt)
library(SeqArray)
library(SeqVarTools)
library(STAAR)
library(STAARpipeline)

###########################################################
#           User Input
###########################################################
## Phenotype file
#phenotype <- read.csv("/home/STAAR/TOPMed_Full_Cohort/topmed_input_model_data_cleaned.csv")
## (sparse) GRM file
#sgrm <- get(load("/path_to_the_file/sGRM.Rdata"))
## file directory for the output file 
#output_path <- "/home/STAAR/TOPMed_Full_Cohort/staar_null_model/"
## output file name
output_name <- "obj_nullmodel.Rdata"

phenotype_fl <- commandArgs(TRUE)[1]
output_path <- commandArgs(TRUE)[2]

phenotype <- read.csv(phenotype_fl)

print(paste("Using file", phenotype_fl))
print(paste("Saving to", output_path))

dir.create(output_path)

###########################################################
#           Main Function 
###########################################################
## fit null model
obj_nullmodel <- fit_nullmodel(is_IPF~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+as.factor(rs35705950),
                               data=phenotype,kins=NULL,use_sparse=NULL,kins_cutoff=0.022,id="sample.id",
                               groups=NULL,family=binomial(link="logit"),verbose=TRUE)
# run an attempt without using PCs to see if the inflation stays. This is to check if this method even adjusts for population effects at all
#obj_nullmodel <- fit_nullmodel(is_IPF~as.factor(rs35705950),
#                               data=phenotype,kins=NULL,use_sparse=NULL,kins_cutoff=0.022,id="sample.id",
#                               groups=NULL,family=binomial(link="logit"),verbose=TRUE)

# try using a linear mixed effect model and use ancestry as group
# note: what will do this is, I think, including the kins argument while keeping family=binomial
#obj_nullmodel <- fit_nullmodel(is_IPF~PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+as.factor(rs35705950),
#                               data=phenotype, kins=NULL, use_sparse=NULL, kins_cutoff=0.022, id="sample.id",
#                               groups="ancestry_cluster", family=gaussian(link="identity"), method.optim="AI", verbose=TRUE)

save(obj_nullmodel,file=paste0(output_path,output_name))

