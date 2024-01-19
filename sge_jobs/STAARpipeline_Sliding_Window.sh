#!/bin/bash

#$ -N SW_TM
#$ -l h_rt=144:00:00
#$ -t 1-573
#$ -l h_vmem=11G
#$ -tc 30
#$ -m n


#export R_LIBS_USER=$HOME/R-3.6.1-MKL
#echo $R_LIBS_USER

# specify the exact R and r libs I need
alias R="/opt/R/R-4.3.1/bin/R"
#export R_LIBS="/groups/garcia/users/kd2630/R/x86_64-pc-linux-gnu-library/4.3"
export R_LIBS="/groups/garcia/users/kd2630/R/x86_64-pc-linux-gnu-library/4.3:/opt/R/R-4.3.1/lib64/R/library"

# path info from .bash_profile
PATH=$PATH:$HOME/.local/bin:$HOME/bin
# add my own cmake
PATH=~/cmake-3.24.2/bin/:~/libxml2-2.9.9/bin/:$PATH
# adding my own
PATH=/opt/R/R-4.3.1/bin/:$PATH
# adding bcftools
PATH=~/bcftools/:$PATH
# adding htslib (mainly for a tabix that works)
PATH=~/htslib/:$PATH

export PATH

which R
which Rscript

basedir=~/noncoding_telo/STAAR/
savedir=/TOPMed_Euro/
dir_geno=${basedir}/${savedir}

Rscript --slave --no-restore --no-save ${basedir}/STAARpipeline-Tutorial-mod/STAARpipeline_Sliding_Window.r ${SGE_TASK_ID} ${dir_geno} > out"${SGE_TASK_ID}".Rout

