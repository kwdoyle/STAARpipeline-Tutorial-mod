Individual_Analysis_Results_Summary <- function(agds_dir,jobs_num,input_path,output_path,individual_results_name,
                                                obj_nullmodel,known_loci=NULL,
                                                method_cond=c("optimal","naive"),
                                                QC_label="annotation/filter",geno_missing_imputation=c("mean","minor"),
                                                alpha=5E-09,manhattan_plot=FALSE,QQ_plot=FALSE){

	## evaluate choices
	method_cond <- match.arg(method_cond)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	## Summarize Individual Analysis Results
	results_individual_analysis_genome <- c()
#	num <- 0
	for(chr in 1:22)
	{
		print(chr)

#		if(chr > 1)
#		{
#			num <- num + jobs_num$individual_analysis_num[chr-1]
#		}

		results_individual_analysis_chr <- c()
		# only get rows where data exists (i.e., job num isn't 0)
#		jobs_num <- jobs_num[which(jobs_num$chr != 0), ]
		#for(kk in 1:jobs_num$individual_analysis_num[chr])
#		for(kk in 1:nrow(jobs_num))
#		{
		# pull files via chr number instead
		#chruse <- jobs_num$chr[kk]
		#print(kk)
		#print(chruse)
		#job_id <- kk + num
		results_individual_analysis <- try(get(load(paste0(input_path,individual_results_name,"_",chr,".Rdata"))))
		if (class(results_individual_analysis)[1] == "try-error") {
			next
		}

		results_individual_analysis_chr <- rbind(results_individual_analysis_chr,results_individual_analysis)
#		}
		results_individual_analysis_genome <- rbind(results_individual_analysis_genome,results_individual_analysis_chr)

		rm(results_individual_analysis_chr)
	}

	# save results
	save(results_individual_analysis_genome,file=paste0(output_path,"results_individual_analysis_genome.Rdata"))
	## Significant findings
	results_sig <- results_individual_analysis_genome[results_individual_analysis_genome$pvalue<alpha,]

	# save significant results
	save(results_sig,file=paste0(output_path,"results_individual_analysis_sig.Rdata"))
	write.csv(results_sig,paste0(output_path,"results_individual_analysis_sig.csv"))

	## manhattan plot
	if(manhattan_plot)
	{
		png(paste0(output_path,"manhattan_MAC_20.png"), width = 9, height = 6, units = 'in', res = 600)

		print(manhattan_plot(results_individual_analysis_genome$CHR, results_individual_analysis_genome$POS, results_individual_analysis_genome$pvalue, col = c("blue4", "orange3"),sig.level=alpha))

		dev.off()
	}

	## Q-Q plot
	observed <- sort(results_individual_analysis_genome$pvalue)
	lobs <- -(log10(observed))

	expected <- c(1:length(observed))
	lexp <- -(log10(expected / (length(expected)+1)))

	rm(results_individual_analysis_genome)
	gc()

	if(QQ_plot)
	{
		png(paste0(output_path,"qqplot_MAC_20.png"), width = 9, height = 9, units = 'in', res = 600)

		par(mar=c(5,6,4,4))
		plot(lexp,lobs,pch=20, cex=1, xlim = c(0, max(lexp)), ylim = c(0, max(lobs)),
		xlab = expression(Expected ~ ~-log[10](italic(p))), ylab = expression(Observed ~ ~-log[10](italic(p))),
		font.lab=2,cex.lab=2,cex.axis=2,font.axis=2)

		abline(0, 1, col="red",lwd=2)

		dev.off()
	}

	## Conditional Analysis
	if(length(known_loci)!=0)
	{
		results_sig_cond <- c()
		for(chr in 1:22)
		{
			if(sum(results_sig$CHR==chr)>=1)
			{
				results_sig_chr <- results_sig[results_sig$CHR==chr,]

				gds.path <- agds_dir[chr]
				genofile <- seqOpen(gds.path)

				results_sig_cond_chr <- Individual_Analysis_cond(chr=chr,individual_results=results_sig_chr,genofile,obj_nullmodel=obj_nullmodel,known_loci=known_loci,variant_type="variant", QC_label=QC_label, geno_missing_imputation=geno_missing_imputation, method_cond=method_cond)

				results_sig_cond <- rbind(results_sig_cond,results_sig_cond_chr)

				seqClose(genofile)
			}
		}
		save(results_sig_cond,file=paste0(output_path,"results_sig_cond.Rdata"))
		write.csv(results_sig_cond,paste0(output_path,"results_sig_cond.csv"))
	}
}
