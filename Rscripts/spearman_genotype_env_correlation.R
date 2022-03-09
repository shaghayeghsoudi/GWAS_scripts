#!/usr/bin/env Rscript

## This script estimates association between covariables and SNPs (without correction for population structure)

#setwd("/data/shaghayegh/H_annuus_WGS_2016/no_population_structure_geno.file_and_downstream_analysis_DP0_maf_0.05/run_file1")


args <- commandArgs(trailingOnly = TRUE) 
start1 <- 1
end1 <- 26


input <- read.table ("no_pop_str_part1", header = FALSE)
input$sub1 <- paste (input[,1],input[,2],sep = "__")

input_ord_data <- input[,3:29]
scafpos <- input[,1:2]

input_good <- data.matrix(input_ord_data)
scafpos_good <- as.matrix (scafpos)



results_out_p <- array (NA, c((nrow (input_good) * (end1 - start1 + 1)),7))
count <- 0

system.time(
#loop through focal snp
for (i in start1:end1){
	
#	loop through all SNPs
	for (j in 1:nrow (input_good)){
		count <- count + 1
		 results_out_p[count,7] <- cor.test(input_good[i,], input_good[j,], method = "spearman", use = "pairwise.complete.obs")$p.value
		results_out_p [count,1] <- scafpos_good[i,1]
		results_out_p [count,2] <- scafpos_good[i,2]
		#results_out_p [count,3] <- scafpos_good[i,3]
		results_out_p [count,4] <- scafpos_good[j,1]
		results_out_p [count,5] <- scafpos_good[j,2]
		#results_out_p [count,6] <- scafpos_good[j,3]
	
	}
}
)

outname_p <- paste ("output_p", start1,"_",end1,".txt",sep = "")
write.table (results_out_p, outname_p, col.names = F, row.names = F, quote = F)

