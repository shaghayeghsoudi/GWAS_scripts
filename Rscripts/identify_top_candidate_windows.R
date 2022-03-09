###
#rm(list = ls())
## This script estimates outlier loci and produce manhattan plots for spearman correlation outputs


## load packages
library(qqman)
library(stringr)

### set working directory and load filal data-set containing loci with assigned bins to them
setwd("data/shaghayegh/H_annuus_WGS_2016/no_population_structure_geno.file_and_downstream_analysis_DP0_maf_0.05/plots")
all_good_chrom<- read.table (file= "var_out_Freebayes_H.annuss_2016_ALL.summary_chrom_good_noDP_filt_maf_0.05_filt_high_calling_rate_spearman", header = FALSE)

colnames(all_good_chrom)<-c("gene_id","chrom","pos","latitude_p","longitude_p","elevation_p","MAT_p", "MWMT_p","MCMT_p","TD_p","MAP_p","AHM_p","SHM_p","DD_0_p","DD5_p","DD_18_p","DD18_p","NFFD_p","bFFP_p","eFFP_p", "FFP_p","PAS_p","EMT_p","EXT_p","Eref_p","CMD_p","MAR_p","RH_p")

## assign types
test_type_env <- array ("0", ncol (all_good_chrom))
#env_type_re <- grep ("_re$", colnames (all_good_chrom))
env_type_pe <- grep ("_p$", colnames (all_good_chrom))

#test_type_env[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54)]<-"envir_re"
test_type_env[min(env_type_pe):(min(env_type_pe)+24)]<-"envir"

out_res_p <- NULL
test_names_env <- colnames (all_good_chrom)[test_type_env == "envir_pe"]
the_i <- which (test_type_env == "envir_pe") 
p_value <- 0.05


for (i in 1:length (the_i)){

	outliers <- all_good_chrom[,the_i[i]] <= p_value
	snps <- all_good_chrom[,the_i[i]] < 10000
	
	outliers_count <- tapply (outliers, list(as.character (all_good_chrom$window_20k)),sum, na.rm = T)
	outliers_count2 <- outliers_count[outliers_count >=1]

	snps_count <- tapply (snps, list(as.character (all_good_chrom$window_20k)),sum, na.rm = T)
	snps_count2 <- snps_count[outliers_count >=1]

	sub_good <- data.frame (names(snps_count2),snps_count2,outliers_count2)
	sub_good$test_name <- test_names_env[i]
	
	out_res_p <- rbind (out_res_p,sub_good)

}

colnames (out_res_p) <- c("window_20K", "snp_count","outlier_count", "test_name")
write.table(out_res_p,file = "out_res_all_spearman_p_threshold_150k_bin",col.names= TRUE, row.names = FALSE)

totsnp1 <- tapply (out_res_p$snp_count,list (out_res_p$test_name),sum)
totout1 <- tapply (out_res_p$outlier_count, list (out_res_p$test_name),sum)
expect1 <- data.frame (totout1 / totsnp1)
expect1$test_name <- row.names (expect1)
merg1 <- merge (out_res_p, expect1, by.x = "test_name", by.y = "test_name", all.x = T)
output <- merg1
colnames (output)[5] <- "expected"

## estimate qbiom and bbinom
output$p4 <- qbinom (0.9999,output$snp_count,output$expected)
output$p8 <- qbinom (0.99999999,output$snp_count,output$expected)
output$pscores <- -1 * log10(1 - pbinom (output$outlier_count,output$snp_count,output$expected))


outgood <- output[which (output$outlier_count > output$p4),]

the_lev <- unique (output$test_name)
write.table (output, "super_outliers_150k_bins_per_variable_H.annuus_2016_raw_ALL_spearman_p.txt", row.names = F, col.names = T, quote = F, sep = "\t")

#pdf ("Top_candidate_windowd20K_sunflowers")
par(mfrow=c(3,2))
for (i in 1:length (the_lev)){
	test1 <- output[output$test_name == the_lev[i],]
	test2 <- outgood[outgood$test_name == the_lev[i],]
	plot (test1$snp_count,test1$outlier_count, pch = 19, main = gsub("_pe","",the_lev[i]), cex = 0.6, xlab = "number of SNPs per 150k window",ylab = "number of outlier SNPs per window")
	points (test2$snp_count,test2$outlier_count, col = "blue", pch = 19, cex = 0.6)
	
}

#dev.off()


