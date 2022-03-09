#!/usr/bin/Rscript 

library(dplyr)

setwd("/data/users/soudi/paper_sunflowers/top_candidates/annuus")
vars<-read.table(file = "var_out_annuus_2018_HA412_maf_0.03_all_covariates_spearman_window5k_attached_Pvalues.txt", header = FALSE)

colnames(vars)<-c("snp_id","chrom","pos","latitude_e","longitude_e","elevation_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","PAS_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e","window")

#colnames(vars)<-#c("snp_id","chrom","pos","latitude_e","longitude_e","elevation_e","MAT_e","MWMT_e","MCMT_e","TD_e","MAP_e","MSP_e","AHM_e","SHM_e","DD_0_e","DD5#_e","DD_18_e","DD18_e","NFFD_e","bFFP_e","eFFP_e","FFP_e","EMT_e","EXT_e","Eref_e","CMD_e","MAR_e","RH_e","window")


out_res<-read.table(file = "H_annuus_recombination_rate_bins.txt", header = TRUE)
loci_good<-vars[,c(2,3)]
chroms<-c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17") 

### assing SNPs to recombination rate bins
colnames(loci_good)<-c("chr","pos")
out_res_50k <- NULL
for (i in 1:17){
  jj<-out_res[out_res$chr==chroms[i],]       ## limit the bin file  to the chrom that is being looped through
  jj<-na.omit(jj)
  breaks <- as.vector(do.call(rbind,jj[,c(2,3)]))
  qq<-loci_good[loci_good$chr==i,] 
  qq<-na.omit(qq)
  ## limit the loci file  to the chrom that is being looped through
  ints_good <- ceiling(findInterval(qq[,2], sort(breaks))/2)
  
  window<-jj[ints_good,4]
  out_res_good <-data.frame(qq,window)
  out_res_50k <- rbind (out_res_50k,out_res_good)
  
}
colnames(out_res_50k)<-c("chr","pos","recomrate")
write.table(out_res_50k, file = "H.annuusHA412_recombination_rate_bins_assigned_to_windows.txt", col.names= TRUE, row.names = FALSE, quote = FALSE)


########
out_res_50k$quantile <- ntile(out_res_50k$recomrate, 5)  

win<-data.frame("window"=(vars[,5]))
out_res_50k<-cbind(out_res_50k,win)

write.table(out_res_50k, file = "out_res_5k_with_5quantile.table", col.names = TRUE, row.names = FALSE, quote= FALSE)

quantile<-data.frame("id"=(out_res_50k[,4]))

all_good_chrom<-cbind(vars,quantile)


### identify top candidates for each quantile

#colnames(all_good_chrom)<-c("gene_id","chrom","pos","latitude_p","longitude_p","elevation_p","MAT_p", "MWMT_p","MCMT_p","TD_p","MAP_p","AHM_p","SHM_p","DD_0_p","DD5_p","DD_18_p","DD18_p","NFFD_p","bFFP_p","eFFP_p", "FFP_p","PAS_p","EMT_p","EXT_p","Eref_p","CMD_p","MAR_p","RH_p")



## assign types
test_type_env <- array ("0", ncol (all_good_chrom))
#env_type_re <- grep ("_re$", colnames (all_good_chrom))
env_type <- grep ("_e$", colnames (all_good_chrom))

#test_type_env[c(4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54)]<-"envir_re"
test_type_env[min(env_type):(min(env_type)+25)]<-"envir"
#test_type_env[min(env_type):(min(env_type)+24)]<-"envir"

the_quantile <- 0.01

#bf <- 10
quantile1 <- quantile (all_good_chrom[, test_type_env == "envir"], the_quantile,na.rm = T)
test_names_env <- colnames (all_good_chrom)[test_type_env == "envir"]
the_i <- which (test_type_env == "envir") 
the_q<-c(1,2,3,4,5)
out_res <- NULL


for (i in 1:length (the_i)){
     
    for(j in 1:length(the_q)){
    all_good_chrom_foc<-all_good_chrom[,c(1,2,3,the_i[i],30,31),] 
    all_good_chrom_foc<-all_good_chrom_foc[all_good_chrom_foc$id==the_q[j],]
    #outliers <- all_good_chrom_foc[,the_i[i]] < quantile1
    outliers <- all_good_chrom_foc[,4] < quantile1
    snps <- all_good_chrom_foc[,4] < 10000
    
    outliers_count <- tapply (outliers, list(as.character (all_good_chrom_foc$window)),sum, na.rm = T)
    outliers_count2 <- outliers_count[outliers_count >=1]
    
    snps_count <- tapply (snps, list(as.character (all_good_chrom_foc$window)),sum, na.rm = T)
    snps_count2 <- snps_count[outliers_count >=1]
    
    sub_good <- data.frame (names(snps_count2),snps_count2,outliers_count2)
    sub_good$test_name <- test_names_env[i]
    sub_good$q_name <- the_q[j]
    
    out_res <- rbind (out_res,sub_good)
    
     }
     
}

colnames (out_res) <- c("windows", "snp_count","outlier_count", "test_name","q_name")
write.table(out_res,file = "out_res_H.annuusHA412_windows_withqID_0.95_quantile_Spearman_3vars.txt",col.names= TRUE, row.names = FALSE, quote= FALSE)


totsnp1 <- tapply (out_res$snp_count,list (out_res$test_name),sum)
totout1 <- tapply (out_res$outlier_count, list (out_res$test_name),sum)
expect1 <- data.frame (totout1 / totsnp1)
expect1$test_name <- row.names (expect1)
merg1 <- merge (out_res, expect1, by.x = "test_name", by.y = "test_name", all.x = T)
output <- merg1
colnames (output)[6] <- "expected"

output$p4 <- qbinom (0.9999, output$snp_count, output$expected)
output$p8 <- qbinom (0.99999999, output$snp_count, output$expected)
#output$pscores <- -1 * log10(1 - pbinom (output$outlier_count, output$snp_count, (output$outlier_count/output$snp_count)))
outgood <- output[which (output$outlier_count > output$p4),]

the_lev <- unique (output$test_name)

recombination_bin<-c("0-20%","20-40%","40-60%","60-80%","80-100%")


for (i in 1:length (the_lev)){
  
  pdf (file = paste("Vertical_test_super_outliers_Spearman_window5K_NONrecomadjusted_variable_H.annuusHA412_",the_lev[i],".pdf", sep = ""),width = 4, height= 14)
  par(mfrow=c(5,1))
  par(mar = c(5.1,4.8,4.1,2.1))
  for(j in 1:length(the_q)){
    output_focal<-output[output$test_name==the_lev[i],]
    output_focal<-output_focal[output_focal$q_name==the_q[j],]  
    outgood_focal<-outgood[outgood$q_name==the_q[j],]
    test1 <- output_focal[output_focal$test_name == the_lev[i],]
    test2 <- outgood_focal[outgood_focal$test_name == the_lev[i],]
    plot (test1$snp_count,test1$outlier_count, pch = 1, main = recombination_bin[j], cex = 0.7, xlab = "number of SNPs per window",ylab = "number of outliers per window", cex.lab=1.5, cex.axis=1.5, ylim=c(0,80), xlim = c(0,150))
    points (test2$snp_count,test2$outlier_count, col = "red", pch = 1, cex = 0.7)
  }
  
  write.table (output, "super_outliers_window5K_withq_id_quantileper_variable_H.annuusHA412_raw_3vars_spearman.txt", row.names = F, col.names = T, quote = F, sep = "\t")
  write.table (outgood, "super_outliers_OUTGOOD_window5K_withq_id_quantileper_variable_H.annuusHA412_raw_3vars_spearman.txt", row.names = F, col.names = T, quote = F, sep = "\t")
  
  dev.off()
}


#### alternative_with_vertical_shape
#for (i in 1:length (the_lev)){
#  
#  pdf (file = paste("super_outliers_BayPass_window5K_NONrecomadjusted_variable_H.annuusHA412_",the_lev[i],".pdf"),width = 3, height= 12)
#  par(mfrow=c(5,1))
#  par(mar = c(5.1,4.8,4.1,2.1))
#  for(j in 1:length(the_q)){
#    output_focal<-output[output$test_name==the_lev[i],]
#    output_focal<-output_focal[output_focal$q_name==the_q[j],]  
#    outgood_focal<-outgood[outgood$q_name==the_q[j],]
#    test1 <- output_focal[output_focal$test_name == the_lev[i],]
#   test2 <- outgood_focal[outgood_focal$test_name == the_lev[i],]
#    plot (test1$snp_count,test1$outlier_count, pch = 1, main = recombination_bin[j], cex = 0.8, xlab = "number of SNPs per 5K window",ylab = "number of outliers per 5K window", cex.lab=1.2, cex.axis=1.5, ylim=c(0,80), xlim = c(0,150))
#    points (test2$snp_count,test2$outlier_count, col = "red", pch = 1, cex = 0.8)
#  }
#  
#  write.table (output, "super_outliers_window5K_withq_id_quantileper_variable_H.annuusHA412_raw_3vars_BayPass.txt", row.names = F, col.names = T, quote = F, sep = "\t")
#  write.table (outgood, "super_outliers_OUTGOOD_window5K_withq_id_quantileper_variable_H.annuusHA412_raw_3vars_BayPass.txt", row.names = F, col.names = T, quote = F, sep = "\t")
#  
#  dev.off()
#}



############################################################################
### assing different thresholds for each recombination rate bin ###
############################################################################
the_lev <- unique (out_res$test_name)
for (i in 1:length (the_lev)){
  
  pdf (file = paste("super_outliers_Spearman_window5K_recomAdjusted_variable_H.annuusHA412_",the_lev[i],".pdf", sep = ""),width = 4, height= 14)
  par(mfrow=c(5,1))
  par(mar = c(5.1,4.8,4.1,2.1))
  
  for(j in 1:length(the_q)){
    
  out_res_focal<-out_res[out_res$test_name==the_lev[i],]
  out_res_focal<-out_res_focal[out_res_focal$q_name==the_q[j],]
  totsnp1 <- tapply (out_res_focal$snp_count,list (out_res_focal$test_name),sum)
  totout1 <- tapply (out_res_focal$outlier_count, list (out_res_focal$test_name),sum)
  expect1 <- data.frame (totout1 / totsnp1)
  expect1$test_name <- row.names (expect1)
  merg1 <- merge (out_res_focal, expect1, by.x = "test_name", by.y = "test_name", all.x = T)
  output <- merg1
  colnames (output)[6] <- "expected"
  
  output$p4 <- qbinom (0.9999, output$snp_count, output$expected)
  output$p8 <- qbinom (0.99999999, output$snp_count, output$expected)
  #output$pscores <- -1 * log10(1 - pbinom (output$outlier_count, output$snp_count, (output$outlier_count/output$snp_count)))
  outgood <- output[which (output$outlier_count > output$p4),]
 
  #the_lev <- unique (output$test_name)
  
  
  #for(j in 1:length(the_q)){
  
  recombination_bin<-c("0-20%","20-40%","40-60%","60-80%","80-100%")
  test1 <- output[output$test_name == the_lev[i],]
  test2 <- outgood[outgood$test_name == the_lev[i],]
  plot (test1$snp_count,test1$outlier_count, pch = 1, main = recombination_bin[j], cex = 0.7, xlab = "number of SNPs per window",ylab = "number of outliers per window", cex.lab=1.5, cex.axis=1.5, ylim=c(0,80), xlim = c(0,150))
  points (test2$snp_count,test2$outlier_count, col = "red", pch = 1, cex = 0.7)
  write.table (output, file = paste("super_outliers_window5K_withq_id_recomadjusted_quantileper_variable_H.annuusHA412_raw_ALL_",the_lev[i],"_recomrate_",the_q[j],".txt", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
  write.table (outgood, file = paste("super_outliers_OutGoodwindow5K_withq_id_recomadjusted_quantileper_variable_H.annuusHA412_raw_ALL_",the_lev[i],"_recomrate_",the_q[j],".txt", sep = ""), row.names = F, col.names = T, quote = F, sep = "\t")
  
  } 
  
dev.off()

}


