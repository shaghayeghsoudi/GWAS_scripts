#!/usr/bin/Rscript
library(corrplot)
#library(ape)

setwd("/data/home/shaghayegh/H_argophyllus_2018_HA412/baypass/maf_0.03/random_10K_files/inversion_removed")

#upload estimate of Omega
omega=as.matrix(read.table(file = "arg_HA412_noinv_random1_anacore_mat_omega.out"))
#pop.names=c("IA1A","MO1A","MO1W","ND1A","ND1W","SD1A","SD1W","SD2A","SD2W","SK1A","SK1W",
#            "IA1W","MK1","MK2","MK3","MK4","MK5","MK6","MK7","MK8","IA2A","IA2W","KS1A","KS1W","KS2A","KS2W","MB1W")



pop.names=c("ARG_01","ARG_02","ARG_03","ARG_04","ARG_05","ARG_06","ARG_07","ARG_08","ARG_09","ARG_10","ARG_11","ARG_12","ARG_13","ARG_14","ARG_15","ARG_16",
"ARG_17","ARG_18","ARG_19","ARG_20","ARG_21","ARG_22","ARG_23","ARG_24","ARG_25","ARG_26","ARG_27","ARG_28","ARG_29","ARG_30")

dimnames(omega)=list(pop.names,pop.names)


## plot covariance matrix ##
## pdf(file="BayPass_covariance_matrix_no_DP_filt_high_calling_rate_random1_run1_10K.pdf")
##
## corrplot(omega, method="circle",mar=c(2,1,2,2)+0.1,
##         main=expression("BayPass Covariance matrix"))

## dev.off()

#Compute and visualize the correlation matrix
cor.mat=cov2cor(omega)

pdf(file="BayPass_correlation_matrix_H_argophyllus_HA412_no.inv.random1_run1_10K.pdf")
corrplot(cor.mat,method="circle",mar=c(2,1,2,2)+0.1,
         main=expression("BayPass Correlation map based on"~hat(Omega)))

dev.off()

#Visualize the correlation matrix as hierarchical clustering tree
 bta14.tree=as.phylo(hclust(as.dist(1-cor.mat**2)))

pdf(file="BayPass_correlation_matrix_H_argophyllus_HA412_no.inv_hierarchical_clustering_tree.pdf")
 plot(bta14.tree,type="p",
      main=expression("Hier. clust. tree based on"~hat(Omega)~"("*d[ij]*"=1-"*rho[ij]*")"))

dev.off()

