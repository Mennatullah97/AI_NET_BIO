

########################Cross Validation using 25000 iteration for DEGS identification###################################


library(limma)
datafile = "E:/dataset/discovery_dataset/unhealthy_behavior/obesity/obesity_discovery_metadataset_75%_symbol.csv"
phenotype_file = "E:/dataset/discovery_dataset/unhealthy_behavior/obesity/obesity_phenotypic.CSV"
pheno_name = "sample_obesity"
control_name = "control_"

metadata_cor = read.csv(datafile, header =T, row.names = 1L)
metadata_cor_pheno = read.csv(phenotype_file, header = T, row.names = 1L)


overlap_coef = matrix(data=NA, nrow=2500, ncol=1)
pheno_normal = metadata_cor[,metadata_cor_pheno$covariate==control_name]
pheno_ubnormal = metadata_cor[,metadata_cor_pheno$covariate==pheno_name]
no_normal= ncol(metadata_cor[, metadata_cor_pheno$covariate==control_name])%/%2
no_ubnormal = ncol(metadata_cor[,metadata_cor_pheno$covariate==pheno_name])%/%2

degs_table = data.frame()
DEG_list = list()

for (i in 1:2500){
  random_normal = pheno_normal[,c(sample(1:ncol(pheno_normal), no_normal, replace=F))]
  remaining_normal = pheno_normal[,-which(colnames(pheno_normal) %in% colnames(random_normal)) ]
  
  random_pheno = pheno_ubnormal[,c(sample(1:ncol(pheno_ubnormal), no_ubnormal, replace=F))]
  remaining_pheno = pheno_ubnormal[,-which(colnames(pheno_ubnormal) %in% colnames(random_pheno)) ]
  
  random = cbind(random_normal, random_pheno)
  remaining =cbind(remaining_normal,remaining_pheno)
  
  
  design_random = model.matrix(~0 +  as.factor(rep(c(control_name, pheno_name),
                                                   times = c(no_normal,no_ubnormal))))
  colnames(design_random) = c(control_name, pheno_name)
  fit_random = lmFit(random, design_random)
  
  ########add the phenotype name here 
  contrast_matrix_random = makeContrasts(behaviour_normal = sample_PTSD - control_PTSD, 
                                         levels = design_random)
  Control_fits_random = contrasts.fit(fit_random, contrast_matrix_random)
  Control_fits_random = eBayes(Control_fits_random)
  
  
  
  design_remaining = model.matrix(~0 +  as.factor(rep(c(control_name, pheno_name), 
                                                      times = c(ncol(remaining_normal),ncol(remaining_pheno)))))
  colnames(design_remaining) = c(control_name, pheno_name)
  fit_remaining = lmFit(remaining, design_remaining)
  
  ########add the pheno type name here 
  contrast_matrix_remaining = makeContrasts(behaviour_normal = sample_PTSD - control_PTSD,
                                            levels = design_remaining)
  Control_fits_remaining = contrasts.fit(fit_remaining, contrast_matrix_remaining)
  Control_fits_remaining = eBayes(Control_fits_remaining)
  
  random_no = length(rownames(topTable(Control_fits_random, number =nrow(random), adjust = 'BH')))
  remaining_no= length(rownames(topTable(Control_fits_remaining, number=nrow(remaining), adjust = 'BH')))
  
  Degs_random_tot = topTable(Control_fits_random, number=nrow(random), adjust = 'BH')
  Degs_random = Degs_random_tot[Degs_random_tot$adj.P.Val < 0.05,]
  
  Degs_remain_tot = topTable(Control_fits_remaining, number=nrow(remaining), adjust = 'BH')
  Degs_remain= Degs_remain_tot[Degs_remain_tot$adj.P.Val < 0.05,]
  
  random_remain_DEGs = c(unique(union(rownames(Degs_random), rownames(Degs_remain))))
  DEG_list = append(DEG_list, random_remain_DEGs)
  
  iteration_degs = rbind(Degs_random,Degs_remain)
  iteration_degs = iteration_degs[which(rownames(iteration_degs) %in% random_remain_DEGs),]
  degs_table = rbind(degs_table,iteration_degs)
  
  intersect = length(unique(union(rownames(Degs_random), rownames(Degs_remain))))
  overlap_coef[i,] =  intersect/(min(length(rownames(Degs_random)))+min(length(rownames(Degs_remain))))
  print(paste("just finished iteration:",i))
  }

setwd("E:/dataset/discovery_dataset/unhealthy_behavior/PTSD")
validation_DEGs = read.csv("PvC_discovery_DEGs.csv", header = T , row.names = 1L)

Degs_crossvalidated = unique(do.call(c,DEG_list))
Degs_crossvalidation = degs_table[which(row.names(degs_table) %in% Degs_crossvalidated),]

DEG_final = intersect(rownames(Degs_crossvalidation), rownames(validation_DEGs))
DEG_tabLe_final = Degs_crossvalidation[which(rownames(Degs_crossvalidation) %in% DEG_final),]

DEGs_up = DEG_tabLe_final[DEG_tabLe_final$logFC > 0 ,  ]
DEGs_down = DEG_tabLe_final[DEG_tabLe_final$logFC < 0 , ]
save(DEG_list, Degs_crossvalidation, DEG_tabLe_final,validation_DEGs, overlap_coef,DEGs_down, 
     DEGs_up,file = "degsdiscov.RDATA")
################################################File saving########################################################
load("degsdiscov.RDATA")
write.csv(DEG_tabLe_final, "cross_discov_DEGs.csv", sep = ",")
write.csv(DEGs_down, "downregulated_degdiscov_CROSS.csv", sep = ",")
write.csv(DEGs_up, "upregulated_degdiscov_CROSS.csv", sep = ",")




############checking the frequency of the overlap coefficient across the 2500 iteration  validation
setwd("E:/dataset/validation_dataset/unhealthy_dataset/PTSD")
load("degsdiscov.RDATA")
library(ggplot2)
library(grid)
library(gridExtra)

overlap_coef = as.data.frame(overlap_coef[!is.na(overlap_coef),])
colnames(overlap_coef) = c("overlap_cofficient")
histo = ggplot(as.data.frame(overlap_coef), aes(x=overlap_cofficient, fill=overlap_cofficient)) + 
  scale_color_manual(values=c("white"))+scale_fill_manual(values=c("green"))+ ylim(c(0,10))+
  xlim(c(0.5,1)) +geom_histogram(aes(y=after_stat(density)), bins=60,colour="black", fill=c("lightblue"))+geom_vline(xintercept=mean(overlap_coef$overlap_cofficient), color="pink",linetype=3,size=2)+
  geom_density(alpha=.2, fill= "brown") +theme_bw()+theme(legend.position = "none", panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),panel.border = element_blank(), panel.background = element_blank())
histo+ggtitle("The overlap coefficient frequence after the 2500 iterations \n Smoking discovery metadatset") 



################finding intersect with discovery, plotting correlation between disocvery and validation logFC 
library(reshape2)
library(corrplot)
setwd("E:/")

disocovery_DEGs = read.csv("cross_discov_DEGs.csv", 
                           header = T, row.names = 1L)
valid_degs = read.csv("cross_valid_DEGs_smoking.csv", header = T, row.names = 1L)
Meta_DEGs = merge(disocovery_DEGs, valid_degs, by = 0)
Meta_DEGs= as.data.frame(cbind(Meta_DEGs$Row.names, Meta_DEGs$logFC.x, Meta_DEGs$logFC.y))
rownamesdegs= Meta_DEGs[,1]
Meta_DEGs = Meta_DEGs[,-1]
Meta_DEGs= as.data.frame(lapply(Meta_DEGs, as.numeric))
rownames(Meta_DEGs)= rownamesdegs

colnames(Meta_DEGs) = c("logFC_Discovery", "logFC_validation")

write.csv(file = "common_degs_obes.csv", Meta_DEGs)
