
###############################finding the correlation between DEGs for smoking,Obesity,PTSD,MDD
library(reshape2)
library(corrplot)
degs_names= read.csv("E:/Paper1_all_files/common_degs_smoking.csv", row.names = 1L)
corelated_degs = degs_names[sign(degs_names$logFC_Discovery) == sign(degs_names$logFC_validation), ]
write.csv(corelated_degs,"common_degs_ptsds.csv")
expression_mat = read.csv("E:/PTSD_network/PTSD_DEGs_expression.csv")
corelated_degs = corelated_degs[which(rownames(corelated_degs) %in% expression_mat$X),]
cor_matrix = cor(corelated_degs, method = "spearman")


corrplot(
  cor_matrix,
  method = "shade",
  type = "lower",
  tl.cex = 0.8,
  tl.offset = 0.4,
  plotCI = "n",
  tl.col = "black",
  tl.srt = 10,
  addCoef.col = "white",
  insig = "p-value"
)
write.csv(corelated_degs,"E:/Paper1_all_files/common_degs_smoking.csv")



###################psych to create adjecency matrix
library("psych")

setwd("E:/Paper1_all_files/")
DEG_dataset = read.csv("PTSD_DEGs_expression.CSV", header = T, row.names = 1L)
corr_degs = read.csv("common_degs_ptsds.csv", header =T, row.names=1L)

DEG_dataset = DEG_dataset[which(rownames(DEG_dataset) %in% rownames(corr_degs)),]



Norm.interest.corr = corr.test( t(DEG_dataset), method="spearman", ci=F)

Norm.interest.corr$p[lower.tri( Norm.interest.corr$p,diag=-
                                  TRUE)]=NA
Pval.adj <- as.data.frame(as.table(Norm.interest.corr$p))


Norm.interest.corr$r [lower.tri( Norm.interest.corr$r,diag=TRUE)]=NA
Correlation =  as.data.frame(as.table(Norm.interest.corr$r))
Cor.table =  na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]
colnames(Cor.table) = c("gene1","gene2","cor","p.adj")
Cor.table.filt = Cor.table [Cor.table[,4] < 0.05,]
Cor.table.filt = Cor.table.filt [Cor.table.filt[,3] > 0.9, ]
write.table(Cor.table.filt[,1:3], "Cor_table_filter_mdd.txt",sep="\t", row.names=F, quote=F)







