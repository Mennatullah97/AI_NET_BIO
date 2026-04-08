####single path for each study in validation dataset 

#############################################preprocessing for datasets#############################################
####Affyamatrix datasets
library(Biobase)
library(affy)
library(affyPLM)
library(oligo)
library(hugene10stv1cdf)
library(pd.hugene.1.1.st.v1)
library(pd.hugene.1.0.st.v1)
library(pd.hg.u133.plus.2)
library(pd.hg.u133a)
library(hgu133a.db)
library(hgu133plus2.db)
library(hgu133acdf)
library(hgu133plus2cdf)
setwd("E:/dataset/validation_dataset/unhealthy_dataset/Smoking/GSE18723_smoking_control_GPL96")
dataset = "E:/dataset/validation_dataset/unhealthy_dataset/Smoking/GSE18723_smoking_control_GPL96"
files = list.files(dataset, pattern = 'CEL')
filesoligo = list.celfiles(dataset, full.name = TRUE)
study ="GSE45803"

affydata= ReadAffy(celfile.path = dataset)
oligoreading = read.celfiles(filesoligo) 
rma_processing = rma(oligoreading) #or rma(affydata) if it wasn't transcriptcluster 

exp_raw = exprs(affydata)
exp = exprs(rma_processing)






#####reading and preprocessing Illumina datasets

library(limma)
library(Biobase)
library(illuminaio)
data =  read.table('E:/dataset/validation_dataset/unhealthy_dataset/MDD/GSE45603_MDD_control_illumina/GSE45603_Discovery_rawdata.txt'
                   ,row.names = 1L,sep = "\t",header = T)
targetinfo = read.csv("E:/dataset/validation_dataset/unhealthy_dataset/MDD/GSE45603_MDD_control_illumina/GSE45603_MDD_control.csv", 
                      header = T)
annot = read.csv("E:/dataset/validation_dataset/unhealthy_dataset/MDD/GSE45603_MDD_control_illumina/GPL10558_HumanHT-12_V4_0_R1_15002873_B.txt"
                 ,header = T, sep = "\t")

setwd("E:/dataset/validation_dataset/unhealthy_dataset/MDD/GSE45603_MDD_control_illumina")


controls = data[grep('', rownames(data))]
detectionpvalues = data[grepl("Detection.Pval", colnames(data))]
data = data[, !grepl("Detection.Pval", colnames(data))]


#creating Elist before normalization
project = new('EListRaw')
project@.Data[[1]] = 'illumina'
project@.Data[[2]] = targetinfo
project@.Data[[3]] = data
project@.Data[[4]] = NULL
project$E = data
project$targets = targetinfo
project$other$Detection = detectionpvalues


#background correction and quantile normaliztion equals to RMA 
exp_norm = neqc(project, offset = 16, robust = T)

#OR two step method 
exp_bgcorrection = backgroundCorrect(project)
exp_norm  = normalizeBetweenArrays(exp_bgcorrection)

###removing controls
Control = exp_norm$genes$Source=="ILMN_Controls"
NoSymbol = exp_norm$genes$Symbol == ""
isexpr <- rowSums(exp_norm$other$Detection <= 0.05) >= 3
exp_filtered <- exp_norm[isexpr, ]
exp_raw = as.data.frame(project$E)
exp = as.data.frame(exp_norm$E)
dim(exp_raw)
dim(exp)




#####reading and preprocessing Agilent datasets 

library(limma)

setwd('E:/dataset/validation_dataset/unhealthy_dataset/MDD/GSE19738_MDD_control_Agilent')
targetinfo = read.csv("E:/dataset/validation_dataset/unhealthy_dataset/MDD/GSE19738_MDD_control_Agilent/GSE19738_MDD_control.csv", 
                      header = T)
study ="GSE19738"
annot_AGI= read.csv( "GPL6480-9577.txt",header = TRUE,sep = '\t', stringsAsFactors = FALSE)


project =read.maimages(files = list.files(getwd(), pattern = ".txt"), source = 'agilent.median',
                       green.only = TRUE,
                       other.columns = 'gIsWellAboveBG')

project.bgcorrect = backgroundCorrect(project, method = 'normexp')

project.bgcorrect.norm = normalizeBetweenArrays(project.bgcorrect, method = 'quantile')

View(project.bgcorrect.norm$other)
Control = project.bgcorrect.norm$genes$ControlType==1L
NoSymbol = is.na(project.bgcorrect.norm$genes$external_gene_name)
IsExpr = rowSums(project.bgcorrect.norm$other$gIsWellAboveBG > 0) >= 4
exp_raw = project$E
exp1 =  project.bgcorrect.norm[!Control & IsExpr, ]
colnames(exp1$E)= targetinfo$sample
exp1 = exp1$E
dim(exp_raw)
dim(exp1)






#####plotting pre and post normalization 
setwd("E:/dataset/discovery_dataset/unhealthy_behavior/smoking/graphs/boxplots")
dev.new(width=3+ncol(exp_raw)/6, height=5)
par(mar=c(7,4,2,1))
title <- paste ("Pre-normalization probeset expression level for", study)
boxplot(log2(exp_raw[,-1]), boxwex=0.7, notch=T, main=title, outline=FALSE
        ,las=2,
        col= seq(1:ncol(exp_raw)))

title <- paste ("RMA-normalization probeset expression level for" ,study)
boxplot(exp1[,-1], boxwex=0.7, notch=T, main=title, outline=FALSE, las=2,
        col= seq(1:ncol(exp1)))

write.csv(file = paste(study, "rma_data.csv"), exp)

#######################################################Annotation################################################
####for Affymatrix 

library(hgu133plus2.db)
library(hgu133a.db)
library(hgu133acdf)
library(hgu133plus2cdf)
library(pd.hugene.1.0.st.v1)
library(pd.hugene.1.1.st.v1)
library(hugene10stv1cdf)
library(hugene10sttranscriptcluster.db)
library(hugene11sttranscriptcluster.db)

exp = read.csv("E:/dataset/validation_dataset/unhealthy_dataset/Smoking/GSE12585_smokers_affy_GPL96/GSE12585,_rma_annotated.csv",
                               header = T, row.names = 1L)
study = "GSE12585"
symbols <- hgu133aSYMBOL
probes <- mappedkeys(symbols)
probe_genes <- data.frame(symbols[probes])
rownames(probe_genes) <- probe_genes$probe_id
##2nd merging the genesymbols df to the metadata####
rdata = merge(probe_genes, exp,by=0)
aggregate_annotate = aggregate(rdata, by = list(rdata$symbol), mean)
exp.data.agg = data.frame(aggregate_annotate[, 5:ncol(aggregate_annotate)], 
                  row.names = aggregate_annotate$Group.1)
setwd("E:/dataset/validation_dataset/unhealthy_dataset/PTSD/")
write.csv(exp.data.agg, file = paste( study, "_rma_annotated.csv"), sep = ",")



######for Ilumina 

#assinging the sample names and annotating 
samplenames = targetinfo$sample
colnames(exp) = samplenames
#choosing gene symbols as the annoation symbols for the probset
annot = annot[,which(colnames(annot) %in% c('Probe_Id','Symbol'))]
annot = annot[which(annot$Probe_Id %in% rownames(exp)),]
annot = annot[match(rownames(exp), annot$Probe_Id),]
exp$Probe_Id = rownames(exp)
exp = exp[,c("Probe_Id",samplenames) ]
data2=merge(exp, annot,by="Probe_Id",all.x=T)
data2=data2[,-1]
data2=data2[ ! is.na(data2$Symbol), ]
x=duplicated(data2$Symbol)  
sum(x)
exp.data=data2[-dim(data2)[2]]
exp.data=apply(exp.data,2, as.numeric)
exp.data.agg= aggregate(exp.data, by=list(data2$Symbol),FUN=mean)
names(exp.data.agg)
rownames(exp.data.agg)=exp.data.agg$Group.1
exp.data.agg=exp.data.agg[- 1]

setwd('E:/dataset/validation_dataset/unhealthy_dataset/MDD/GSE39653_MDD_control_illumina/')
write.csv(exp.data.agg, file = paste( study,"_rma-annotated.csv"), sep= ",")




#####for Agilent

colnames(annot_AGI)[1] = 'AgilentID'
annot_AGI = annot_AGI[which(annot_AGI$AgilentID %in% exp1$genes$ProbeName),]
annot_AGI = annot_AGI[match(exp1$genes$ProbeName, annot_AGI$AgilentID),]
table(exp1$genes$ProbeName == annot_AGI$AgilentID) # check that annots are aligned


exp1$genes$AgilentID = annot_AGI$AgilentID
exp1$genes$GeneName= annot_AGI$GENE_SYMBOL

exp_data = as.data.frame(exp1$E)
samplenames = targetinfo$sample
colnames(exp_data) = samplenames
AgilentID= exp1$genes$AgilentID
exp_data =cbind(AgilentID, exp_data)

data2=merge(exp_data, annot_AGI,by="AgilentID",all.x=T)
data2=data2[,-1]
data2=data2[ ! is.na(data2$GENE_SYMBOL), ]
x=duplicated(data2$GENE_SYMBOL)  
sum(x)
exp.data=data2[-dim(data2)[2]]
exp.data=apply(exp.data,2, as.numeric)
exp.data.agg= aggregate(exp.data, by=list(data2$GENE_SYMBOL),FUN=mean)
names(exp.data.agg)
rownames(exp.data.agg)=exp.data.agg$Group.1
exp.data.agg=exp.data.agg[- 1]
exp.data.agg= exp.data.agg[,which(colnames(exp.data.agg) %in% samplenames)]


# setwd('E:/dataset/validation_dataset/unhealthy_dataset/MDD/GSE39653_MDD_control_illumina/')
write.csv(exp.data.agg, file = paste( study,"_rma-annotated.csv"), sep= ",")




 ###############################################exploratory analysis################################################
##########First: PCA analysis 
library(genefilter)
library(devtools)
library("factoextra")
library(ggplot2)
library(ggfortify)
library(scatterplot3d) 
targetinfo = read.csv("GSE18723_smoking_control.csv", header = T)
exp.data.agg= read.csv("GSE13850 _rma_annotated.csv", header = T, row.names = 1L)
study = "GSE13850"

mat_var = cbind(as.matrix(exp.data.agg), rowVars(as.matrix(exp.data.agg)))
mat_var1= mat_var[order(mat_var[,ncol(mat_var)],decreasing=T),]

######getting the top500 varieted genes 
mat_var1 = mat_var1[1:500,-ncol(mat_var1) ]


rpca = prcomp(t(mat_var1))

#####2D plot
colors_1= rep("gray",dim(targetinfo)[1])
colors_1[which(targetinfo$covariate=="sample_MDD")]="purple"
autoplot(rpca, data = targetinfo , color = "covariate", frame = T, frame.type = "norm")
#####3D plot


plot3d= scatterplot3d(rpca$x[,1:3], pch = 16, type="h",  color=colors_1)
legend("top", legend = c(unique(targetinfo$covariate)), col = c( "gray", "purple")
      , pch = 16,inset = -0.15 ,xpd = TRUE, horiz = TRUE)
############Second: Density plot
dat = stack(as.data.frame(exp.data.agg))
ggplot(dat, aes(x=values, fill=ind)) + geom_density(alpha=0.3)+ labs(title =paste("Density plot for" , study))





#############################################batch-effect removal using ComBat###################################
setwd("E:/dataset/validation_dataset/unhealthy_dataset/Smoking/")
ptsd1= "metadataset_MDD1.csv"
ptsd2 = "E-MTAB-5278 _rma_annotated.csv"
ptsd3= "E-MTAB-5332 _rma_annotated.csv"
ptsd4="GSE12585 _rma_annotated.csv"
ptsd5= "GSE13850 _rma_annotated.csv"
ptsd6 = "GSE18723 _rma_annotated.csv"
ptsd7= "GSE47415 _rma-annotated.csv"
bat1 = read.csv(ptsd1, header = T, row.names = 1L)
bat2 =read.csv(ptsd2, header= T, row.names = 1L)
bat3= read.csv(ptsd3, header = T, row.names = 1L)
bat4= read.csv(ptsd4, header = T, row.names = 1L)
bat5= read.csv(ptsd5, header = T, row.names = 1L)
bat6= read.csv(ptsd6, header = T, row.names = 1L)
bat7= read.csv(ptsd7, header = T, row.names = 1L)


#######prestep: merging & removing the 25th the dataset for each phenotype
library(genefilter)

mdata = merge(bat1, exp.data.agg, by=0)
mdata1 = merge(bat3,bat7, by=0)
mdata3 = merge(bat5, bat6, by=0)
mdata4 = merge(mdata, mdata1, by="Row.names")
mdata5 = merge(mdata3,mdata4, by = "Row.names")
rownames(mdata) = mdata$Row.names
mdata = mdata[,-1]
metadata_1 = merge(metadata_2, bat4, by=0)
rownames(metadata_1) = metadata_1$Row.names
metadata=metadata_1[,-1]
store_rownames = rownames(metadata)

# metadata= apply(metadata, 2,as.numeric)
# rowvars = rowVars(metadata)
# quartile_25_rma = quantile(rowvars, 0.25)
# rownames(metadata) = store_rownames
# metadata = metadata[rowvars > quartile_25_rma,]
write.csv(mdata, file = "metadataset_MDD.csv")




#####ComBat batch removal
library(Biobase)
library(sva)

setwd("E:/dataset/validation_dataset/unhealthy_dataset/MDD")
metadataset = read.csv("metadataset_MDD.csv", header = T, row.names = 1L)
pheno = read.csv("pheno_validation_MDD.csv", header = T, row.names = 1L)
metadata = as.matrix(metadataset)

x =setdiff(colnames(metadata),rownames(pheno))

phenotypicdata = data.frame(labelDescription = c( "covariate", "batch"),
                            row.names = c("covariate", "batch" ))

phenodata = new("AnnotatedDataFrame", data = pheno[,1:2], 
                varMetadata = phenotypicdata)

Unhealthy_eset = ExpressionSet(assayData = metadata, phenoData = phenodata)
model =  model.matrix(~ as.factor(pheno$covariate), data = pData(Unhealthy_eset))
val_data = ComBat(dat = exprs(Unhealthy_eset), batch = phenodata$batch, 
                                          mod = model, par.prior = TRUE, prior.plots = T)



####plotting the difference before and after batch-effect removal
library(Rfast)
library(ggpubr)

var_frma_only=  colVars(exprs(Unhealthy_eset))
var_frma_combat= colVars(val_data)
compare_var = rbind(cbind(var_frma_only, rep(c("rma"), ncol(val_data)))
                    ,cbind(var_frma_combat,rep(c("ComBat+rma"), ncol(val_data))))
colnames(compare_var) = c("variance", "group")
####boxplot 
var_boxplot= ggplot(as.data.frame(compare_var),
                    aes(x = reorder(factor(group), -as.numeric(variance)), 
                        y = as.numeric(variance))) + 
  geom_boxplot(fill = c("darkred", "lightgreen"), alpha = 0.5) +
  theme_bw() + theme(axis.title.y = element_text(size = 14),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5), 
                     axis.text.y = element_text(size = 14), 
                     plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
  labs(x = "Method", y = "Array Variance", title = "Batch-removal effect in \n MDD validation Meta-dataset using ComBat") +
  stat_compare_means(method = "wilcox.test", paired = TRUE, 
                     label = "p.format", label.y = 2.5, label.x.npc = "center") +
  scale_x_discrete(labels = c("rma only", "rma + combat"))
var_boxplot

####plotMDS
library(limma)

colors_1= rep("blue",dim(pheno)[1])
colors_1[which(pheno$batch==1)]="lightblue"
colors_1[which(pheno$batch== 2)]="green"
colors_1[which(pheno$batch== 3)]="darkred"
colors_1[which(pheno$batch== 4)]="red"
colors_1[which(pheno$batch== 5)]="lightgreen"
# colors_1[which(pheno$batch== 6)]="purple"



shapes =rep(16, dim(pheno)[1]) 
shapes[which(pheno$covariate== "sample_MDD")]= 17

par(mfrow = c(1, 2))
plotMDS(Unhealthy_eset, col= colors_1, pch = shapes)
title("MDD validation metadataset before batch-effect removal ")
legend("bottomright", legend = c(paste("batch",sort(unique(pheno$batch)))), col = c( "lightblue", "green", "darkred",
                                                                               "red","lightgreen","blue")
      , pch = c(16) ,inset = 0 ,xpd = TRUE, horiz = F, cex =0.8)
plotMDS(val_data, col=colors_1, pch = shapes)
title("MDD validation metadataset after ComBat-batch-effect removal ")
legend("topleft", legend = c(paste("batch",unique(pheno$batch))), col =c( "lightblue", "green", "darkred",
                                                                             "red","lightgreen","blue")
       , pch = c(16) ,inset = 0 ,xpd = TRUE, horiz = F, cex = 0.8)








#################################DEGs identification for each phenotype using limma#####################################
library(limma)
val_data = read.csv("E:/dataset/validation_dataset/unhealthy_dataset/PTSD/metadataset_ptsd.csv", header = T, 
                    row.names = 1L)
pheno = read.csv("E:/dataset/validation_dataset/unhealthy_dataset/PTSD/pheno_validation_PTSD.csv", header = T, 
                 row.names = 1L)
design_matrix = model.matrix(~0 + pheno$covariate, data = pheno)
colnames(design_matrix)= c("control_smoking","sample_smoking")

head(design_matrix)
contrast_matrix = makeContrasts(PvC= sample_smoking-control_smoking,
                                levels = design_matrix)
fit = lmFit(val_data, design_matrix)
fit2 = contrasts.fit(fit, contrast_matrix)
fit2 = eBayes(fit2)
results = decideTests(fit2)
DEGs = topTable(fit2, n=nrow(val_data), adjust = 'BH')
validation_DEGs= DEGs[DEGs$adj.P.Val < 0.05,]
write.csv(validation_DEGs, "SvC_validation_DEGs.csv")





##########################preparing the deg dataset from the metadata
setwd("E:/dataset/validation_dataset/unhealthy_dataset/Smoking/")
setwd("E:/PTSD_network/")
degs_names = read.csv("common_degs_PTSD.csv", header = TRUE, row.names = 1L)
disocovery_dat = read.csv("smoking_discovery_metadataset_75%_symbol.csv", header = TRUE, row.names = 1L)
validation_dat =  read.csv("metadataset_smoking.csv", header = TRUE, row.names = 1L)
pheno_val = read.csv("pheno_validation_smoking.csv", header = TRUE)
pheno_dis = read.csv("phenotype_smok_discov.csv", header = TRUE)
discovery = disocovery_dat[rownames(degs_names),]
validation = validation_dat[rownames(degs_names),]

DEG_dataset = as.data.frame(t(cbind(discovery,validation)))

state = c(pheno_dis$covariate,pheno_val$covariate)
DEG_dataset["state"] =state


write.csv(DEG_dataset, "smoking_DEGs_expression.csv")
exp_deg=read.csv("smoking_DEGs_expression.csv", header=TRUE)
exp_deg= as.data.frame(t(exp_deg))
write.csv(exp_deg, "smoking_DEGs_expression_1.csv")
















