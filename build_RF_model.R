

library(randomForest)
library(ROCR)
library(genefilter)
library(Hmisc)
library(randomForestExplainer)
library(rfPermute)
library(limma)
library(PRROC)


#Set working directory and filenames for Input/output
setwd("E:/PTSD_network/")


degs_names = read.table("Cor_table_filter_ptsd.txt", header = TRUE)
disocovery_dat = read.csv("PTSD_discovery_metadataset_75%_symbol.csv", header = TRUE, row.names = 1L)
validation_dat =  read.csv("metadataset_ptsd.csv", header = TRUE, row.names = 1L)
discov_pheno =read.csv("PTSD_phenotypic.csv", header = TRUE, row.names = 1L)
valid_pheno = read.csv("pheno_validation_PTSD.csv", header = TRUE, row.names = 1L)


degs=c(degs_names[,1],degs_names[,2])
degs = unique(degs)
data_import=cbind(disocovery_dat[degs,],validation_dat[degs,])
clindata= rbind(discov_pheno[,1:2], valid_pheno)

predictor_data = as.data.frame(t(data_import))
target= clindata[,1]
target[target==0]="sample_PTSD"
target[target==1]="control_PTSD"
target=as.factor(target)

#Run RandomForests
tmp = as.vector(table(target))
num_classes = length(tmp)
min_size = tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes = rep(min_size,num_classes)
rf_output=randomForest(x=predictor_data, y=target, importance = TRUE, ntree = 50001,sampsize = sampsizes, proximity=TRUE, na.action = na.omit)



#Save RF classifier with save()
save(rf_output, file="RF_model_Ptsd")
load("RF_model_PTSD")

RF_predictor_names=rownames(rf_output$importance)



min_depth_frame = min_depth_distribution(rf_output)
save(min_depth_frame, file = "min_depth_frame.rda")
plot_min_depth_distribution(min_depth_frame, mean_sample= "top_trees", k=20 ,
                            main = "Distrubutuion of PTSD Modules minimum depth mean" )

importance_frame = measure_importance(rf_output)
save(importance_frame, file = "importance_frame.rda")
load("importance_frame.rda")
rf_degs = importance_frame[ importance_frame$p_value < 0.05,]
deg_expresseion = read.csv("E:/PTSD_network/PTSDD_DEGs_expression.csv", row.names = 1L)
deg_expresseion = deg_expresseion[which(rownames(deg_expresseion) %in% rf_degs$variable),]
write.csv(deg_expresseion,"PTSD_RF_DEGS_expression.csv")
common_degs = read.csv("common_degs_PTSD.csv", row.names = 1L, header = TRUE)
common_degs = common_degs[which(rownames(common_degs) %in% rownames(deg_expresseion)),]
write.csv(common_degs, "common_degs_PTSD1.csv")


####ploting the function
load("E:/MDD_network/importance_frame.rda")
load("E:/MDD_network/min_depth_frame.rda")
plotProximity(rf_output, legend.type = "label", point.size = NULL, circle.size = NULL, group.alpha = 0.5)
plotProximity(rf_output, point.size = NULL, circle.size = NULL, group.alpha = 0.5)
plotProximity(rf_output, group.type = "hull")
plotProximity(rf_output, group.type = "contour")


plot(rf_output, main = "Learning curve of RF model for PTSD")
legend("topright", c("error for samples", "OOB error", "error for controls"), lty = c(1,1,1), col = c("green", "black", "red"))


plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes",
                          main = "Multi-way importance plot for PTSD DEGs",min_no_of_trees = 50, no_of_labels = 40)
plot_multi_way_importance(importance_frame, no_of_labels = 30,main = "Multi-way importance plot for PTSD DEGs",
                          x_measure = "accuracy_decrease", y_measure = "gini_decrease", size_measure = "p_value")

summary(rf_output)
plotMDS(rf_output$proximity)

#######################################
#Load saved model - this will save time if re-running and you want to skip RF run

#Get importance measures
rf_importances=importance(rf_output, scale=FALSE)

#Determine performance statistics
confusion=rf_output$confusion
sensitivity=(confusion[2,2]/(confusion[2,2]+confusion[2,1]))*100
specificity=(confusion[1,1]/(confusion[1,1]+confusion[1,2]))*100
overall_error=rf_output$err.rate[length(rf_output$err.rate[,1]),1]*100
class1_error=paste(rownames(confusion)[1]," error rate= ",confusion[1,3], sep="")
class2_error=paste(rownames(confusion)[2]," error rate= ",confusion[2,3], sep="")
overall_accuracy=100-overall_error

#Prepare stats for output to file
sens_out=paste("sensitivity=",sensitivity, sep="")
spec_out=paste("specificity=",specificity, sep="")
err_out=paste("overall error rate=",overall_error,sep="")
acc_out=paste("overall accuracy=",overall_accuracy,sep="")
misclass_1=paste(confusion[1,2], rownames(confusion)[1],"misclassified as", colnames(confusion)[2], sep=" ")
misclass_2=paste(confusion[2,1], rownames(confusion)[2],"misclassified as", colnames(confusion)[1], sep=" ")

#Prepare confusion table for writing to file
confusion_out=confusion[1:2,1:2]
confusion_out=cbind(rownames(confusion_out), confusion_out)

#Print results to file
outfile ="importance_MDD.txt"
write.table(rf_importances[,4],outfile, sep="\t", quote=FALSE, col.names=FALSE)
write("confusion table", file=outfile, append=TRUE)
write.table(confusion_out,file=outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE, append=TRUE)
write(c(sens_out,spec_out,acc_out,err_out,class1_error,class2_error,misclass_1,misclass_2), file=outfile, append=TRUE)

#Produce MDS plot
target_labels=as.vector(target)
target_labels[target_labels=="sample_MDD"]="S"
target_labels[target_labels=="control_MDD"]="C"
MDSplot(rf_output, target, k=2, xlab="", ylab="", pch=target_labels, palette=c("Purple", "green"), main="MDS plot")
dev.off()

#Create ROC curve plot and calculate AUC
#Can use Relapse/non-relapse vote fractions as predictive variable
#The ROC curve will be generated by stepping up through different thresholds
predictions=as.vector(rf_output$votes[,2])
pred=prediction(predictions,target)
#First calculate the AUC value
perf_AUC=performance(pred,"auc")
AUC=perf_AUC@y.values[[1]]


#plot the actual ROC curve

perf_ROC=performance(pred,"tpr","fpr")
plot(perf_ROC, main="ROC for OOB Sample of PTSD", colorize= TRUE)

text(0.5,0.5,paste("AUC = ",format(AUC, digits=5, scientific=FALSE)))
dev.off()


#Save case predictions to file:
case_predictions=cbind(clindata,target,rf_output$predicted,rf_output$votes)
case_pred_outfile="trainset_CasePredictions.txt"
write.table(case_predictions,file=case_pred_outfile, sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
