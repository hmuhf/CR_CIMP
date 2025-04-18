####boruta+random forest
type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
CIMP_neg <- "C2"
CIMP_pos <- "C1"
mid_dir <- "cancer_specific_condition_data"
##
library(caret)
library(magrittr)
library(randomForest)
library(Boruta)
library(dplyr)
beta_res <- read.table(paste0( mid_dir, "/", condition, "/beta_res.txt"), sep = "\t", as.is = T, header= T, check.names = F)
select_probes <- read.table(paste0(condition,"/",type_output, "_selectProbes.txt"), header = F, as.is = T)
subtype<-read.table(paste0(condition, "/", type_output,"_subtype of Cluster_",N,".txt"),as.is = T,header= T)
subtype$sample <- rownames(subtype)
subtype$clusterCS[subtype$clusterCS==CIMP_neg] <- "CIMP_neg"
subtype$clusterCS[subtype$clusterCS==CIMP_pos] <- "CIMP_pos"
beta_res_select <- beta_res[select_probes$V1, ]
dat_matrix <- as.data.frame(t(beta_res_select))
dat_matrix$y <- factor(subtype[rownames(dat_matrix), "clusterCS"], levels = c("CIMP_neg", "CIMP_pos"))
#boruta 
set.seed(777)
boruta <- Boruta(x=dat_matrix[, -which(colnames(dat_matrix)=="y")], y=dat_matrix$y, pValue=0.05, mcAdj=T, maxRuns=300)
boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(boruta, withTentative = F), Type="Boruta_without_tentative")
select_feature <- boruta.finalVarsWithTentative$Item
#save
write.table(select_feature, file = "all_sample_boruta_select_feature.txt", quote = F, col.names = F, row.names = F)
#
dat_matrix <- dat_matrix[, c(select_feature, "y")]
CIMP_neg_matrix <- dat_matrix[dat_matrix$y=="CIMP_neg", ]
CIMP_pos_matrix <- dat_matrix[dat_matrix$y=="CIMP_pos", ]
#heatmap of select CpGs
library(ComplexHeatmap)
row_annotation <- rowAnnotation(
	type = c(rep("neg", dim(CIMP_neg_matrix)[1]), rep("pos", dim(CIMP_pos_matrix)[1])),
	show_legend = T
	)

pdat <- rbind(CIMP_neg_matrix[, -which(colnames(CIMP_neg_matrix)=="y")], 
	      CIMP_pos_matrix[, -which(colnames(CIMP_pos_matrix)=="y")])

result_hc <- Heatmap(pdat, cluster_columns = F, cluster_rows = F, 
			 show_row_names = F, show_column_names = F,
			 left_annotation = row_annotation,
			 show_heatmap_legend = T,
			 column_title = ""
			 # column_split = CIMP_info$type,
			 )
pdf(file = paste0("Resistant_Sample_heatmap.pdf"),width = 10,height = 6)
print(result_hc)
dev.off()

#cross-validation
CVgroup <- function(k,datasize,seed){
	set.seed(seed)
	n <- rep(1:k,ceiling(datasize/k))[1:datasize]
	temp <- sample(n,datasize)
	x <- 1:k
	dataseq <- 1:datasize
	cvlist <- lapply(x,function(x) dataseq[temp==x])
	return(cvlist)
}
k=3
result <- c()
time <- 1
seed=123
for(i in 1:10){
cvlist_1 <- CVgroup(k=k,datasize=nrow(CIMP_pos_matrix),seed=seed)
cvlist_2 <- CVgroup(k=k,datasize=nrow(CIMP_neg_matrix),seed=seed)
#print(cvlist)
for(j in 1:k){
	for(m in 1:k){
		resistant_test <- CIMP_pos_matrix[cvlist_1[[j]],]
		resistant_train <- CIMP_pos_matrix[-cvlist_1[[j]],]
		sensitive_test <- CIMP_neg_matrix[cvlist_2[[m]],]
		sensitive_train <- CIMP_neg_matrix[-cvlist_2[[m]],]
		test <- rbind(resistant_test,sensitive_test)
		test <- as.data.frame(test)
		train <- rbind(resistant_train,sensitive_train)
		train <- as.data.frame(train)
		train$y <- factor(train$y)
		print(dim(train))
		print(dim(test))
		#
		otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
		pred <- predict(otu_train.forest,test,type="prob")
		each_result <- data.frame(prediction=pred[,2],real_data=test$y,time=time)
		result <- rbind(result,each_result)
		time <- time+1
	}
}
seed <- seed+1
print(seed)
}
result <- as.data.frame(result)

write.csv(result, file = "boruta_rf_train_result.csv", row.names = F)
#ROC
result <-read.csv("boruta_rf_train_result.csv")
library(pROC)
roc_curve <- roc(result[, 2], result[, 1])
auc_each <- round(auc(roc_curve),2)
print(auc_each)
pdf(file = "boruta_cv_roc.pdf")
plot(roc_curve)
dev.off()
#final model using selected CpGs with all samples
train <- dat_matrix
print(dim(train))
#
otu_train.forest <- randomForest(y~., data = train, importance = TRUE)
save(otu_train.forest, file = "train_model.rda")
