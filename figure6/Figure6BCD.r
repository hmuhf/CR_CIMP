library(ggplot2)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(pheatmap)
library(magrittr)
condition <- "Temozolomide_LGG"
setwd("figures/figure6")
#
dat <- read.table("CIMP_validation/GSE104293_valid.txt", sep = "\t", quote = "", header = T, as.is = T)
rownames(dat) <- dat$ID_REF
dat <- dat[, -1]
print(dim(dat))
print(dat[1:3,1:3])
#low probes filtering
txt_files <- list.files("CIMP_validation/detective_p_value", pattern = "*.txt")
filename <- paste0("CIMP_validation/detective_p_value/", txt_files)
data_list <- lapply(filename,read.table,sep='\t',header=TRUE)
det_p <- c()
for (i in 1:length(data_list)){
	det_p <- cbind(det_p, as.numeric(data_list[[i]][, 3]))
}
det_p <- as.data.frame(det_p)
colnames(det_p) <- gsub('-[0-9][0-9][0-9][0-9].txt','',txt_files)
rownames(det_p) <- data_list[[1]][, 1]
Pfailed <- det_p > 0.01
failedProbes_pvalue <-rownames(Pfailed)[rowMeans(Pfailed)>0.2]
anno_450k <- read.csv("humanmethylation450_15017482_v1-2.csv")
failedProbes_XY <- anno_450k$Name[anno_450k$CHR%in%c("X","Y")]
dat_filter <- dat[!rownames(dat)%in%c(failedProbes_XY, failedProbes_pvalue), ]
#
sd_value <- unlist(apply(dat_filter,1,sd))
sd_info <- data.frame(ID = rownames(dat_filter), sd_value)
sd_info <- sd_info[order(sd_info$sd_value, decreasing = T),]
write.table(sd_info, file = "GSE104293_sd_info.txt",row.names = F, col.names = T, sep = "\t", quote = F)
#
cutoff_control <- 1
sd_cutoff <- sort(sd_value,decreasing = T)[length(sd_value)*(0.01*cutoff_control)]
print(sd_cutoff)
sd_top1_cpgs <- sd_info[sd_info$sd_value > sd_cutoff, "ID"]
boruta_select_cpgs <- read.table("all_sample_boruta_select_feature.txt", header = F, as.is = T)
beta_select <- dat_filter[rownames(dat_filter)%in%(sd_info[sd_info$sd_value > sd_cutoff, "ID"]), ]
#
type_output <- paste0("GSE104293_sdTop(",cutoff_control, "_per)")
print("consensusClass:")
title <- paste(type_output,"_Consensus",sep  = "")
print(dim(beta_select))
beta_select <- as.matrix(beta_select)
print(class(beta_select))
print(head(beta_select[, 1:4]))
gradient_colors <- colorRampPalette(c("#29348B", "#9C1E28"))(10)
result_cs <- ConsensusClusterPlus(beta_select,maxK = 6,clusterAlg = "km",distance = "euclidean",innerLinkage = "average",finalLinkage ="average", title = title,seed = 123456,plot = 'pdf',tmyPal = gradient_colors)
#
N <- 2
clusterCS <- result_cs[[N]][["consensusClass"]]
consensus_type <- data.frame(clusterCS=factor(paste0("C", clusterCS)))
rownames(consensus_type) <- names(clusterCS)
# obtain the samples order of clustering of consensus.
cluster_order <- result_cs[[N]][["consensusTree"]][["order"]]
beta_select <- beta_select[,cluster_order]
write.table(colnames(beta_select), file = paste0(type_output, "_colOrder of Cluster_", N, ".txt"), col.names = F,  row.names = F)
print("probe clustering:")
anno_col <- list(clusterCS=c(C1 = "#A6CEE3", C2 = "#1F78B4"))
#anno_col <- list(clusterCS=c(C1 = "#3C83C5", C2 = "#EF7B52"))
title <- paste0(type_output, "\t", "Probes: ", dim(beta_select)[1])
select_probes <- rownames(beta_select)
column_annotation <- HeatmapAnnotation(
                                clusterCS = consensus_type[colnames(beta_select), "clusterCS"],
                                show_legend = T,
                                col = anno_col
                                )
result_hc <- Heatmap(beta_select, cluster_columns = F, cluster_rows = T,
                                                 show_row_names = F, show_column_names = F,
                                                 # column_split = consensus_type[colnames(beta_select), "clusterCS"],
                                                 top_annotation = column_annotation,
                                                 heatmap_legend_param = list(title = "beta", at = c(0, 0.5, 1)),
                                                 col = c("#1D4A9E", "white", "#E7171B"),
						 column_title = "")
select_probes <- select_probes[result_hc@row_order]
write.table(select_probes, file = paste(type_output,"_selectProbes.txt",sep = ""), row.names = F, col.names = F,quote = F)
pdf(file = paste0(type_output, "_clustering of Cluster_", N, ".pdf"), width = 5, height = 6)
print(result_hc)
dev.off()
filename <- paste0(type_output, "_subtype of Cluster_", N , ".txt")
write.table(consensus_type, file = filename, col.names = T, row.names = T)
#boruta
select_probes <- read.table("all_sample_boruta_select_feature.txt", header = F, as.is = T)
select_probes <- select_probes$V1
print(length(select_probes))
#
dat_select <- dat_filter[select_probes, ]
#
load("train_model.rda")
#
#prediction based on train_model
library(randomForest)
independent_dat <- t(cbind(dat_select))
pred <- predict(otu_train.forest, independent_dat, type="prob")
#
library(pROC)
type_actual <- read.table(paste0(type_output, "_subtype of Cluster_", N , ".txt"), header = T, as.is = T)
type_actual$sample <- rownames(type_actual)
type_actual$clusterCS[type_actual$clusterCS=="C1"] <- "CIMP_neg"
type_actual$clusterCS[type_actual$clusterCS=="C2"] <- "CIMP_pos"
type_actual <- type_actual[rownames(pred), ]
roc_curve <- roc(type_actual$clusterCS, pred[,2])
#AUC
auc <- round(auc(roc_curve),2)
print(auc)
pdf(file = "GSE104293_roc.pdf")
plot(roc_curve)
dev.off()
#
pred <- predict(otu_train.forest, independent_dat, type="response")
pred_info <- data.frame(sample = names(pred), type = as.character(pred))
write.csv(pred_info, file = "GSE104293_predict_CIMP.csv")
acc <- 14/16
print(paste0("ACC: ",round(acc, 2)))#0.88
spc <- 4/6#0.67
sen <- 10/10#1
#
pred_info <- pred_info[order(pred_info$type), ]
dat_select <- dat_select[, pred_info$sample]
anno_col <- list(predict_type=c(CIMP_neg = "#A6CEE3", CIMP_pos = "#1F78B4"),
		 actual_type = c(CIMP_neg ="#3C83C5", CIMP_pos = "#EF7B52"))
column_annotation <- HeatmapAnnotation(
                                predict_type = pred_info$type,
                                actual_type = type_actual[pred_info$sample, "clusterCS"],
				show_legend = T,
                                col = anno_col
                                )
result_hc <- Heatmap(dat_select, cluster_columns = F, cluster_rows = T,
                                                 show_row_names = F, show_column_names = F,
                                                 # column_split = consensus_type[colnames(beta_select), "clusterCS"],
                                                 top_annotation = column_annotation,
                                                 heatmap_legend_param = list(title = "beta", at = c(0, 0.5, 1)),
						 col = c("#1D4A9E", "white", "#E7171B"),
                                                 column_title = "")
pdf("Figure6D.pdf", width = 5, height = 6)
print(result_hc)
dev.off()
#
library(openxlsx)
clinic_info <-read.xlsx("CIMP_validation/GSE104293.clinical.processed.xlsx")
type_actual <- read.table(paste0(type_output, "_subtype of Cluster_", N , ".txt"), header = T, as.is = T)
type_actual$sample <- rownames(type_actual)
type_actual$clusterCS[type_actual$clusterCS=="C1"] <- "CIMP_neg"
type_actual$clusterCS[type_actual$clusterCS=="C2"] <- "CIMP_pos"
#
age_neg <- clinic_info[clinic_info$sample%in%type_actual$sample[type_actual$clusterCS=="CIMP_neg"], "age"] %>% as.numeric()
age_pos <- clinic_info[clinic_info$sample%in%type_actual$sample[type_actual$clusterCS=="CIMP_pos"], "age"] %>% as.numeric()
violin_plot <- function(title_name, gene_normal, gene_tumor){
	library(ggpubr)
	pdat <- data.frame(value = c(gene_normal, gene_tumor), type = factor(c(rep("CIMP_neg", length(gene_normal)), rep("CIMP_pos", length(gene_tumor))), levels = c("CIMP_neg","CIMP_pos")))
	colour_man <- c("#3C83C4","#EF7B52")
	p1 <- ggviolin(pdat, x = "type", y = "value", width=0.5, fill="type",size=0.3, palette = colour_man, add = "boxplot", add.params = list(size=0.5, fill = "white"), legend = "none") +
	theme(text = element_text(size = 15),
          plot.title = element_text(hjust=0.5))+
          stat_compare_means(method = "wilcox.test", size = 6)
	return(p1)
}
pdf("Figure6C.pdf", width = 5, height = 6)
p1 <- violin_plot("", age_neg, age_pos)
print(p1)
dev.off()
#
surv_info <- merge(type_actual, clinic_info, by = "sample", all.x = T, all.y =F)
surv_info$pfs_month <- as.numeric(surv_info$pfs)
surv_info$pfs_status <- as.numeric(surv_info$pfs_status)
surv_info$clusterCS <- factor(surv_info$clusterCS, levels = c("CIMP_neg","CIMP_pos"))
library(survival)
library(survminer)
pdf("Figure6B.pdf", width = 8, height = 6)
fit <- survfit(Surv(pfs_month, pfs_status) ~ clusterCS, data = surv_info)
colour_man <- c("#3C83C5","#EF7B52")
p1 <- ggsurvplot(
                   data = surv_info,
                   fit,
                   pval = TRUE,
                   conf.int = TRUE,
                   risk.table = TRUE,
                   risk.table.col = "strata",
                   linetype = "strata",
                   surv.median.line = "hv",
                   ggtheme = theme_bw(),
                   palette = colour_man,
                   font.main = 18,     #
                   font.x = 15,        #
                   font.y = 15,        #
                   font.tickslab = 15
                   )
print(p1)
dev.off()
