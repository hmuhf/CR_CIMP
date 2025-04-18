setwd("figures/figure4")
source("/data/home/houfei/nuaa/301/program/function.R")
type_condition <- "drug_cancer"
condition <- "Gemcitabine_PAAD"
type_output <- "sdTop(1_per)"
N <- 2
CIMP_neg <- "C1"
CIMP_pos <- "C2"
library(ggplot2)
library(ggpubr)
library(GSVA)
library(magrittr)
library(ComplexHeatmap)
library(patchwork)
#
mid_dir <- "cancer_specific_condition_data"
#
print("read the data: ")
exp_mat<-read.table(paste0("/expression-FPKM/", condition, "/exp_mat.txt"), sep = "\t", as.is = T, header = T,check.names = F)
#
#
print(head(exp_mat[, 1:4]))
exp_mat <- 2^exp_mat-1
print(head(exp_mat[, 1:4]))
fpkmToTpm <- function(fpkm)
{
	exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
exp_mat <- apply(exp_mat, 2, fpkmToTpm)
print(dim(exp_mat))
print(head(exp_mat[, 1:4]))
#
dir <- paste0(type_condition, "/overlap_name/methylation/", condition, "/")
file <- list.files(dir, pattern = "*_overlap.txt")
print(file)
data_list <- lapply(paste0(dir, file), read.table, sep = "\t", header = F, as.is = T, na.strings = "")
sample_info <- list()
for (i in 1:length(data_list)){
	tmp <- data_list[[i]]
	tmp[, "cancer"] <- unlist(strsplit(file[i], split = "_"))[1]
	sample_info[[i]] <- tmp
}
sample_info <- do.call(rbind, sample_info)
colnames(sample_info)[1] <- "sampleID"
rownames(sample_info) <- sample_info$sampleID
#msigdb signature
setName <- "DDR"
setwd(paste0("signature/MsigDB/",setName))
filename <- list.files(pattern = "*.txt")
signature_list <- lapply(filename, function(x){gene = read.table(x, header = F, as.is = T); gene = gene$V1; return(gene[3:length(gene)])})
signature_name <- strsplit(filename, split = ".txt") %>% unlist()
names(signature_list) <- signature_name
print(length(signature_list))
#
setwd("figures/figure4")
print(paste("cluster num:",N))
consensus_type<-read.table(paste(condition, "/",type_output, "_subtype of Cluster_",N,".txt",sep = ""),as.is = T,header= T)
consensus_type$name <- rownames(consensus_type)
col_order <- read.table(paste(condition, "/",type_output, "_colOrder of Cluster_",N,".txt",sep = ""),header = F,as.is = T)
col_order <- col_order$V1
overlap_sample <- intersect(rownames(consensus_type), colnames(exp_mat))
na_sample <- setdiff(rownames(consensus_type), overlap_sample)
print(na_sample)
exp_mat <- exp_mat[, overlap_sample]
print(dim(exp_mat))
consensus_type <- consensus_type[overlap_sample, ]
print(dim(consensus_type))
col_order <- col_order[col_order%in%overlap_sample]
print(length(col_order))
exp_cluster_all <- exp_mat[, col_order]
#4.
print("caculate the score of ssGSEA:")
res_ssGSEA <- ssgseaParam(as.matrix(exp_cluster_all),
	signature_list,
	minSize     = 2,
	normalize = TRUE)
res_ssGSEA <- gsva(res_ssGSEA)
print(head(res_ssGSEA[,1:4]))
write.table(res_ssGSEA, file = paste0(condition, "_ddr_ssgsea.txt"), col.names = T, row.names = T, sep = "\t")
#
pValue <- list()
ssGSEA_list <- lapply(1:N, function(x){return(res_ssGSEA[, rownames(consensus_type)[consensus_type$clusterCS == paste0("C", x)]])})
plot_list <- list()
for(i in 1:dim(res_ssGSEA)[1]){
	score_list <- lapply(1:N, function(x){return(ssGSEA_list[[x]][i, ] %>% as.numeric())})
	pdat <- data.frame(value = unlist(score_list), subtype = c(rep(paste0("C", 1:N), times = lapply(1:N, function(x){return(length(score_list[[x]]))}) %>% unlist())))
	
	pdat$subtype<- gsub(CIMP_neg, 'CIMP_neg', pdat$subtype)
	pdat$subtype<- gsub(CIMP_pos, 'CIMP_pos', pdat$subtype)
	pdat$subtype<- factor(pdat$subtype, levels = c("CIMP_neg", "CIMP_pos"))
	name_title <- rownames(res_ssGSEA)[i]
	colour_man <- c("#00BFC4","#F8766D")
	pl <- ggboxplot(pdat, 
						x = "subtype", 
						y = "value", 
						fill = "subtype", 
						palette = colour_man, 
						width = 0.3, 
						size = 0.2,
						title = name_title,
						legend="none",
						add = "jitter",
						add.params=list(size=0.5))+
						theme(text = element_text(size = 10), plot.title = element_text(size = 10))+
						stat_compare_means(method="wilcox.test", size = 3)
	plot_list[[i]] <- pl
	
	
	p <- wilcox.test(score_list[[1]], score_list[[2]])$p.value
	pValue[[i]] <- p
}
pl_sum <- plot_list[[1]]
for(i in 2:length(plot_list)){
	pl_sum <- pl_sum + plot_list[[i]]
}
pl_sum <- pl_sum + plot_layout(ncol = 5)
pdf(file = paste0(condition, "_ddr_c", N, ".pdf"), width = 12, height = 14)
print(pl_sum)
dev.off()
pValue_info <- data.frame(signature = rownames(res_ssGSEA), p_value = unlist(pValue))
pValue_info <- pValue_info[order(pValue_info$p_value), ]
write.table(pValue_info, file = paste0(condition,"_ddr_test_c",N,".txt"), col.names = T, row.names = F, sep = "\t", quote = F)
# 5 visulization
if(N == 2){
	consensus_type$clusterCS <- gsub(CIMP_neg, 'CIMP_neg', consensus_type$clusterCS)
	consensus_type$clusterCS <- gsub(CIMP_pos, 'CIMP_pos', consensus_type$clusterCS)
	consensus_type$clusterCS <- factor(consensus_type$clusterCS, levels = c("CIMP_neg", "CIMP_pos"))
}
if(N==2){anno_col <- list(clusterCS=c(CIMP_neg = "#3C83C5",CIMP_pos = "#EF7B52"))}
column_annotation <- HeatmapAnnotation(
		clusterCS = consensus_type[colnames(exp_cluster_all), "clusterCS"],
		show_legend = F,
		col = anno_col
		)
pdat <- as.matrix(res_ssGSEA)
pdat <- apply(res_ssGSEA, 1, scale) %>% t()
rownames(pdat) <- rownames(res_ssGSEA)
colnames(pdat) <- colnames(res_ssGSEA)
result_hc <- Heatmap(pdat, cluster_columns = F, cluster_rows = F, 
					 show_row_names = T, show_column_names = F,  
					 row_names_gp = gpar(fontsize = 6), 
					 column_split = consensus_type[colnames(exp_cluster_all), "clusterCS"],
					 # column_names_rot = 45,
					 column_names_gp = gpar(fontsize = 8),
					 top_annotation = column_annotation,
					 column_title = type_output,
					 heatmap_legend_param = list(title = "ssGSEA"), col = c("#1D4A9E", "white", "#E7171B"))
pdf(file = "Figure4D.pdf", width = 5, height = 3)
print(result_hc)
dev.off()
