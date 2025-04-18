setwd("figures/figure4")
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
#FPKM to TPM
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

setName <- "Cell_cycle"
setwd(paste0("signature/MsigDB/",setName))
filename <- list.files(pattern = "*.txt")
signature_list <- lapply(filename, function(x){gene = read.table(x, header = F, as.is = T); gene = gene$V1; return(gene[3:length(gene)])})
signature_name <- strsplit(filename, split = ".txt") %>% unlist()
names(signature_list) <- signature_name
print(length(signature_list))
#
setwd("figures/figure4")
print(paste("cluster num:",N))
consensus_type<-read.table(paste(condition,"/",type_output, "_subtype of Cluster_",N,".txt",sep = ""),as.is = T,header= T)
consensus_type$name <- rownames(consensus_type)
col_order <- read.table(paste(condition,"/",type_output, "_colOrder of Cluster_",N,".txt",sep = ""),header = F,as.is = T)
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
#4.ssgsea score
print("caculate the score of ssGSEA:")
res_ssGSEA <- ssgseaParam(as.matrix(exp_cluster_all),
		signature_list,
		minSize     = 2,
		normalize = TRUE)
res_ssGSEA <- gsva(res_ssGSEA)
write.table(res_ssGSEA, file = paste0(condition, "_res_ssGSEA.txt"), col.names = T, row.names = T, sep = "\t")
#
pValue <- list()
ssGSEA_list <- lapply(1:N, function(x){return(res_ssGSEA[, rownames(consensus_type)[consensus_type$clusterCS == paste0("C", x)]])})
plot_list <- list()
print(class(ssGSEA_list))
print(ssGSEA_list[[1]])
#
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
for(i in 1:1){
	score_list <- lapply(1:N, function(x){return(ssGSEA_list[[x]] %>% as.numeric())})
	pdat <- data.frame(value = unlist(score_list), subtype = c(rep(paste0("C", 1:N), times = lapply(1:N, function(x){return(length(score_list[[x]]))}) %>% unlist())))
	if(N == 2){
		pdat$subtype<- gsub(CIMP_neg, 'CIMP_neg', pdat$subtype)
		pdat$subtype<- gsub(CIMP_pos, 'CIMP_pos', pdat$subtype)
		pdat$subtype<- factor(pdat$subtype, levels = c("CIMP_neg", "CIMP_pos"))
	}
	
	name_title <- substr(rownames(res_ssGSEA)[i], 10, nchar(rownames(res_ssGSEA)[i]))
	if(N == 2){
		neg_value <- pdat[pdat$subtype=="CIMP_neg", "value"] %>% as.numeric()
		pos_value <- pdat[pdat$subtype=="CIMP_pos", "value"] %>% as.numeric()
		pdf(paste0(name_title, ".pdf"), width = 4, height = 6)
		p1 <- single_box("", neg_value, pos_value)
		print(p1)
		dev.off()
		pdf(paste0("Figure4C.pdf"), width = 4, height = 6)
		p2 <- violin_plot("", neg_value, pos_value)
		print(p2)
		dev.off()
	}
}
