setwd("figures/figure4")
type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
CIMP_neg <- "C2"
CIMP_pos <- "C1"
library(ggplot2)
library(ggpubr)
library(GSVA) 
library(magrittr)
library(ComplexHeatmap)
library(patchwork)
#
mid_dir <- "cancer_specific_condition_data"
print("here")
setwd(workPath)
print("there")
#
print("read the data: ")
exp_mat <- read.table(paste0("/expression-FPKM/", condition, "/exp_mat.txt"), sep = "\t", as.is = T, header = T, check.names = F)
fpkmToTpm <- function(fpkm)
	{
		exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}
exp_mat <- 2^exp_mat -1
exp_mat <- apply(exp_mat, 2, fpkmToTpm)
exp_mat <- log(exp_mat+1, 2)
#sample information
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
setName <- "HALLMARK"
setwd(paste0("signature/MsigDB/",setName))
filename <- list.files(pattern = "*.txt")
signature_list <- lapply(filename, function(x){gene = read.table(x, header = F, as.is = T); gene = gene$V1; return(gene[3:length(gene)])})
signature_name <- strsplit(filename, split = ".txt") %>% unlist()
names(signature_list) <- signature_name
##hallmark classification
cellular_component <- paste0("HALLMARK_",c("APICAL_JUNCTION", "APICAL_SURFACE", "PEROXISOME"))
development <- paste0("HALLMARK_", c("ADIPOGENESIS", "ANGIOGENESIS", "EPITHELIAL_MESENCHYMAL_TRANSITION", "MYOGENESIS", "SPERMATOGENESIS", "PANCREAS_BETA_CELLS"))
DNA_damage <- paste0("HALLMARK_", c("DNA_REPAIR", "UV_RESPONSE_DN", "UV_RESPONSE_UP"))
immune <- paste0("HALLMARK_", c("ALLOGRAFT_REJECTION", "COAGULATION", "COMPLEMENT", "INTERFERON_ALPHA_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "IL6_JAK_STAT3_SIGNALING", "INFLAMMATORY_RESPONSE"))
metabolic <- paste0("HALLMARK_", c("BILE_ACID_METABOLISM", "CHOLESTEROL_HOMEOSTASIS", "FATTY_ACID_METABOLISM", "GLYCOLYSIS", "HEME_METABOLISM", "OXIDATIVE_PHOSPHORYLATION", "XENOBIOTIC_METABOLISM"))
pathway <- paste0("HALLMARK_", c("APOPTOSIS", "HYPOXIA", "PROTEIN_SECRETION", "UNFOLDED_PROTEIN_RESPONSE", "REACTIVE_OXYGEN_SPECIES_PATHWAY"))
proliferation <- paste0("HALLMARK_", c("E2F_TARGETS", "G2M_CHECKPOINT", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "P53_PATHWAY", "MITOTIC_SPINDLE"))
signaling <- paste0("HALLMARK_", c("ANDROGEN_RESPONSE", "ESTROGEN_RESPONSE_EARLY", "ESTROGEN_RESPONSE_LATE", "IL2_STAT5_SIGNALING", "KRAS_SIGNALING_UP", "KRAS_SIGNALING_DN", "MTORC1_SIGNALING", "NOTCH_SIGNALING", "PI3K_AKT_MTOR_SIGNALING", "HEDGEHOG_SIGNALING", "TGF_BETA_SIGNALING", "TNFA_SIGNALING_VIA_NFKB", "WNT_BETA_CATENIN_SIGNALING"))
##
print(paste("cluster num:",N))
consensus_type<-read.table(paste(condition, "/", type_output, "_subtype of Cluster_",N,".txt",sep = ""),as.is = T,header= T)
consensus_type$name <- rownames(consensus_type)
col_order <- read.table(paste(condition, "/", type_output, "_colOrder of Cluster_",N,".txt",sep = ""),header = F,as.is = T)
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
#ssGSEA scores
print("caculate the score of ssGSEA:")
res_ssGSEA <- ssgseaParam(as.matrix(exp_cluster_all),
			signature_list,
			minSize     = 2,
			normalize = TRUE)
res_ssGSEA <- gsva(res_ssGSEA)
pValue <- list()
direct <- list()
ssGSEA_list <- lapply(1:N, function(x){return(res_ssGSEA[, rownames(consensus_type)[consensus_type$clusterCS == paste0("C", x)]])})
plot_list <-list()
for(i in 1:dim(res_ssGSEA)[1]){
	score_list <- lapply(1:N, function(x){return(ssGSEA_list[[x]][i, ] %>% as.numeric())})
	pdat <- data.frame(value = unlist(score_list), subtype = c(rep(paste0("C", 1:N), times = lapply(1:N, function(x){return(length(score_list[[x]]))}) %>% unlist())))
	if(N == 2){
		pdat$subtype<- gsub(CIMP_neg, 'CIMP_neg', pdat$subtype)
		pdat$subtype<- gsub(CIMP_pos, 'CIMP_pos', pdat$subtype)
		pdat$subtype<- factor(pdat$subtype, levels = c("CIMP_neg", "CIMP_pos"))
	}
	if(N == 2){
		colour_man <- c("#00BFC4","#F8766D")
		pl <- ggboxplot(pdat, 
						x = "subtype", 
						y = "value", 
						fill = "subtype", 
						palette = colour_man, 
						width = 0.3, 
						size = 0.12,
						title = substr(rownames(res_ssGSEA)[i], 10, nchar(rownames(res_ssGSEA)[i])), 
						legend="none",
						add = "jitter",
						add.params=list(size=0.2))+
						theme(text = element_text(size = 8), plot.title = element_text(size = 8))+
						stat_compare_means(method="wilcox.test", size = 2)
		plot_list[[i]] <- pl
	}
	if(N == 2){
		p <- wilcox.test(score_list[[1]], score_list[[2]])$p.value
		pValue[[i]] <- p
		if(mean(score_list[[(as.numeric(substr(CIMP_neg,2,2)))]]) > mean(score_list[[(as.numeric(substr(CIMP_pos,2,2)))]]))
		{
			direct[[i]] <- "neg"
		}
		else {
			direct[[i]] <- "pos"
		}
		
	}
}
pl_sum <- plot_list[[1]]
for(i in 2:length(plot_list)){
	pl_sum <- pl_sum + plot_list[[i]]
}
pl_sum <- pl_sum + plot_layout(ncol = 10)
pdf(file = paste0(condition, "_ind_C", N, ".pdf"), width = 14, height = 9.5)
print(pl_sum)
dev.off()
pValue_info <- data.frame(signature = rownames(res_ssGSEA), p_value = unlist(pValue),
direct = unlist(direct))
pValue_info$FDR <- p.adjust(pValue_info$p_value, method = "fdr")
consensus_type$clusterCS <- gsub(CIMP_pos, 'CIMP_pos', consensus_type$clusterCS)
consensus_type$clusterCS <- gsub(CIMP_neg, 'CIMP_neg', consensus_type$clusterCS)
consensus_type$clusterCS <- factor(consensus_type$clusterCS, levels = c("CIMP_pos", "CIMP_neg"))
if(N==2){anno_col <- list(clusterCS=c(CIMP_neg = "#3C83C5", CIMP_pos = "#EF7B52"))}
column_annotation <- HeatmapAnnotation(
		clusterCS = consensus_type[colnames(exp_cluster_all), "clusterCS"],
		show_legend = T,
		col = anno_col
		)
pdat <- as.matrix(res_ssGSEA)
pdat <- apply(res_ssGSEA, 1, scale) %>% t()
rownames(pdat) <- rownames(res_ssGSEA)
colnames(pdat) <- colnames(res_ssGSEA)
pdat <- pdat[c(cellular_component, development, DNA_damage, immune, metabolic, pathway, proliferation, signaling), ]
print(head(pdat[, 1:2], 20))
rownames(pValue_info) <- pValue_info$signature
pValue_info <- pValue_info[rownames(pdat), ]
print(head(pValue_info))
write.csv(pValue_info, file = "HALLMARK_pvalue_info.csv", row.names = F)
print(pValue_info[pValue_info$FDR < 0.05, "signature"])
pdat <- pdat[pValue_info[pValue_info$FDR < 0.05, "signature"], ]
print(dim(pdat))
result_hc <- Heatmap(pdat, cluster_columns = F, cluster_rows = F, 
					 show_row_names = T, show_column_names = F,
					 row_names_gp = gpar(fontsize = 11), 
					 column_split = consensus_type[colnames(exp_cluster_all), "clusterCS"],
					 top_annotation = column_annotation,
					 column_title = "HALLMARK GENE SET",
					 heatmap_legend_param = list(title = "ssGSEA"), col = c("#1D4A9E", "white", "#E7171B"))
pdf(file = paste0("Figure4A.pdf"),width = 12,height = 10)
print(result_hc)
dev.off()
