rm(list = ls())
library(ComplexHeatmap)
library(pheatmap)
type_output <- "sdTop(1_per)"
N <- 2
condition <- "Temozolomide_LGG"
setwd("figures/figure1")
#Heatmap
anno_col <- list(Cluster=c(C1 = "#A6CEE3", C2 = "#1F78B4"))
beta_select <- read.table("beta_select.txt", header = T, sep = "\t", check.names = F)
consensus_type <- read.table(paste0(type_output, "_subtype of Cluster_", N , ".txt"))
col_order <- read.table(paste0(type_output, "_colOrder of Cluster_", N, ".txt"))
beta_select <- beta_select[, col_order$V1]
column_annotation <- HeatmapAnnotation(
				Cluster = consensus_type[colnames(beta_select), "clusterCS"],
				# cancerType = sample_info[colnames(beta_select), "cancer"],
				show_legend = T,
				col = anno_col
				)
result_hc <- Heatmap(beta_select, cluster_columns = F, cluster_rows = T, 
						 show_row_names = F, show_column_names = F,
						 # column_split = consensus_type[colnames(beta_select), "clusterCS"],
						 top_annotation = column_annotation,
						 heatmap_legend_param = list(title = "beta", at = c(0, 0.5, 1)),
						 column_title = "", col = c("#1D4A9E","white","#E7171B"))
pdf(file = paste0("Figure1B.pdf"), width = 5.5, height = 4)
print(result_hc)
dev.off()
