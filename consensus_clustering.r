library(ggplot2)
library(ComplexHeatmap)
library(ConsensusClusterPlus)
library(pheatmap)
library(magrittr)
#
condition <- "Temozolomide_LGG"
type_condition <- "drug_cancer"
mid_dir <- "cancer_specific_condition_data"
print(paste("condition:",condition,sep = ""))
##############################################################################################
print("read the data:")
subPath <- ""
# cancer type information
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
#
print("read the data")
beta_res_select <- read.table(paste0(mid_dir, "/", condition, "/beta_res", subPath, ".txt"), sep = "\t", as.is = T, header= T, check.names = F)
##############################################################################################
sd_value <- unlist(apply(beta_res_select,1,sd))
sd_info <- data.frame(ID = rownames(beta_res_select), sd_value)
sd_info <- sd_info[order(sd_info$sd_value, decreasing = T),]
#
for(cutoff_control in c(1)){
	sd_cutoff <- sort(sd_value,decreasing = T)[length(sd_value)*(0.01*cutoff_control)]
	print(sd_cutoff)
	beta_select <- beta_res_select[rownames(beta_res_select)%in%(sd_info[sd_info$sd_value > sd_cutoff, "ID"]), ]
	write.table(beta_select, file = "beta_select.txt", row.names = T, col.names = T, sep = "\t", quote = F)
	#
	type_output <- paste0("sdTop(",cutoff_control, "_per)")
	if(dim(beta_select)[1] < 2) {print("beta_select num is too small!!!, next!"); next}
	if(dim(beta_select)[1] > 10000){print("beta_select num is too large!!!, next!"); next}
	################################################
	print("consensusClass:")
	title <- paste(type_output,"_Consensus",sep  = "")
	print(dim(beta_select))
	beta_select <- as.matrix(beta_select)
	print(class(beta_select))
	print(head(beta_select[, 1:4]))
	gradient_colors <- colorRampPalette(c("#29348B", "#9C1E28"))(10)
	#Figure 1A
	result_cs <- ConsensusClusterPlus(beta_select,maxK = 6,
								  clusterAlg = "km",distance = "euclidean",innerLinkage = "average",finalLinkage ="average",
								  tmyPal = gradient_colors, 
								  title = title,seed = 123456,plot = 'pdf')
	#Consensus matrix 
	write.table(result_cs[[2]]$consensusMatrix, file = "consensusMatrix_C2.txt", col.names = F, row.names = F, sep = "\t", quote =  F)
	cluster_cutoff_function <- function(N){
		clusterCS <- result_cs[[N]][["consensusClass"]]
		consensus_type <- data.frame(clusterCS=factor(paste0("C", clusterCS)))
		rownames(consensus_type) <- names(clusterCS)
		# obtain the samples order of clustering of consensus.
		cluster_order <- result_cs[[N]][["consensusTree"]][["order"]]
		beta_select <- beta_select[,cluster_order]
		write.table(colnames(beta_select), file = paste0(type_output, "_colOrder of Cluster_", N, ".txt"), col.names = F,  row.names = F)
		#2. hierarchical clustering of probes.(col order keep the same as consensus order)
		print("probe clustering:")
		if(N==2){anno_col <- list(clusterCS=c(C1 = "#A6CEE3", C2 = "#1F78B4"))}
		title <- paste0(type_output, "\t", "Probes: ", dim(beta_select)[1])
		#########
		select_probes <- rownames(beta_select)
		#########
		column_annotation <- HeatmapAnnotation(
			clusterCS = consensus_type[colnames(beta_select), "clusterCS"],
			cancerType = sample_info[colnames(beta_select), "cancer"],
			show_legend = T,
			col = anno_col
			)
		result_hc <- Heatmap(beta_select, cluster_columns = F, cluster_rows = T, 
					 show_row_names = F, show_column_names = F,
					 # column_split = consensus_type[colnames(beta_select), "clusterCS"],
					 top_annotation = column_annotation,
					 heatmap_legend_param = list(title = "beta", at = c(0, 0.5, 1)),
					 column_title = "")
		select_probes <- select_probes[result_hc@row_order]
		write.table(select_probes, file = paste(type_output,"_selectProbes.txt",sep = ""), row.names = F, col.names = F,quote = F)
		#
		pdf(file = paste0(type_output, "_clustering of Cluster_", N, ".pdf"), width = round(8*dim(beta_select)[2]/72, 1), height = 6)
		print(result_hc)
		dev.off()
		filename <- paste0(type_output, "_subtype of Cluster_", N , ".txt")
		write.table(consensus_type, file = filename, col.names = T, row.names = T)
		print("complete!!")
	}
	lapply(2,cluster_cutoff_function)#
}
