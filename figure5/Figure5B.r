setwd("figures/figure5")
type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
drug_id <- "AZD6482_2169"
CIMP_neg <- "C2"
CIMP_pos <- "C1"
GDSC <- "GDSC2"
library(magrittr)
library(ggpubr)
#1 data preparing
print(paste0("condition: ", condition))
drug <- unlist(strsplit(condition, split = "_"))[1]
subtype <- read.table(paste0(condition,"/",type_output,"_subtype of Cluster_",N,".txt"), as.is = T, header= T)
drugPredict <- read.csv(file = paste0("calcPhenotype_Output/DrugPredictions.csv"), check.names = F)
cell_name_list <- lapply(1:N, function(x){return(intersect(drugPredict[, 1], rownames(subtype)[subtype$clusterCS==paste0("C", x)]))})
#
cell_res_list <- lapply(1:N, function(x){return(drugPredict[drugPredict[, 1]%in%cell_name_list[[x]], drug_id] %>% as.numeric())})
pdat <- data.frame(response = unlist(cell_res_list), subtype = rep(paste0("C", 1:N), times = unlist(lapply(1:N, function(x){return(length(cell_res_list[[x]]))}))))
#
pdat$subtype <- gsub(CIMP_neg, 'CIMP_neg', pdat$subtype)
pdat$subtype <- gsub(CIMP_pos, 'CIMP_pos', pdat$subtype)
pdat$subtype <- factor(pdat$subtype, levels = c("CIMP_neg", "CIMP_pos"))
#
neg_value <- pdat[pdat$subtype=="CIMP_neg", "response"] %>% as.numeric()
pos_value <- pdat[pdat$subtype=="CIMP_pos", "response"] %>% as.numeric()
#function
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
single_box <- function(title_name, gene_normal, gene_tumor){
library(ggpubr)
pdat <- data.frame(value = c(gene_normal, gene_tumor), type = c(rep("CIMP_neg", length(gene_normal)), rep("CIMP_pos", length(gene_tumor))))
pdat$type <- factor(pdat$type, levels = c("CIMP_neg", "CIMP_pos"))
# colour_man <- c("#00BFC4", "#F8766D")
colour_man <- c("#3C83C5","#EF7B52")
p1 <- ggboxplot(pdat,
				x = "type",
				y = "value",
				fill = "type",
				palette = colour_man,
				width = 0.5,
				legend = "none",
				title = title_name,
				add = "jitter",
				#outlier.shape = NA,
				add.params = list(size=0.5)
				# bxp.errorbar = TRUE
				)+
				theme(text = element_text(size = 15),
				plot.title = element_text(hjust=0.5))+
				stat_compare_means(method = "wilcox.test", size = 6)
return(p1)
}
p1 <- single_box("", neg_value, pos_value)
pdf(file = paste0(condition,"_",drug_id,"_prediction_", GDSC, ".pdf"),width = 4, height = 6)
print(p1)
dev.off()
p2 <- violin_plot("", neg_value, pos_value)
pdf(file = paste0("Figure5B.pdf"),width = 4, height = 6)
print(p2)
dev.off()
