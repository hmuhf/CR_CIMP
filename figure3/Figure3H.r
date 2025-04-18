setwd("figures/figure3")
type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
gene <- "MEOX2"
CIMP_neg <- "C2"
CIMP_pos <- "C1"
library(magrittr)
print(paste0("condition: ", condition))
drug <- unlist(strsplit(condition, split = "_"))[1]
#expression comparison
#methylation comparison
#correlation between expression and methylation
exp_mat <- read.table(paste0("expression-FPKM/", condition, "/exp_mat.txt"), sep = "\t", as.is = T, header = T, na.strings = "", check.names = F)
#
fpkmToTpm <- function(fpkm)
	{
		exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}
exp_mat <- 2^exp_mat - 1#
exp_mat <- apply(exp_mat, 2, fpkmToTpm)
exp_mat <- log(exp_mat+1, 2)
#
meth_mat <- read.table(paste0(condition, "/promoter_island/cor/",gene,"_cpgs_beta.txt"), sep = "\t", as.is = T, header= T, check.names = F)
#exp
exp_value <- exp_mat[gene, ] %>% as.numeric()
names(exp_value) <- colnames(exp_mat)
#methylation
meth_value <- meth_mat %>% apply(., 2, mean) %>% unlist() %>% as.numeric()
names(meth_value) <- colnames(meth_mat)
#subtype 
subtype <- read.table(paste(condition,"/",type_output,"_subtype of Cluster_",N,".txt",sep = ""),as.is = T,header= T)
#
meth_mat_neg <- meth_mat[,colnames(meth_mat)%in%rownames(subtype)[subtype$clusterCS==CIMP_neg]]
meth_mat_pos <- meth_mat[,colnames(meth_mat)%in%rownames(subtype)[subtype$clusterCS==CIMP_pos]]
meth_mat_neg_cpgs <- meth_mat_neg %>% apply(., 1, mean) %>% unlist() %>% as.numeric()
meth_mat_pos_cpgs <- meth_mat_pos %>% apply(., 1, mean) %>% unlist() %>% as.numeric()
meth_mat_cpgs <- cbind(meth_mat_neg_cpgs, meth_mat_pos_cpgs)
colnames(meth_mat_cpgs) <- c("neg_value", "pos_value")
rownames(meth_mat_cpgs) <- rownames(meth_mat)
write.table(meth_mat_cpgs, file = paste0(gene, "_meth_mat_cpgs.txt"), col.names = T, row.names = T,sep = "\t", quote = F)
#
meth_neg_value <- meth_value[rownames(subtype)[subtype$clusterCS==CIMP_neg]]
meth_pos_value <- meth_value[rownames(subtype)[subtype$clusterCS==CIMP_pos]]
#
exp_neg_value <- exp_value[rownames(subtype)[subtype$clusterCS==CIMP_neg]]
exp_pos_value <- exp_value[rownames(subtype)[subtype$clusterCS==CIMP_pos]]
#visulization
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
					# xlab("PAAD\n(num(T)=178; num(N)=167)")+
					# ylab("cibersort score")+
					theme(text = element_text(size = 15),
					plot.title = element_text(hjust=0.5))+
					stat_compare_means(method = "wilcox.test", size = 6)
	return(p1)
}
vis_cor <- function(title_name, value_x, value_y){
	library(ggpubr)
	dat <- data.frame(value_x, value_y)
	color_man <- "#0572BE"
	fill_man <- "#ABC9E6"
	p1 <- ggscatter(dat, x = "value_x", y = "value_y",
			add = "reg.line",
			conf.int = TRUE,
			color = color_man,
			add.params = list(color = color_man,
                            fill = fill_man),
			ggtheme = theme_minimal(),
			title = title_name)+
			border()
	p2 <- p1+stat_cor(method = "pearson", color = color_man, size = 6, fontface="bold")+
			xlab("Methylation")+
			ylab("Expression")+
			theme(plot.title = element_text(hjust=0.5),
					text = element_text(size = 15),
					,axis.text=element_text(size=15, color = "black"))
	return(p2)
}
pdf(paste0(gene, "_meth_comp.pdf"), width = 4, height = 6)
picture <- single_box("", meth_neg_value, meth_pos_value)
print(picture)
dev.off()
pdf(paste0(gene, "_meth_comp_violin.pdf"), width = 4, height = 6)
p2 <- violin_plot("", meth_neg_value, meth_pos_value)
print(p2)
dev.off()
#
pdf(paste0(gene, "_exp_comp.pdf"), width = 4, height = 6)
picture <- single_box(gene, exp_neg_value, exp_pos_value)
print(picture)
dev.off()
pdf(paste0(gene, "_exp_comp_violin.pdf"), width = 4, height = 6)
p2 <- violin_plot(gene, exp_neg_value, exp_pos_value)
print(p2)
dev.off()
#cor
exp_value <- exp_value[names(meth_value)]
print(length(exp_value))
print(length(meth_value))
print(exp_value)
print("split")
print(meth_value)
pdf(paste0("Figure3H.pdf"), width = 4, height = 4)
picture <- vis_cor(gene, meth_value, exp_value)
print(picture)
dev.off()
