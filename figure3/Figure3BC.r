setwd("figures/figure3")
type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
CIMP_neg <- "C2"
CIMP_pos <- "C1"
#library(maftools)
library(magrittr)
library(ggpubr)
subtype <- read.table(paste(type_output,"_subtype of Cluster_",N,".txt",sep = ""),as.is = T, header= T)
subtype$sample <- rownames(subtype)
subtype$clusterCS[subtype$clusterCS==CIMP_neg] <- "CIMP_neg"
subtype$clusterCS[subtype$clusterCS==CIMP_pos] <- "CIMP_pos"
print(head(subtype))
#immune_score from literature
immune_score <- read.table("TCGA_signature_score.txt", header = T, sep = "\t", na.strings = "NA", fill = T, quote = "")
print(colnames(immune_score))
immune_score$TCGA.Participant.Barcode <- paste0(immune_score$TCGA.Participant.Barcode, "-01A")
print(head(immune_score))
#
feature <- c("SNV.Neoantigens", "Aneuploidy.Score")
for(i in feature){
	value_CIMP_neg <- immune_score[immune_score$TCGA.Participant.Barcode%in%subtype[subtype$clusterCS=="CIMP_neg", "sample"], i] %>% na.omit()
	value_CIMP_pos <- immune_score[immune_score$TCGA.Participant.Barcode%in%subtype[subtype$clusterCS=="CIMP_pos", "sample"], i] %>% na.omit()
	print(paste0("CIMP_neg num-", i,": ", length(value_CIMP_neg)))
	print(paste0("CIMP_pos num-", i,": ", length(value_CIMP_pos)))
	value <- c(value_CIMP_neg, value_CIMP_pos)
	group <- c(rep("CIMP_neg", length(value_CIMP_neg)), rep("CIMP_pos", length(value_CIMP_pos)))
	pdat <- data.frame(value, group)
	ylab <- "score"
	if(i=="SNV.Neoantigens"){
		value_CIMP_neg = log(value_CIMP_neg, 2)
	value_CIMP_pos = log(value_CIMP_pos, 2)
	}
	#plot function
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
	pdf(file = paste0(i,"_score",".pdf"),width = 4,height = 6)
	p1 <- single_box("", value_CIMP_neg, value_CIMP_pos)
	print(p1)
	dev.off()
	pdf(file = paste0("Figure2BC.pdf"),width = 4,height = 6)
	p2 <- violin_plot("", value_CIMP_neg, value_CIMP_pos)
			print(p2)
	dev.off()
}
