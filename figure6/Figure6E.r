library(magrittr)
condition <- "Temozolomide_LGG"
setwd("figures/figure6")
#########1
probes <- read.table("all_sample_boruta_select_feature.txt", header = F, as.is = T)
probes <- probes$V1
#########
exp_mat <- read.table(paste0("expression-FPKM/", condition, "/exp_mat.txt"), sep = "\t", as.is = T, header = T, na.strings = "", check.names = F)
#FPKM to tpm
fpkmToTpm <- function(fpkm)
	{
		exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}
exp_mat <- 2^exp_mat - 1#
exp_mat <- apply(exp_mat, 2, fpkmToTpm)
exp_mat <- log(exp_mat+1, 2)
#
beta_res <- read.table(paste0(condition, "/beta_res.txt"), sep = "\t", as.is = T, header= T, check.names = F)
#
cpg_map_gene <- read.csv("cpg_map_gene.csv")
#gene name correction
cpg_map_gene$Gene[cpg_map_gene$Gene%in%c("BRUNOL5")] <- "CELF5"
cpg_map_gene <- cpg_map_gene[-which(cpg_map_gene$Gene%in%c("HLA-H", "NCRNA00092")), ]
#
record <- c()
for(i in 1:dim(cpg_map_gene)[1]){
	cpg <- cpg_map_gene$CpG[i]
	gene <- cpg_map_gene$Gene[i]
	cpg_value <- beta_res[cpg, ] %>% as.numeric()
	gene_value <- exp_mat[gene, colnames(beta_res)] %>% as.numeric()
	cor_info <- cor.test(cpg_value,gene_value, alternative = "two.sided", method = "pearson")
	p <- cor_info$p.value
	cor_value <- cor_info$estimate
	record <- rbind(record, c(cpg, gene, cor_value, p, cpg_map_gene$Group[i], cpg_map_gene$Island[i]))
}
colnames(record) <- c("cpg","gene","cor","p","group","island")
#
record <- as.data.frame(record)
record$gene_symbol <- record$gene
record$gene_symbol[record$gene_symbol%in%"CCDC109B"] <- "MCUB"
write.csv(record, file = "cpg_gene_cor.csv", row.names = F)
#
cpg_gene_cor <- read.csv("cpg_gene_cor.csv")
cpg_gene_cor$direction <- ifelse(cpg_gene_cor$cor > 0.1 & cpg_gene_cor$p < 0.05 , "Positive", ifelse(cpg_gene_cor$cor < -0.1 & cpg_gene_cor$p < 0.05, "Negative", "Not"))
pdf(paste0("Figure6E.pdf"), width = 7, height = 6)
picture <- volcano_dif_exp("", cpg_gene_cor)
print(picture)
dev.off()


