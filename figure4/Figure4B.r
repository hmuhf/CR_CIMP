setwd("figures/figure4")
##radarchart
library(ggradar)
library(dplyr)
library(scales)
library(tibble)
condition <- "Temozolomide_LGG"
CIMP_neg <- "C2"
CIMP_pos <- "C1"
type_condition <- "drug_cancer"
type_output <- "sdTop(1_per)"
N<-2
gene <- "CD274"
#
immune_score <- read.table("TCGA/TCGA_signature_score.txt", header = T, sep = "\t", na.strings = "NA", fill = T, quote = "")
print(colnames(immune_score))
immune_score$TCGA.Participant.Barcode <- paste0(immune_score$TCGA.Participant.Barcode, "-01A")
subtype <- read.table(paste(condition,"/",type_output,"_subtype of Cluster_",N,".txt",sep = ""),as.is = T, header= T)
subtype$sample <- rownames(subtype)
subtype$clusterCS[subtype$clusterCS==CIMP_neg] <- "CIMP_neg"
subtype$clusterCS[subtype$clusterCS==CIMP_pos] <- "CIMP_pos"
#TMB

#neoantigens
value_neg <- immune_score[immune_score$TCGA.Participant.Barcode%in%subtype[subtype$clusterCS=="CIMP_neg", "sample"], "SNV.Neoantigens"] %>% na.omit()
value_pos <- immune_score[immune_score$TCGA.Participant.Barcode%in%subtype[subtype$clusterCS=="CIMP_pos", "sample"], "SNV.Neoantigens"] %>% na.omit()
value <- c(value_neg, value_pos)
value <- rescale(value)
neg_len <- length(value_neg)
neoantigen <- c(mean(value[1:neg_len]), mean(value[(neg_len+1):length(value)]))
#lymphoctype infiltration
value_neg <- immune_score[immune_score$TCGA.Participant.Barcode%in%subtype[subtype$clusterCS=="CIMP_neg", "sample"], "Lymphocyte.Infiltration.Signature.Score"] %>% na.omit()
value_pos <- immune_score[immune_score$TCGA.Participant.Barcode%in%subtype[subtype$clusterCS=="CIMP_pos", "sample"], "Lymphocyte.Infiltration.Signature.Score"] %>% na.omit()
value <- c(value_neg, value_pos)
value <- rescale(value)
neg_len <- length(value_neg)
lymphoctype <- c(mean(value[1:neg_len]), mean(value[(neg_len+1):length(value)]))
#leukocyte fraction
value_neg <- immune_score[immune_score$TCGA.Participant.Barcode%in%subtype[subtype$clusterCS=="CIMP_neg", "sample"], "Leukocyte.Fraction"] %>% na.omit()
value_pos <- immune_score[immune_score$TCGA.Participant.Barcode%in%subtype[subtype$clusterCS=="CIMP_pos", "sample"], "Leukocyte.Fraction"] %>% na.omit()
value <- c(value_neg, value_pos)
value <- rescale(value)
neg_len <- length(value_neg)
leukocyte <- c(mean(value[1:neg_len]), mean(value[(neg_len+1):length(value)]))
#CD274
exp_mat <- read.table(paste0("/expression-FPKM/", condition, "/exp_mat.txt"), sep = "\t", as.is = T, header = T, na.strings = "", check.names = F)
#FPKM to TPM
fpkmToTpm <- function(fpkm)
	{
		exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}
exp_mat <- 2^exp_mat - 1#
exp_mat <- apply(exp_mat, 2, fpkmToTpm)
exp_mat <- log(exp_mat+1, 2)
cell_CIMP_neg <- intersect(rownames(subtype)[subtype$clusterCS=="CIMP_neg"], colnames(exp_mat))
cell_CIMP_pos <- intersect(rownames(subtype)[subtype$clusterCS=="CIMP_pos"], colnames(exp_mat))
exp_CIMP_neg <- exp_mat[, cell_CIMP_neg]
exp_CIMP_pos <- exp_mat[, cell_CIMP_pos]
gene_name <- gene
value_neg <- exp_CIMP_neg[gene_name, ] %>% as.numeric()
value_pos <- exp_CIMP_pos[gene_name, ] %>% as.numeric()
print(c(mean(value_neg), mean(value_pos)))
#print(value)
value <- c(value_neg, value_pos)
value <- rescale(value)
neg_len <- length(value_neg)
CD274 <- c(mean(value[1:neg_len]), mean(value[(neg_len+1):length(value)]))

#MHC score
#MHC_score
gene_set <- c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'NLRC5', 'PSMB9', 'PSMB8', 'B2M')
gene_set_mat <- exp_mat[gene_set, ]
mhc_score <- apply(gene_set_mat, 2, mean) %>% as.numeric()
names(mhc_score) <- colnames(exp_mat)
value_neg <- mhc_score[names(mhc_score)%in%cell_CIMP_neg]
value_pos <- mhc_score[names(mhc_score)%in%cell_CIMP_pos]
value <- c(value_neg, value_pos)
value <- rescale(value)
neg_len <- length(value_neg)
mhc <- c(mean(value[1:neg_len]), mean(value[(neg_len+1):length(value)]))
#cyt score
gene_set <- c("GZMA", "PRF1")
gene_set_mat <- exp_mat[gene_set, ]
cyt_score <- apply(gene_set_mat, 2, geometry_mean) %>% as.numeric()
names(cyt_score) <- colnames(exp_mat)
value_neg <- cyt_score[names(cyt_score)%in%cell_CIMP_neg]
value_pos <- cyt_score[names(cyt_score)%in%cell_CIMP_pos]
value <- c(value_neg, value_pos)
value <- rescale(value)
neg_len <- length(value_neg)
cyt <- c(mean(value[1:neg_len]), mean(value[(neg_len+1):length(value)]))
#TMB
tmb <- read.csv("Temozolomide_LGG_tmb.csv")
print(head(tmb))
tmb$SAMPLE_ID <- substring(tmb$Tumor_Sample_Barcode, 1, 16)
value_neg <- tmb[tmb$SAMPLE_ID%in%cell_CIMP_neg, "total_perMB"] %>% na.omit()
value_pos <- tmb[tmb$SAMPLE_ID%in%cell_CIMP_pos, "total_perMB"] %>% na.omit()
value <- c(value_neg, value_pos)
value <- rescale(value)
neg_len <- length(value_neg)
tmb <- c(mean(value[1:neg_len]), mean(value[(neg_len+1):length(value)]))

#
print("here")
data <- data.frame(neoantigen, lymphoctype, leukocyte, CD274, mhc, cyt, tmb)
data <- cbind(c("CIMP-","CIMP+"),data)
colnames(data) <- c("group", "neoantigen", "lymphoctype", "leukocyte", "CD274", "mhc", "cyt", "tmb")
write.table(data, file = "radar_dat_rescale.txt", sep = "\t", quote = F, col.names = T, row.names = F)
pdf("Figure4B.pdf")
ggradar(data,values.radar=c("0","0.5","1"), background.circle.colour="white",group.colours=c("#3C83C5","#EF7B52"))
dev.off()
#
