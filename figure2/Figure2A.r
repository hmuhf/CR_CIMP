setwd("figures/figure2")
type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
CIMP_neg <- "C2"
CIMP_pos <- "C1"
library(RColorBrewer)
library(magrittr)
library(survival)
library(survminer)
#
print("read the data:")
print(paste0("condition: ", condition))
drug_info <- unlist(strsplit(condition, split = "_"))
#
survival_info <- read.table(paste0(condition, "/survival.txt"), header = T, as.is = T, sep = "\t")
rownames(survival_info) <- survival_info$sampleID
subtype <- read.table(paste0(condition, "/", type_output,"_subtype of Cluster_",N,".txt"),as.is = T,header= T)
subtype$ID <- rownames(subtype)
subtype_order <- read.table(paste0(condition, "/", type_output, "_colOrder of Cluster_", N, ".txt"), header = F, as.is = T)
#
survival_spe <- survival_info[subtype_order$V1, ]
print(survival_spe)
print(head(subtype))
survival_spe <- merge(survival_spe, subtype, by.x = "sampleID", by.y = "ID", all.x = T, all.y = F)
print(head(survival_spe))
print(tail(survival_spe))
if(N == 2){
	survival_spe$clusterCS <- gsub(CIMP_neg, 'CIMP_neg', survival_spe$clusterCS)
	survival_spe$clusterCS <- gsub(CIMP_pos, 'CIMP_pos', survival_spe$clusterCS)
	survival_spe$clusterCS <- factor(survival_spe$clusterCS, levels = c("CIMP_neg", "CIMP_pos"))
}
#
fit <- survfit(Surv(OS.time, OS) ~ clusterCS, data = survival_spe)
print(fit)
if(N == 2){
	colour_man <- c("#3C83C5","#EF7B52")
}
pl <-  ggsurvplot(
	   data = survival_spe,
	   fit,
	   pval = TRUE, 
	   conf.int = TRUE,
	   risk.table = TRUE, 
	   risk.table.col = "strata", 
	   linetype = "strata", 
	   surv.median.line = "hv",
	   ggtheme = theme_bw(),
	   palette = colour_man,
	   font.main = 18,    
	   font.x = 15,      
	   font.y = 15,     
	   font.tickslab = 15
	   )
pdf(file = paste("Figure2A.pdf",sep = ""),width = 8,height = 6)
print(pl)
dev.off()
