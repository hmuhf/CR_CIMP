type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
library(magrittr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
#color
low_value <- "#E7171B"
high_value <- "#1D4A9E"
setwd("figures/figure1")
#
select_probes <- read.table(paste0(condition,"/",type_output, "_selectProbes.txt"), header = F, as.is = T)
#
map_file <- read.table("cpg_map_gene_promoter_island.txt", header = T, as.is = T, sep = "\t") #promoter island cpgs with gene annotation information
promoter_info <- read.csv("promoter_info.csv", header =T, as.is = T)
select_probes <- select_probes$V1[select_probes$V1%in%promoter_info$IlmnID]
print(paste0("select promoter probes: ", length(select_probes)))
select_symbol <- map_file[map_file$ID%in%select_probes, "SYMBOL"] %>% unique() %>% na.omit()
write.table(select_symbol, file = paste0(condition, "_select_genes.txt"), col.names = F, row.names = F, quote = F)
select_genes <- map_file[map_file$ID%in%select_probes, "ENTREZ"] %>% unique() %>% na.omit()
print(paste0("select_probes: ", length(select_probes)))
print(paste0("select_genes:  ", length(select_genes)))
#GO
print("GO_BP: ")
ego_BP <- enrichGO(
gene          = select_genes,
OrgDb         = org.Hs.eg.db,
ont           = "BP",
pAdjustMethod = "BH",
pvalueCutoff  = 1,
qvalueCutoff = 1,
minGSSize = 10,
readable=TRUE)
print("batplot:")
pdf(file = paste(condition, "_BP_dotplot.pdf",sep = ""),width = 8,height = 5)
p <- dotplot(ego_BP, color = "pvalue", showCategory = 20, label_format = 60)+
scale_color_continuous(low=low_value,high=high_value)
print(p)
dev.off()
write.csv(ego_BP,file = paste(condition,"_go_BP.csv",sep = ""))
print("KEGG")
ekegg <- enrichKEGG(
gene          = select_genes,
keyType       = "kegg",
organism      = 'hsa',
pvalueCutoff  = 1,
qvalueCutoff  = 1,
pAdjustMethod = "BH",
minGSSize = 10)
pdf(file = paste("Figure1E.pdf",sep = ""),width = 8,height = 5)
p <- dotplot(ekegg, color = "pvalue", showCategory = 20, label_format = 60)+
scale_color_continuous(low=low_value, high=high_value)+ scale_size_continuous(range = c(3, 6))
print(p)
dev.off()
write.csv(ekegg ,file = paste(condition, "_kegg.csv",sep = ""))
