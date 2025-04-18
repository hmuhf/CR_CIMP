setwd("figures/figure3")
type_condition <- "drug_cancer" 
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
cancer <- "LGG"
CIMP_neg <- "C2"
CIMP_pos <- "C1"
library(magrittr)
library(RColorBrewer)
library(ggpubr)
library(ggplot2)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(ggrepel)
print(paste0("condition: ", condition))
drug <- unlist(strsplit(condition, split = "_"))[1]
exp_mat_fpkm <- read.table(paste0("expression-FPKM/Temozolomide/exp_mat.txt"), sep = "\t", as.is = T, header = T, na.strings = "", check.names = F)
exp_mat <- read.table(paste0("expression/", condition, "/exp_mat.txt"), sep = "\t", as.is = T, header = T, na.strings = "", check.names = F)
exp_mat <- 2^exp_mat - 1 #raw count
print(colnames(exp_mat) %>% head())
print(dim(exp_mat))
#
subtype <- read.table(paste(condition,"/", type_output,"_subtype of Cluster_",N,".txt",sep = ""),as.is = T,header= T)
print(rownames(subtype) %>% head())
exp_mat <- exp_mat[, colnames(exp_mat)%in%rownames(subtype)]#
print(dim(exp_mat))
cell_CIMP_neg <- intersect(rownames(subtype)[subtype$clusterCS==CIMP_neg], colnames(exp_mat))
cell_CIMP_pos <- intersect(rownames(subtype)[subtype$clusterCS==CIMP_pos], colnames(exp_mat))
print(length(cell_CIMP_neg))
print(length(cell_CIMP_pos))
exp_CIMP_neg <- exp_mat[, cell_CIMP_neg]
exp_CIMP_pos <- exp_mat[, cell_CIMP_pos]
exp_CIMP <- cbind(exp_CIMP_neg,exp_CIMP_pos)
#
print(dim(exp_CIMP))
print(head(exp_CIMP[,1:4]))
#2.edgeR
print("caculate the difference:")
print(dim(exp_CIMP_neg))
print(dim(exp_CIMP_pos))
#
group<-factor(c(rep("CIMP-",dim(exp_CIMP_neg)[2]),rep("CIMP+",dim(exp_CIMP_pos)[2])),levels = c("CIMP-","CIMP+"))
DGElist_data <- DGEList(counts = exp_CIMP, group = group)
keep <- filterByExpr(DGElist_data)#filter low expression genes.
y <- DGElist_data[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
design<-model.matrix(~group)
y <- estimateDisp(y, design)
#
fit <- glmFit(y,design) #
lrt <- glmLRT(fit, coef = 2) 
result <- topTags(lrt, n = nrow(DGElist_data$counts))
DMPs <- result$table
#
print("plot the volcano:")
logFC_cutoff <- 1
DMPs$change = as.factor(ifelse(DMPs$FDR < 0.05 &
								 abs(DMPs$logFC)>logFC_cutoff,
							   ifelse(DMPs$logFC>logFC_cutoff,'UP','DOWN'),'NOT'))
write.table(DMPs, file = paste(condition,"_DMPs_edgeR","_C",N,".txt",sep = ""), col.names = T, row.names = T, sep = "\t")
this_title<-paste0(condition,'\n','up genes: ',
				   nrow(DMPs[DMPs$change == 'UP',]),'; down genes: ',
				   nrow(DMPs[DMPs$change == 'DOWN',]))
pdf("Figure3F.pdf",width = 4,height = 4)
DMPs$symbol <- rownames(DMPs)
label_show <- DMPs[DMPs$symbol%in%c("STAC", "MEOX2"), ]
volcano = ggplot(data = DMPs,aes(x = logFC, y = -log10(FDR), color = change))+
  geom_point(alpha = 0.4,size = 1.75)+
theme(legend.background=element_rect(fill="transparent"),axis.line = element_line(color = "black"),
	panel.background=element_rect(fill="transparent"),
	plot.title = element_text(hjust = 0.5))+	  
xlab("Log2(FC)")+
  ylab("-log10(FDR)")+
  ggtitle(this_title)+
  theme(plot.title = element_text(size = 20,hjust = 0.5))+
  scale_colour_manual(values = c("#1D4A9E", "#B3B3B3", "#E7171B"))+geom_hline(yintercept = 1.3,linetype = "dotted",size = 1)+geom_vline(xintercept = c(-1,1),size=1,linetype = "dotted")+
  geom_point(size = 1.75, shape = 1, data = label_show) +
			geom_label_repel(
			aes(label = symbol),
			data = label_show,
			color="black")
print(volcano)
dev.off()
print("volcano plot complete!")
#
fpkmToTpm <- function(fpkm)
	{
		exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}
exp_mat <- 2^exp_mat_fpkm - 1
exp_mat <- apply(exp_mat, 2, fpkmToTpm)
exp_mat <- log(exp_mat+1, 2)
#
exp_CIMP_tpm <- exp_mat[, colnames(exp_CIMP)]
up_gene <- rownames(DMPs)[DMPs$change=="UP"]
down_gene <- rownames(DMPs)[DMPs$change=="DOWN"]
exp_tpm_heatmap <- exp_CIMP_tpm[c(up_gene, down_gene), ]

if(N==2){anno_col <- list(clusterCS=c(C1 = "#00BFC4",C2 = "#F8766D"))}
if(N==3){anno_col<-list(clusterCS=c(C1 = "#619CFF",C2 = "#00BA38",C3 = "#F8766D"))}
column_annotation <- HeatmapAnnotation(
		clusterCS = subtype[colnames(exp_tpm_heatmap), "clusterCS"],
		show_legend = F,
		col = anno_col
		)
row_annotation <- rowAnnotation(
		Direction = c(rep("up", length(up_gene)), rep("down", length(down_gene))),
		show_legend = F,
		col = anno_col
		)
pdat <- as.matrix(exp_tpm_heatmap)
pdat <- apply(exp_tpm_heatmap, 1, scale) %>% t()
rownames(pdat) <- rownames(exp_tpm_heatmap)
colnames(pdat) <- colnames(exp_tpm_heatmap)
result_hc <- Heatmap(pdat, cluster_columns = F, cluster_rows = F, 
					 show_row_names = F, show_column_names = F,
					 row_names_gp = gpar(fontsize = 7), 
					 top_annotation = column_annotation,
					 left_annotation = row_annotation,
					 column_title = type_output,
					 heatmap_legend_param = list(title = "z-score")
					 )
pdf("Figure3E.pdf", width = 8, height = 6)
print(result_hc)
dev.off()
# 5. 
print("differential expression:")
up_gene<-rownames(DMPs[DMPs$change=='UP',])
down_gene<-rownames(DMPs[DMPs$change=='DOWN',])
up_gene_id<-bitr(up_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
up_gene_id<-up_gene_id$ENTREZID
print(length(up_gene))
print(length(up_gene_id))
down_gene_id<-bitr(down_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
down_gene_id<-down_gene_id$ENTREZID
print(length(down_gene))
print(length(down_gene_id))
gene_id<-data.frame(gene=c(up_gene_id,down_gene_id),type = c(rep("up",length(up_gene_id)),rep("down",length(down_gene_id))))
case<-list(c("up"),c("down"), c("up", "down"))
for(i in 3){
	print(case[[i]])
	gene <- gene_id[gene_id$type%in%case[[i]],"gene"]
	print("GO_BP")
	ego_BP <- enrichGO(
	gene          = gene,
	OrgDb         = org.Hs.eg.db,
	ont           = "BP",
	pAdjustMethod = "BH",
	pvalueCutoff  = 1,
	qvalueCutoff = 0.05,
	minGSSize = 10,
	readable=TRUE)
	pdf(file = paste(condition,"_",paste0(case[[i]], collapse = "_"),"_BP_barplot_edgeR","_C",N,".pdf",sep = ""),width = 12,height = 6)
	p <- barplot(ego_BP,showCategory = 20,color = "p.adjust")
	print(p)
	dev.off()
	#
	pdf(file = paste(condition,"_",paste0(case[[i]], collapse = "_"),"_BP_dotplot_edgeR","_C",N,".pdf",sep = ""),width = 8,height = 8)
	p <- dotplot(ego_BP, color = "p.adjust", showCategory = 15)
	print(p)
	dev.off()
	write.csv(ego_BP,file = paste(condition,"_",paste0(case[[i]], collapse = "_"),"_function_go_BP_edgeR.csv",sep = ""))
	#KEGG
	print("KEGG")
	ekegg <- enrichKEGG(
	gene          = gene,
	keyType       = "kegg",
	organism      = 'hsa',
	pvalueCutoff  = 1,
	qvalueCutoff  = 0.05,
	pAdjustMethod = "BH",
	minGSSize = 10)
	if(dim(ekegg)[1]==0) {
		next
	}
	pdf(file = paste(condition,"_",paste0(case[[i]], collapse = "_"),"_kegg_barplot_edgeR.pdf",sep = ""),width = 12,height = 6)
	p <- barplot(ekegg,showCategory = 10,color = "p.adjust")
	print(p)
	dev.off()
	#
	low_value <- "#E7171B" 
	high_value <-"#1D4A9E"
	pdf("Figure3J.pdf",width = 8,height = 6)
	p <- dotplot(ekegg, color = "p.adjust", showCategory = 15, label_format = 60)+
		scale_color_continuous(low=low_value, high=high_value)+ scale_size_continuous(range = c(3, 6))
	print(p)
	dev.off()
	write.csv(ekegg ,file = paste(condition,"_", paste0(case[[i]], collapse = "_"), "_function_kegg_edgeR.csv",sep = ""))
}
