setwd("figures/figure3")
type_condition <- "drug_cancer" 
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
cancer <- "LGG"
CIMP_neg <- "C2"
CIMP_pos <- "C1"
library(maftools)
library(magrittr)
library(ggpubr)
fileList <- read.table(paste0("Mutation_MAF/",cancer,"/maf/fileList/MANIFEST.txt"), 
					 header = T, as.is =  T, sep = "\t")
fileList <- fileList[fileList$state=="validated", ]
#
if(cancer == "PAAD") {
skip_sample <- c(5,33,52,56,79,94,95,134,141,146)
}
if(cancer == "LGG"){
skip_sample <- c(118, 120, 194, 265, 364, 380)
}
if(cancer == "BLCA"){
skip_sample <- c(309)
}
if(cancer == "STAD"){
skip_sample <- c(238,344,405)
}
if(cancer == "LUAD"){
skip_sample <- c(150,549)
}
###PAAD: non synonymous-mutation: 5,33,52,56,79,94,95,134,141,146
# ###LGG:non synonymous-mutation: 118, 120, 194, 265, 364, 380
# ###BLCA: non synonymous-mutation: 309
subtype <- read.table(paste(type_output,"_subtype of Cluster_",N,".txt",sep = ""),as.is = T, header= T)
sample_annotation <- data.frame(sample = rownames(subtype), subtype = subtype$clusterCS)
if(N==2){
sample_annotation$subtype <- gsub(CIMP_neg, 'CIMP_neg', sample_annotation$subtype)
sample_annotation$subtype <- gsub(CIMP_pos, 'CIMP_pos', sample_annotation$subtype)
sample_annotation$subtype <- factor(sample_annotation$subtype, levels = c("CIMP_neg", "CIMP_pos"))
}
#
j <- 1
maf_list <- list()
for(i in 1:dim(fileList)[1]){
# print(i)
if(i%in%skip_sample) next()
maf_ind <- read.maf(paste0("Mutation_MAF/",cancer,"/maf/fileList/", fileList[i, 'filename']),
					verbose = F)
tumor_sample <- substring(maf_ind@data$Tumor_Sample_Barcode[1] %>% as.character(), 1, 16)
if(tumor_sample%in%rownames(subtype)){
  print(i)  
  maf_list[[j]] <- maf_ind
  j <- j+1
}
}
maf_all <- merge_mafs(maf_list, verbose = T)
maf_data <- maf_all@data
#

maf_sample <- data.frame(Tumor_Sample_Barcode = maf_data$Tumor_Sample_Barcode %>% as.character() %>% unique(),
					   code_pre = maf_data$Tumor_Sample_Barcode %>% as.character() %>% unique() %>% 
						 substring(., 1, 16))
sample_annotation <- merge(sample_annotation, maf_sample, by.x = 'sample', by.y = 'code_pre')
if(N==2){
subcolors <- c("#3C83C5","#EF7B52")
names(subcolors) <- c("CIMP_neg", "CIMP_pos")
}
annocolors <- list(subtype=subcolors)
vc_cols = RColorBrewer::brewer.pal(n = 7, name = 'Paired')
names(vc_cols) = c(
'Frame_Shift_Del',
'Missense_Mutation',
'Nonsense_Mutation',
'Multi_Hit',
'Frame_Shift_Ins',
'Splice_Site',
'In_Frame_Del'
)
pdf(file = paste0("Figure3D.pdf"),width = 10,height = 8)
oncoplot(maf_all, annotationDat = sample_annotation, clinicalFeatures = "subtype", sortByAnnotation = T,
	   annotationColor = annocolors, top = 10, colors = vc_cols)
dev.off()
##
tmb <- tmb(maf_all, logScale = FALSE)
write.csv(tmb, file = paste0(condition,"_tmb.csv"))
tmb <- merge(tmb, sample_annotation, by = 'Tumor_Sample_Barcode')
pdat <- tmb
print(head(pdat))
title <- "tmb"
#
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
if(N==2){
neg_value <- log(unlist(pdat[pdat$subtype=="CIMP_neg", "total_perMB"]) %>% as.numeric(), 2)
pos_value <- log(unlist(pdat[pdat$subtype=="CIMP_pos", "total_perMB"]) %>% as.numeric(), 2)
pl <- violin_plot("", neg_value, pos_value)
}	

if(N==2){
width_value <- 4
}
pdf(file = paste("Figure3A.pdf",sep = ""),width = width_value,height = 6)
print(pl)
dev.off()
