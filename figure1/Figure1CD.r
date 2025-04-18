condition <- "Temozolomide_LGG"
setwd("figures/figure1")
anno_file <- read.csv("humanmethylation450_15017482_v1-2.csv", header =  T, as.is = T, check.names = F, na.strings = "")
probes <- read.table(paste0(condition,"/sdTop(1_per)_selectProbes.txt"), header = F, as.is = T)
probes <- probes$V1
#
gene_region<-c("TSS1500","TSS200","5'UTR","1stExon","Body","3'UTR","IGR")
island_region <- c("Island", "N_Shelf", "S_Shelf", "N_Shore", "S_Shore", "openSea")
#
anno_file$UCSC_RefGene_Group[is.na(anno_file$UCSC_RefGene_Group)] <- "IGR"
anno_file$Relation_to_UCSC_CpG_Island[is.na(anno_file$Relation_to_UCSC_CpG_Island)] <- "openSea"

cpg_select <- read.table("filtered_cpgs.txt",header = T,as.is = T)#filtered cpgs as background
anno_file <- anno_file[anno_file$IlmnID%in%cpg_select[,1], ]
probes_anno <- anno_file[anno_file$IlmnID%in%probes, ]
#
TSS1500_back <- dim(anno_file[grepl(gene_region[1], anno_file$UCSC_RefGene_Group), ])[1]
TSS200_back <- dim(anno_file[grepl(gene_region[2], anno_file$UCSC_RefGene_Group), ])[1]
UTR5_back <- dim(anno_file[grepl(gene_region[3], anno_file$UCSC_RefGene_Group), ])[1]
Exon1st_back <- dim(anno_file[grepl(gene_region[4], anno_file$UCSC_RefGene_Group), ])[1]
Body_back <- dim(anno_file[grepl(gene_region[5], anno_file$UCSC_RefGene_Group), ])[1]
UTR3_back <- dim(anno_file[grepl(gene_region[6], anno_file$UCSC_RefGene_Group), ])[1]
IGR_back <- dim(anno_file[grepl(gene_region[7], anno_file$UCSC_RefGene_Group), ])[1]
#
Island_back <- dim(anno_file[grepl(island_region[1], anno_file$Relation_to_UCSC_CpG_Island), ])[1]
N_Shelf_back <- dim(anno_file[grepl(island_region[2], anno_file$Relation_to_UCSC_CpG_Island), ])[1]
S_Shelf_back <- dim(anno_file[grepl(island_region[3], anno_file$Relation_to_UCSC_CpG_Island), ])[1]
N_Shore_back <- dim(anno_file[grepl(island_region[4], anno_file$Relation_to_UCSC_CpG_Island), ])[1]
S_Shore_back <- dim(anno_file[grepl(island_region[5], anno_file$Relation_to_UCSC_CpG_Island), ])[1]
openSea_back <- dim(anno_file[grepl(island_region[6], anno_file$Relation_to_UCSC_CpG_Island), ])[1]
#
#
TSS1500_num <- dim(probes_anno[grepl(gene_region[1], probes_anno$UCSC_RefGene_Group), ])[1]
TSS200_num <- dim(probes_anno[grepl(gene_region[2], probes_anno$UCSC_RefGene_Group), ])[1]
UTR5_num <- dim(probes_anno[grepl(gene_region[3], probes_anno$UCSC_RefGene_Group), ])[1]
Exon1st_num <- dim(probes_anno[grepl(gene_region[4], probes_anno$UCSC_RefGene_Group), ])[1]
Body_num <- dim(probes_anno[grepl(gene_region[5], probes_anno$UCSC_RefGene_Group), ])[1]
UTR3_num <- dim(probes_anno[grepl(gene_region[6], probes_anno$UCSC_RefGene_Group), ])[1]
IGR_num <- dim(probes_anno[grepl(gene_region[7], probes_anno$UCSC_RefGene_Group), ])[1]
#
Island <- dim(probes_anno[grepl(island_region[1], probes_anno$Relation_to_UCSC_CpG_Island), ])[1]
N_Shelf <- dim(probes_anno[grepl(island_region[2], probes_anno$Relation_to_UCSC_CpG_Island), ])[1]
S_Shelf <- dim(probes_anno[grepl(island_region[3], probes_anno$Relation_to_UCSC_CpG_Island), ])[1]
N_Shore <- dim(probes_anno[grepl(island_region[4], probes_anno$Relation_to_UCSC_CpG_Island), ])[1]
S_Shore <- dim(probes_anno[grepl(island_region[5], probes_anno$Relation_to_UCSC_CpG_Island), ])[1]
openSea <- dim(probes_anno[grepl(island_region[6], probes_anno$Relation_to_UCSC_CpG_Island), ])[1]
#
promoter_num <- TSS1500_num+TSS200_num+UTR5_num+Exon1st_num
Shelf <- N_Shelf+S_Shelf
Shore <- N_Shore+ S_Shore
promoter_back <- TSS1500_back+TSS200_back+UTR5_back+Exon1st_back
Shelf_back <- N_Shelf_back+S_Shelf_back
Shore_back <- N_Shore_back+ S_Shore_back
# count <- c(TSS1500_num, TSS200_num, UTR5_num, Exon1st_num,Body_num,UTR3_num,IGR_num,
           # TSS1500_back, TSS200_back, UTR5_back, Exon1st_back, Body_back,UTR3_back,IGR_back)
count <- c(promoter_num,Body_num,UTR3_num,IGR_num,
           promoter_back, Body_back,UTR3_back,IGR_back)
ratio <- c(promoter_num/dim(probes_anno)[1],Body_num/dim(probes_anno)[1], UTR3_num/dim(probes_anno)[1],
           IGR_num/dim(probes_anno)[1], promoter_back/dim(anno_file)[1], Body_back/dim(anno_file)[1],
           UTR3_back/dim(anno_file)[1], IGR_back/dim(anno_file)[1])

group <- c(rep("top sd 1%", 4), rep("background", 4))
variable <- factor(c("Promoter", gene_region[5:7]), levels = c("Promoter", gene_region[5:7]))
pdf(paste0("Figure1C.pdf"))
print(variable)
print(count)
#color<- c('#ef476f','#ffd166','#06d6a0','#118ab2')
color <- c('#118ab2', '#ef476f')
p <- count_bar("", group, variable, ratio, color)
print(p)
dev.off()
#
count <- c(Island,Shore,Shelf,openSea, Island_back,Shore_back,Shelf_back,openSea_back)
ratio <- c(Island/dim(probes_anno)[1], Shore/dim(probes_anno)[1], Shelf/dim(probes_anno)[1], openSea/dim(probes_anno)[1],
            Island_back/dim(anno_file)[1], Shore_back/dim(anno_file)[1], Shelf_back/dim(anno_file)[1],
            openSea_back/dim(anno_file)[1])

group <- c(rep("top sd 1%", 4), rep("background", 4))
variable <- factor(rep(c("Island","Shore","Shelf","openSea"), 2), levels  = c("Island","Shore","Shelf","openSea"))
#function plot
count_bar <- function(title_name, group, variable, count, man_color){
	library(ggpubr)
	pdat <- data.frame(group, variable, count)
	colour_man <- c("#3C83C4","#EF7B52")
	p1 <- ggbarplot(pdat, x = "variable", y = "count", fill = "group",
					width = 0.3, title = title_name, 
					# label = T, lab.col = "white", lab.pos = "in",
					label = F, 
					#position = position_fill())+
					position = position_dodge(width = 0.3))+ 
					xlab("")+
					ylab("ratio")+
					theme(text = element_text(size = 15), 
					plot.title = element_text(hjust=0.5))+
	  scale_fill_manual(values = man_color)
					
	return(p1)
}	
pdf(paste0("Figure1D.pdf"))
print(variable)
print(count)
p <- count_bar("", group, variable, ratio, color)
print(p)
dev.off()
#hypergeometric test
overlap <- promoter_num
back_fun <- promoter_back
back_all <- dim(anno_file)[1]
set_int <- dim(probes_anno)[1]
q<-overlap#
m<-back_fun# 
n<-back_all-m #
k<-set_int #
p<-phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p<-signif(p,2)
print(paste0("promoter test: ", p))
promoter_info <- c(set_int, back_all, back_fun, overlap, p)
#
overlap <- Body_num
back_fun <- Body_back
back_all <- dim(anno_file)[1]
set_int <- dim(probes_anno)[1]
q<-overlap# 
m<-back_fun# 
n<-back_all-m #
k<-set_int #
p<-phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p<-signif(p,2)
print(paste0("Body test: ", p))
Body_info <- c(set_int, back_all, back_fun, overlap, p)
#
overlap <- UTR3_num
back_fun <- UTR3_back
back_all <- dim(anno_file)[1]
set_int <- dim(probes_anno)[1]
q<-overlap# 
m<-back_fun# 
n<-back_all-m #
k<-set_int #
p<-phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p<-signif(p,2)
print(paste0("UTR3 test: ", p))
UTR3_info <- c(set_int, back_all, back_fun, overlap, p)
#
overlap <- IGR_num
back_fun <- IGR_back
back_all <- dim(anno_file)[1]
set_int <- dim(probes_anno)[1]
q<-overlap# 
m<-back_fun# 
n<-back_all-m #
k<-set_int #
p<-phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p<-signif(p,2)
print(paste0("IGR test: ", p))
IGR_info <- c(set_int, back_all, back_fun, overlap, p)
#
overlap <- Island
back_fun <- Island_back
back_all <- dim(anno_file)[1]
set_int <- dim(probes_anno)[1]
q<-overlap# 
m<-back_fun# 
n<-back_all-m #
k<-set_int #
p<-phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p<-signif(p,2)
print(paste0("island test: ", p))
island_info <- c(set_int, back_all, back_fun, overlap, p)
#
overlap <- Shore
back_fun <- Shore_back
back_all <- dim(anno_file)[1]
set_int <- dim(probes_anno)[1]
q<-overlap# 
m<-back_fun# 
n<-back_all-m #
k<-set_int #
p<-phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p<-signif(p,2)
print(paste0("shore test: ", p))
shore_info <- c(set_int, back_all, back_fun, overlap, p)
#
overlap <- Shelf
back_fun <- Shelf_back
back_all <- dim(anno_file)[1]
set_int <- dim(probes_anno)[1]
q<-overlap# 
m<-back_fun# 
n<-back_all-m #
k<-set_int #
p<-phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p<-signif(p,2)
print(paste0("shelf test: ", p))
shelf_info <- c(set_int, back_all, back_fun, overlap, p)
#
overlap <- openSea
back_fun <- openSea_back
back_all <- dim(anno_file)[1]
set_int <- dim(probes_anno)[1]
q<-overlap# 
m<-back_fun# 
n<-back_all-m #
k<-set_int #
p<-phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p<-signif(p,2)
print(paste0("openSea test: ", p))
opensea_info <- c(set_int, back_all, back_fun, overlap, p)
result <- rbind(promoter_info, Body_info, UTR3_info, IGR_info, island_info, shore_info, shelf_info, opensea_info)
colnames(result) <- c("interesting", "back", "back_fun", "overlap", "p")

#
write.table(result, file = paste0(condition, "_location_test.txt"), sep = "\t", quote = F, row.names = F, col.names = T)
