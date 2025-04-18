setwd("figures/figure2")
type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
CIMP_neg <- "C2"
CIMP_pos <- "C1"

library(ggplot2)
library(RColorBrewer)
library(magrittr)
# library(ComplexHeatmap)
library(ggpubr)
library(patchwork)
#
print("read the data:")
print(paste0("condition: ", condition))
drug_info <- unlist(strsplit(condition, split = "_"))
clinic <- read.table(paste0(condition, "/phenotype.txt"), header = T, as.is = T, sep = "\t")
colnames(clinic) <- c("id", "age", "gender", "race", "M_stage", "N_stage", "T_stage", "histologic_grade", "cancer_status", "Tumor_stage")
colnames(clinic) <- c("id", "age", "gender", "race", "M_stage", "N_stage", "T_stage", "histologic_grade", "cancer_status", "Tumor_stage")
subtype<-read.table(paste0(condition, "/", type_output,"_subtype of Cluster_",N,".txt"),as.is = T,header= T)
subtype_order <- read.table(paste0(condition, "/", type_output, "_colOrder of Cluster_", N, ".txt"), header = F, as.is = T)
#
print("clinic prepare:")
clinic_spe <- clinic[clinic$id%in%rownames(subtype), ]
rownames(clinic_spe) <- clinic_spe$id
clinic_spe <- clinic_spe[subtype_order$V1, ]
print(dim(clinic_spe))
#
clinic_spe$gender <- factor(clinic_spe$gender,levels = c("female","male"))
gender <- clinic_spe$gender
#
clinic_spe$race <- factor(clinic_spe$race, levels = c("asian", "black or african american", "white"))
race <- clinic_spe$race
#
clinic_spe$M_stage[grep("M1",clinic_spe$M_stage)] <- "M1"
clinic_spe$M_stage <- factor(clinic_spe$M_stage, levels = c("M0", "M1"))
M_stage <- clinic_spe$M_stage
clinic_spe$N_stage[grep("N1",clinic_spe$N_stage)] <- "N1"
clinic_spe$N_stage[grep("N2",clinic_spe$N_stage)] <- "N2"
clinic_spe$N_stage[grep("N3",clinic_spe$N_stage)] <- "N3"
clinic_spe$N_stage <- factor(clinic_spe$N_stage, levels = c("N0", "N1", "N2", "N3"))
N_stage <- clinic_spe$N_stage
clinic_spe$T_stage[grep("T1",clinic_spe$T_stage)] <- "T1"
clinic_spe$T_stage[grep("T2",clinic_spe$T_stage)] <- "T2"
clinic_spe$T_stage[grep("T3",clinic_spe$T_stage)] <- "T3"
clinic_spe$T_stage[grep("T4",clinic_spe$T_stage)] <- "T4"
clinic_spe$T_stage <- factor(clinic_spe$T_stage, levels = c("T0","T1", "T2", "T3", "T4"))
T_stage <- clinic_spe$T_stage
clinic_spe$Tumor_stage[clinic_spe$Tumor_stage %in% c("stage i", "stage ia", "stage ib", "stage ic")] <- "stage I"
clinic_spe$Tumor_stage[clinic_spe$Tumor_stage %in% c("stage ii", "stage iia", "stage iib", "stage iic")] <- "stage II"
clinic_spe$Tumor_stage[clinic_spe$Tumor_stage %in% c("stage iii", "stage iiia", "stage iiib", "stage iiic")] <- "stage III"
clinic_spe$Tumor_stage[clinic_spe$Tumor_stage %in% c("stage iv", "stage iva", "stage ivb", "stage ivc")] <- "stage IV"
clinic_spe$Tumor_stage <- factor(clinic_spe$Tumor_stage, levels = c("stage 0", "stage I", "stage II",
"stage III", "stage IV"))
Tumor_stage <- clinic_spe$Tumor_stage
print(table(Tumor_stage))
#
clinic_spe$histologic_grade <- factor(clinic_spe$histologic_grade, levels = c("G1", "G2", "G3", "G4"))
# clinic_spe$histologic_grade <- factor(clinic_spe$histologic_grade, levels = c("Low Grade", "High Grade"))
print(table(clinic_spe$histologic_grade))
clinic_spe$cancer_status <- factor(clinic_spe$cancer_status, levels = c("TUMOR FREE", "WITH TUMOR"))
##
##
clinic_list <- lapply(1:N, function(x){return(clinic_spe[rownames(subtype)[subtype$clusterCS == paste0("C", x)], ])})
##
feature_vector <- c("gender", "race", "M_stage", "N_stage", "T_stage", "histologic_grade", "cancer_status",
"Tumor_stage")
feature_value_list <- list(c("female","male"), c("asian", "black or african american", "white"), c("M0", "M1"), c("N0", "N1", "N2", "N3"), c("T0","T1", "T2", "T3", "T4"), c("G1", "G2", "G3", "G4"), c("TUMOR FREE", "WITH TUMOR"), c("stage 0", "stage I", "stage II", "stage III", "stage IV"))
# feature_value_list <- list(c("female","male"), c("asian", "black or african american", "white"), c("M0", "M1"), c("N0", "N1", "N2", "N3"), c("T0","T1", "T2", "T3", "T4"), c("Low Grade", "High Grade"), c("TUMOR FREE", "WITH TUMOR"), c("stage 0", "stage I", "stage II", "stage III", "stage IV"))
plot_list <- list()
for(i in 1:length(feature_vector)){
	print(paste0("feature: ", feature_vector[i]))
	print(paste0("feature value: ", paste0(feature_value_list[[i]], collapse = ",")))
	feature <- feature_vector[i]
	feature_value <- feature_value_list[[i]]
	feature_CIMP_list <- lapply(1:N, function(x){return(lapply(1:length(feature_value), function(y){return(which(clinic_list[[x]][, feature]==feature_value[y]) %>% length())}) %>% unlist())})
	print(feature_CIMP_list)
	
	mat <- matrix(unlist(feature_CIMP_list), length(feature_value), N, dimnames = list(feature = feature_value, subtype = paste0("C", 1:N)))
	print(mat)
	if(length(unique(as.numeric(mat))) == 1){
		if(unique(as.numeric(mat)) == 0){
			next
		}
	}
	p <- fisher.test(mat)$p.value %>% round(., 3)
	# p <- "ns"
	print(paste0("fisher.test: ", p))
	# p <- 1
	#visulization
	pdat <- data.frame(num = unlist(feature_CIMP_list), 
					   subtype = rep(paste0("C", 1:N), each = length(feature_value)),
					   feature = rep(feature_value, N))
	print(pdat)
	pdat <- pdat[which(pdat$num!=0), ]
	print(pdat)
	if(N == 2){
		pdat$subtype <- gsub(CIMP_neg, 'CIMP_neg', pdat$subtype)
		pdat$subtype <- gsub(CIMP_pos, 'CIMP_pos', pdat$subtype)
		pdat$subtype <- factor(pdat$subtype, levels = c("CIMP_neg", "CIMP_pos"))
	}
	title <- paste0(feature_vector[i], ": ", p) 
	if(length(unique(pdat$feature))==2){
	 pl <-ggbarplot(pdat, x = "subtype", y = "num", fill = "feature",
							   width = 0.3, title = title,
							   #label = T, lab.col = "white", lab.pos = "in",
							   position = position_fill())+
							   xlab("")+
							   ylab("frequency")+
							   theme(text = element_text(size = 15))+
							   scale_fill_manual(values = c("#118ab2", "#ef476f"))
			plot_list[[i]] <- pl
	}
	if(length(unique(pdat$feature))==3){
			 pl <-ggbarplot(pdat, x = "subtype", y = "num", fill = "feature",
							   width = 0.3, title = title,
							   #label = T, lab.col = "white", lab.pos = "in",
							   position = position_fill())+
							   xlab("")+
							   ylab("frequency")+
							   theme(text = element_text(size = 15))+
							   scale_fill_manual(values = c("#118ab2", "#06d6a0","#ef476f"))
			plot_list[[i]] <- pl
			}

	if(N == 2){
		width_value = 4
	}
	#Figure2B,C
	pdf(file = paste0(feature_vector[i], ".pdf"), width = width_value, height = 5)
	print(pl)
	dev.off()
}
##
age_list <- lapply(1:N, function(x){return(clinic_list[[x]][, "age"] %>% na.omit())})
pdat <- data.frame(age = unlist(age_list), subtype = rep(paste0("C", 1:	N), times = lapply(1:N, function(x){return(length(age_list[[x]]))}) %>% unlist()))
if(N == 2){
	pdat$subtype <- gsub(CIMP_neg, 'CIMP_neg', pdat$subtype)
	pdat$subtype <- gsub(CIMP_pos, 'CIMP_pos', pdat$subtype)
	pdat$subtype <- factor(pdat$subtype, levels = c("CIMP_neg", "CIMP_pos"))
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
if(N == 2){
	colour_man <- c("#00BFC4","#F8766D")
	neg_value <- pdat[pdat$subtype=="CIMP_neg", "age"] %>% as.numeric()
	pos_value <- pdat[pdat$subtype=="CIMP_pos", "age"] %>% as.numeric()
	p1 <- single_box("", neg_value, pos_value)
	p2 <- violin_plot("", neg_value, pos_value)
	pdf("age.pdf", width = 4, height = 6)
	print(p1)
	dev.off()
	pdf("Figure2D.pdf", width = 4, height = 6)
			print(p2)
			dev.off()
	plot_list[[length(feature_vector)+1]] <- p1
}
