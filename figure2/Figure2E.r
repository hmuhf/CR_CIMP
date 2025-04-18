setwd("figures/figure2")
##
type_condition <- "drug_cancer"
condition <- "Temozolomide_LGG"
type_output <- "sdTop(1_per)"
N <- 2
drug_id <- "Temozolomide_1375"
GDSC <- "GDSC2"
CIMP_neg <- "C2"
CIMP_pos <- "C1"
library(magrittr)
#1 data preparing
print(paste0("condition: ", condition))
subtype <- read.table(paste0(condition,"/", type_output,"_subtype of Cluster_",N,".txt"), as.is = T, header= T)
#3 plot
library(ggpubr)
drugPredict <- read.csv(file = paste0(condition,"calcPhenotype_Output/DrugPredictions.csv"), check.names = F)

cell_name_list <- lapply(1:N, function(x){return(intersect(drugPredict[, 1], rownames(subtype)[subtype$clusterCS==paste0("C", x)]))})
cell_res_list <- lapply(1:N, function(x){return(drugPredict[drugPredict[, 1]%in%cell_name_list[[x]], drug_id] %>% as.numeric())})
#
pdat <- data.frame(response = unlist(cell_res_list), subtype = rep(paste0("C", 1:N), times = unlist(lapply(1:N, function(x){return(length(cell_res_list[[x]]))}))))
print("there")
if(N == 2){
	pdat$subtype <- gsub(CIMP_neg, 'CIMP_neg', pdat$subtype)
	pdat$subtype <- gsub(CIMP_pos, 'CIMP_pos', pdat$subtype)
	pdat$subtype <- factor(pdat$subtype, levels = c("CIMP_neg", "CIMP_pos"))
}
if(N == 2){
	colour_man <- c("#3C83C5","#EF7B52")
	neg_value <- pdat[pdat$subtype=="CIMP_neg", "response"] %>% as.numeric()
	pos_value <- pdat[pdat$subtype=="CIMP_pos", "response"] %>% as.numeric()
	p1 <- single_box("", neg_value, pos_value)
	p2 <- violin_plot("", neg_value, pos_value)


}
if(N == 2){
	width_value = 4
}
pdf(file = paste0(condition, "_prediction_", GDSC, "_C",N,".pdf"),width = width_value, height = 6)
print(p1)
dev.off()
pdf(file = paste0("Figure2E.pdf"),width = width_value, height = 6)
print(p2)
dev.off()	

