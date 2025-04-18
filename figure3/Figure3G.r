library(trackViewer)
cpgs_info <- read.table("MEOX2_promoter_cpg_info.txt", header= F) #position information
beta_value <- read.table("MEOX2_meth_mat_cpgs.txt", header = T) #βvalue 
beta_value$dif <- beta_value$pos_value - beta_value$neg_value
beta_value <- beta_value[cpgs_info$V1, ]
cpg_pos <- cpgs_info$V3
sample.gr <- GRanges("chr7", IRanges(cpg_pos, width=1, names= cpgs_info$V1))
sample.gr$score <- beta_value$dif*10  #The score must exceed 1, so we multiply it by 10
features <- GRanges("chr7", IRanges(c(15726437, 15725511, 15666371, 15650837), 
                                    width=c(1500, 926, 172, 1399),
                                    names=c("TSS_Upstream","1stExon(5'UTR)", "Exon2","Exon3(3'UTR)")))
sample.gr$color <- sample.int(15, length(cpg_pos), replace = F)
sample.gr$border <- sample(c("gray80", "gray30"), length(cpg_pos), replace = T)
sample.gr$alpha <- sample(100:255, length(cpg_pos), replace = TRUE)/255

features$fill <- c("#FF8833", "#51C6E6", "#DFA32D", "#FFE4E1")
features$height <- c(0.04, 0.03, 0.02, 0.01)
gr <- GRanges("chr7", IRanges(15650837, 15727937, names="MEOX2"))
xaxis <- c(15650837, 15725511, 15726437, 15726637, 15727937)
rescale <- data.frame(
	from.start = c(15650837,15652236, 15666371,15666543, 15725511, 15726437),
	from.end = c(15652236, 15666371,15666543, 15725511,15726437, 15727937),
	to.start = c(15650837, 15658547, 15673967, 15677051,15681677, 15697097),
	to.end=c(15658547, 15673967,15677051, 15681677, 15697097,15727937)
)
pdf("Figure3G.pdf")
lolliplot(sample.gr, features, ranges=gr, rescale = rescale, xaxis = xaxis)
dev.off()
