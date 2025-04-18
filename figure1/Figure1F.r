setwd("figures/figure1")
#
library(ggvenn)
total_gene <- read.table("v22_protein_coding_genesymbol.txt", header = F) #All protein-code genes
total_gene <- total_gene$V1
insgene=read.table("Temozolomide_LGG_select_genes.txt", header = F) #HVGs
insgene <- insgene$V1
#fungene=read.table("CGC_oncogene_symbol.txt",header = F) 
fungene=read.table("EMT_genesymbol.txt",header = F) 
fungene <- fungene$V1
print(length(insgene))
print(length(fungene))
#
ins_protein <- intersect(insgene, total_gene)
fun_protein <- intersect(fungene, total_gene)
print(length(ins_protein))
print(length(fun_protein))
#
q <- length(intersect(ins_protein, fun_protein))
m <- length(fun_protein)
n <- length(total_gene)-m
k <- length(ins_protein)
print(paste0("overlap genes：", q))
print(paste0("funSet genes: ", m))
print(paste0("all genes minus funSet genes: ", n))
print(paste0("insGenes: ", k))
p <- phyper(q-1, m, n, k,log.p = F,lower.tail = F)
p <- signif(p,2)
print(paste0("hyper.test: ", p))
venn_data<-list(
  ins_gene = ins_protein,
  fun_gene = fun_protein
)
p1<-ggvenn(venn_data, c("ins_gene", "fun_gene"))+ggtitle(paste0("hyper.test: ", p))
pdf("Figure1F.pdf")
print(p1)
dev.off()
#
