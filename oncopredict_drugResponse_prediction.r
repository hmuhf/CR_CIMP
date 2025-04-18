#oncopredict predicts drug response
drugRes_predict_stat <- function(type_condition, condition, type_output, N, CIMP_neg, CIMP_pos, GDSC){
  #type_condition <- "drug_cancer"
  #condition <- "Temozolomide_LGG"
  #type_output <- "sdTop(1_per)"
  #N <- 2
  #GDSC <- "GDSC2"
  #CIMP_neg <- "C2"
  #CIMP_pos <- "C1"
  library(magrittr)
  #1 data preparing
  print(paste0("condition: ", condition))
  drug <- unlist(strsplit(condition, split = "_"))[1]
  data_select <- "promoter_island"
  exp_mat <- read.table(paste0(condition, "/exp_mat.txt"), sep = "\t", as.is = T, header = T, na.strings = "", check.names = F)
  #FPKM to tpm
  fpkmToTpm <- function(fpkm)
	{
		exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	}
  exp_mat <- 2^exp_mat -1
  exp_mat <- apply(exp_mat, 2, fpkmToTpm)
  exp_mat <- log(exp_mat+1, 2)
  subtype <- read.table(paste0(type_output,"_subtype of Cluster_",N,".txt"), as.is = T, header= T)
  print(rownames(subtype) %>% head())
  exp_CIMP <- exp_mat[, colnames(exp_mat)%in%rownames(subtype)]
  #2 oncoPredict predicting
  library(oncoPredict)
  trainPhe <- readRDS(paste0("/Training\ Data/",GDSC,"_Res.rds"))
  trainExp <- readRDS(paste0("Training\ Data/",GDSC,"_Expr (RMA Normalized and Log Transformed).rds"))
  testExp <- exp_CIMP
  testExp <- as.matrix(testExp)
  result <- calcPhenotype(trainingExprData = trainExp,
                          trainingPtype = trainPhe,
                          testExprData = testExp,
                          batchCorrect = "eb",
                          minNumSamples = 10,
                          printOutput = F,
                          removeLowVaringGenesFrom = "homogenizeData"
						 )
}
