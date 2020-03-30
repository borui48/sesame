setwd("~/zhoulab/labjournal/xiaob2")
setwd("/Volumes/zhoulab/labjournal/xiaob2")
z_IG <- readRDS("./z_probes/z_probes_IG.rds")
z_IR <- readRDS("./z_probes/z_probes_IR.rds")
z_II <- readRDS("./z_probes/z_probes_II.rds")

load("~/zhou_lab/projects/20200202_TCGA_HM450/merged_mapping.rda")
samples_name <- searchIDATprefixes("~/zhou_lab/projects/20200202_TCGA_HM450/IDATs")
sset_tumor_example <- lapply(samples_name[c(merged.mapping.tumor[,3][1])],readIDATpair)
oobG <- sset_tumor_example[[1]]@oobG
oobR <- sset_tumor_example[[1]]@oobR

# IG
IG <- sset_tumor_example[[1]]@IG
z_score_IG <- (rowSums(IG) - mean(rowSums(IG)))/sd(rowSums(IG))
z_score_IG_matrix <- as.matrix(z_score_IG)
z_score_IG_matrix <- cbind(rownames(z_score_IG_matrix),z_score_IG_matrix)

#calculate detection p-value of IG
p_value_IG <- function(z_x, z_probe, raw_x,ps,pf){
  pzs <- ecdf(z_probe[z_x[1],])(z_x[2])
  M <- ecdf(oobG)(raw_x[z_x[1],"M"])
  U <- ecdf(oobG)(raw_x[z_x[1],"U"])
  pzf <- 1 - apply(cbind(M,U),1,max)
  c(pzs*ps/(pzs*ps + pzf*pf))
}

# IR
IR<- sset_tumor_example[[1]]@IR
z_score_IR <- (rowSums(IR) - mean(rowSums(IR)))/sd(rowSums(IR))
z_score_IR_matrix <- as.matrix(z_score_IR)
z_score_IR_matrix <- cbind(rownames(z_score_IR_matrix),z_score_IR_matrix)

#calculate detection p-value of IR
p_value_IR <- function(z_x, z_probe, raw_x,ps,pf){
  pzs <- ecdf(z_probe[z_x[1],])(z_x[2])
  M <- ecdf(oobR)(raw_x[z_x[1],"M"])
  U <- ecdf(oobR)(raw_x[z_x[1],"U"])
  pzf <- 1 - apply(cbind(M,U),1,max)
  c(pzs*ps/(pzs*ps + pzf*pf))
}

# II
II <- sset_tumor_example[[1]]@II
z_score_II <- (rowSums(II) - mean(rowSums(II)))/sd(rowSums(II))
z_score_II_matrix <- as.matrix(z_score_II)
z_score_II_matrix <- cbind(rownames(z_score_II_matrix),z_score_II_matrix)

#calculate detection p-value of II
p_value_II <- function(z_x, z_probe, raw_x,ps,pf){
  pzs <- ecdf(z_probe[z_x[1],])(z_x[2])
  M <- ecdf(oobG)(raw_x[z_x[1],"M"])
  U <- ecdf(oobR)(raw_x[z_x[1],"U"])
  pzf <- 1 - apply(cbind(M,U),1,max)
  c(pzs*ps/(pzs*ps + pzf*pf))
}

# run test
system.time(result <- apply(z_score_IG_matrix, 1, p_value_IG, z_IG, IG, 0.5, 0.5))
system.time(result <- apply(z_score_IR_matrix, 1, p_value_IR, z_IR, IR, 0.5, 0.5))
system.time(result <- apply(z_score_II_matrix, 1, p_value_II, z_II, II, 0.5, 0.5))

# save result to rds files
saveRDS(result,file = "p_value_IG_ecdf.rds")
saveRDS(result,file = "p_value_IR_ecdf.rds")
saveRDS(result,file = "p_value_II_ecdf.rds")
#scp xiaob2@respublica.research.chop.edu:/home/xiaob2/zhoulab/labjournal/xiaob2/p_value_IG_ecdf.rds ~/Desktop/

# check if its right
IG_p_value[order(IG_p_value,decreasing=FALSE)[1:335]]
plot(density(IG_p_value))
pvalue <- readRDS("pvalue.rds")
