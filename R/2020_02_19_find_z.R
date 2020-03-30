# By Borui Xiao
# directory: "~/zhou_lab/projects/20200202_TCGA_HM450/IDATs"
library(sesame)

#searchIDATprefixes("~/zhou_lab/projects/20200202_TCGA_HM450/IDATs") # get all samples' name in the folder(10102 pairs)
#lapply("~/zhou_lab/projects/20200202_TCGA_HM450/IDATs/9996247061_R06C02",readIDATpair) # return a list
#readIDATpair("~/zhou_lab/projects/20200202_TCGA_HM450/IDATs/9996247061_R06C02") # return a S4 class
load("~/zhou_lab/projects/20200202_TCGA_HM450/merged_mapping.rda")
#merged.mapping.normal[,3]

samples_name <- searchIDATprefixes("~/zhou_lab/projects/20200202_TCGA_HM450/IDATs")
sset_example <- lapply(samples_name[1],readIDATpair)
IG_example <- sset_example[[1]]@IG
result_list <- list()
result_matrix <- matrix(NA, nrow = dim(merged.mapping.normal)[1])

# newest code
sset <- lapply(samples_name[c(merged.mapping.normal[,3])],readIDATpair)

findz <- function(x){
  probes <- x@oobR
  c(rowSums(probes) - mean(rowSums(probes)))/sd(rowSums(probes))
}

result <- lapply(sset,findz)


result_matrix <- matrix(NA, nrow = dim(sset[[1]]@oobR)[1])

#time1 <- proc.time()
system.time(for (i in 1:749){
  result_matrix <- cbind(result_matrix,result[[i]])
})
#print(proc.time() - time1)
result_matrix <- result_matrix[,-1]
saveRDS(result_matrix, file = "z_probes_oobR.rds")


