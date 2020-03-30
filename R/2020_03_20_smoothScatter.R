#saveRDS(rowSums(sset_tumor_example[[1]]@IG),file ="intensity_IG_raw.rds")
#saveRDS(rowSums(sset_tumor_example[[1]]@IR),file ="intensity_IR_raw.rds")
#saveRDS(rowSums(sset_tumor_example[[1]]@II),file ="intensity_II_raw.rds")

setwd("/Volumes/zhoulab/labjournal/xiaob2")

intensity_IG_raw <- readRDS("raw_total_intensity/intensity_IG_raw.rds")
intensity_IR_raw <- readRDS("raw_total_intensity/intensity_IR_raw.rds")
intensity_II_raw <- readRDS("raw_total_intensity/intensity_II_raw.rds")

p_value_IG_gaussian <- readRDS("p_value_gaussian/p_value_IG_gaussian.rds")
p_value_IR_gaussian <- readRDS("p_value_gaussian/p_value_IR_gaussian.rds")
p_value_II_gaussian <- readRDS("p_value_gaussian/p_value_II_gaussian.rds")

p_value_IG_ecdf <- readRDS("p_value_ecdf/p_value_IG_ecdf.rds")
p_value_IR_ecdf <- readRDS("p_value_ecdf/p_value_IR_ecdf.rds")
p_value_II_ecdf <- readRDS("p_value_ecdf/p_value_II_ecdf.rds")

p_value_gaussian <- c(p_value_IG_gaussian,p_value_IR_gaussian,p_value_II_gaussian)
p_value_gaussian <- p_value_gaussian[order(names(p_value_gaussian))]

p_value_ecdf <- c(p_value_IG_ecdf,p_value_IR_ecdf,p_value_II_ecdf)
p_value_ecdf <- p_value_ecdf[order(names(p_value_ecdf))]

p_value_original <- 1 - readRDS("pvalue.rds")@pval
p_value_IG_original <- p_value_original[names(p_value_IG_gaussian)]
p_value_IR_original <- p_value_original[names(p_value_IR_gaussian)]
p_value_II_original <- p_value_original[names(p_value_II_gaussian)]

plot(p_value_gaussian, p_value_original, col = c("red", "blue"))
plot(p_value_gaussian, p_value_ecdf)
plot(p_value_ecdf, p_value_original)

smoothScatter(p_value_original,p_value_gaussian,nbin=128)
smoothScatter(p_value_gaussian,p_value_ecdf,nbin=64)
smoothScatter(p_value_IG_ecdf,p_value_IG_original)
smoothScatter(p_value_IR_ecdf,p_value_IR_original)
smoothScatter(p_value_II_ecdf,p_value_II_original)

# IG, IR, II ecdf pvalue vs IG, IR, II  original pvalue
plot(p_value_IG_ecdf,p_value_IG_original)
plot(p_value_IR_ecdf,p_value_IR_original)
plot(p_value_II_ecdf,p_value_II_original)

#IG, IR, II raw total intensity vs IG, IR, II ORIGINAL p value
par(mfrow=c(1,3))
plot(intensity_IG_raw,p_value_IG_original,xlab = "IG raw total intensity", ylab = "IG p value")
plot(intensity_IR_raw,p_value_IR_original,xlab = "IR raw total intensity", ylab = "IR p value")
plot(intensity_II_raw,p_value_II_original,xlab = "II raw total intensity", ylab = "II p value")
#IG, IR, II raw total intensity vs IG, IR, II ECDF p value
par(mfrow=c(1,3))
plot(intensity_IG_raw,p_value_IG_ecdf,xlab = "IG raw total intensity", ylab = "IG p value of e.c.d.f")
plot(intensity_IR_raw,p_value_IR_ecdf,xlab = "IR raw total intensity", ylab = "IR p value of e.c.d.f")
plot(intensity_II_raw,p_value_II_ecdf,xlab = "II raw total intensity", ylab = "II p value of e.c.d.f")
#IG, IR, II raw total intensity vs IG, IR, II GAUSSIAN p value
par(mfrow=c(1,3))
plot(intensity_IG_raw,p_value_IG_gaussian,xlab = "IG raw total intensity", ylab = "IG p value of gaussian distribution")
plot(intensity_IR_raw,p_value_IR_gaussian,xlab = "IR raw total intensity", ylab = "IR p value of gaussian distribution")
plot(intensity_II_raw,p_value_II_gaussian,xlab = "II raw total intensity", ylab = "II p value of gaussian distribution")





dev.off()
