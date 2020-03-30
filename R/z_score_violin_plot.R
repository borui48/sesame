setwd("/Volumes/zhoulab/labjournal/xiaob2")
z_IG <- readRDS("./z_probes/z_probes_IG.rds")
library(matrixStats)
# plot one probe 
ggplot(NULL,aes(x = rownames(z_IG)[2],y= z_IG[2,]))+ 
  geom_violin()+
  theme_bw()+# white background+
  theme(
    panel.grid.major = element_blank(),# get rid of grid
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90, vjust=0.5)
    )+
  ylab("Within Sample Intensity z-score")+xlab("probe")

jpeg(file = "results_Value_2.jpg",width =1600,height = 2000,units = "px",res =300)
dev.off()


# show some of examples of probes
#orders <- order(rowMedians(z_IG))[1:20]
data_plot <- data.frame(
  probe = factor(rep(c(rownames(z_IG)[1:20]),each = 749)),
  intensity = c(z_IG[1,], z_IG[2,], z_IG[3,],z_IG[4,],z_IG[5,],z_IG[6,],
                z_IG[7,],z_IG[8,],z_IG[9,],z_IG[10,],z_IG[11,],z_IG[12,],z_IG[13,],
                z_IG[14,],z_IG[15,],z_IG[16,],z_IG[17,],z_IG[18,],z_IG[19,],z_IG[20,])
)

tail(data_plot)

ggplot(data_plot,aes(x = reorder(probe,intensity,median),y= intensity))+ 
  geom_violin()+
  theme_bw()+# white background+
  theme(
    panel.grid.major = element_blank(),# get rid of grid
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90, vjust=0.5)
  )+
  ylab("Within Sample Intensity z-score")+xlab("probe")+

# ecdf plot
e1 <-ecdf(z_IG[1,])
e2 <-ecdf(z_IG[3,])
e3 <-ecdf(z_IG[4,])
e4 <-ecdf(z_IG[5,])
e5 <- ecdf(z_IG[6,])

plot(e1,main="ECDF of z-score",xlab="z-score",ylab = "probability",col="red",xlim = c(-2,1),cex=0.15)
par(new=TRUE)
plot(e2,xlab = "", ylab = "",main="",col='green',xlim = c(-2,1),cex=0.15)
par(new=TRUE)
plot(e3,xlab = "", ylab = "",main="",col='yellow',xlim = c(-2,1),cex=0.15)
par(new=TRUE)
plot(e4,xlab = "", ylab = "",main="",col='black',xlim = c(-2,1),cex=0.15)
par(new=TRUE)
plot(e5,xlab = "", ylab = "",main="",col='brown',xlim = c(-2,1),cex=0.15)
legend("right",c(rownames(z_IG)[1],rownames(z_IG)[3],rownames(z_IG)[4],rownames(z_IG)[5],rownames(z_IG)[6]),
       fill = c("red","green","yellow","black","brown"),cex=0.5)
dev.off()
