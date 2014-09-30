library(ggplot2)
library(reshape2)

d <- read.csv('../data/barreslab_rnaseq.csv')

#Getting rid of gene description
d <- d[-2]
#Setting gene names as row names
rownames(d) <- d[,1]
d <- d[-1]

#Raw data density plot
d_formatted <- melt(d)
colnames(d_formatted) <- c('cell_type','fpkm')
ggplot(d_formatted, aes(x=fpkm, fill=cell_type)) + geom_density()

#log10 data density plot
d_log10 <- log10(d)
d_flog10 <- melt(d_log10)
colnames(d_flog10) <- c('cell_type','fpkm')
ggplot(d_flog10, aes(x=fpkm, fill=cell_type, group=cell_type)) + geom_density(alpha=0.3)

#Removing genes with a max 0.4 (lowest) value
d_nomin <- d[apply(d,1,max) > 0.4,]

#New log10 plot (no min values)
d_nomin_log10 <- log10(d_nomin)
d_fnlog10 <- melt(d_nomin_log10)
colnames(d_fnlog10) <- c('cell_type','fpkm')
#ggplot(d_fnlog10, aes(x=fpkm, fill=cell_type)) + geom_density(alpha=0.3)
ggplot(d_fnlog10, aes(x=fpkm, fill=cell_type)) + geom_histogram()

#Example of plotting log10 data of a single cell type
d_astrocyte <- melt(d$Astrocytes)
d10_astrocyte <- log10(d_astrocyte)
ggplot(d10_astrocyte, aes(x=value)) + geom_histogram()

#PRINCIPAL COMPONENT PLOTS (using no_min data)
pcs <- svd(d_nomin)
u_df <- data.frame(pcs$u)
ggplot(u_df,aes(x=X1,y=X2)) + geom_point()
ggplot(u_df,aes(x=X2,y=X3)) + geom_point()
ggplot(u_df,aes(x=X3,y=X4)) + geom_point()

#Heatmap of first 2000 genes
first2000 <- as.matrix(d_nomin[1:2000,])
heatmap(first2000)

#Function for making a qq plot from a (hard coded) distribution
my_qqplot <- function(x){
  m <- mean(x)
  sd <- sd(x)
  n <- length(x)
  
  probs <- (1:n)/(n+1)
  
  quants <- qchisq(probs, m)
  
  plot(sort(quants), sort(x))
  abline(0,1)
}
