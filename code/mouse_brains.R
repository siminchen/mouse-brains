library(ggplot2)
library(reshape2)
library(preprocessCore)

d <- read.table('formatted_data.table',header=T)

set.seed <- 071738561

#==============================================
#Preliminary QC - censoring & log transformation
#==============================================

#Setting gene names as row names
rownames(d) <- d[,1]

#log10 transformation
#not purely necessary, as the quantile normalization will
#squash the difference made here, but good for display purposes
d[,3:ncol(d)] <- log10(d[,3:ncol(d)])

#getting rid of genes for which the minimum expression by any
# cell type is 1e-10
d <- d[ apply(d[,3:ncol(d)],1,max) > -500, ]

#Quantile Normalization
values <- d[,3:ncol(d)]
normalized_values <- as.matrix(values)
#crappy code - but I'm not sure how to index an element out 
# of the result of an 'apply' function like this
for(i in 1:dim(values)[1]){
  normalized_values[i,] <- qqnorm(values[i,], plot.it=F)$x
}
nd <- d
nd[,3:ncol(nd)] <- normalized_values

#Setting low values of the log10 data to -1e10
#(doing this after normalization to not introduce false ties)
values[values < -10] <- -10
d[,3:ncol(d)] <- values
#==============================================
#Plotting Data
#==============================================
#Getting rid of gene descriptions for reformatting
# and plotting
plot_d <- d[-(1:2)]
plot_nd <- nd[-(1:2)]

plot_d <- melt(plot_d)
colnames(plot_d) <- c('cell_type','fpkm')
ggplot(plot_d, 
    aes(x=fpkm, fill=cell_type, group=cell_type)) + 
    geom_density(alpha=0.3) +
    facet_wrap(~ cell_type) +
    labs(title = "Before Quantile Normalization")

plot_nd <- melt(plot_nd)
colnames(plot_nd) <- c('cell_type','fpkm')
ggplot(plot_nd, 
    aes(x=fpkm, fill=cell_type, group=cell_type)) + 
    geom_density(alpha=0.3) +
    facet_wrap(~ cell_type) +
    labs(title = "After Quantile Normalization")
#==============================================
#QQ plots
#==============================================
#Selecting a column from the ORIGINAL data matrix
qqd <- read.table('formatted_data.table',header=T)
#log10 transformation
qqd[,3:ncol(d)] <- log10(qqd[,3:ncol(d)])

qqtest <- qqd$Astrocytes2
qqtest <- qqtest[is.finite(qqtest)]
qqtest <- sort(qqtest)

#Straight q-q plot(norm)
qqnorm <- rnorm(length(qqtest), mean(qqtest), sd(qqtest))
qqnorm <- sort(qqnorm)

plot(qqtest, qqnorm, main='Single Cell-Type QQplot (Normal)')
abline(mean(qqtest), 1)

#Censoring values below 1e-3
qqtest <- qqtest[qqtest > -3]
qqnorm <- rnorm(length(qqtest), mean(qqtest), sd(qqtest))
qqnorm <- sort(qqnorm)

plot(qqtest, qqnorm, main='Single Cell-Type QQplot w/ Censoring (Normal)')
abline(mean(qqtest), 1)
#==============================================
#Saving normalized data
#==============================================
write.table(nd, 'formatted_normalized_data.table')


#==============================================
#CODE BELOW NOT IN USE CURRENTLY
#==============================================



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

