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
d[,3:ncol(d)] <- log10(d[,3:ncol(d)])

#getting rid of genes for which the minimum expression by any
# cell type is 1e-10
d <- d[ apply(d[,3:ncol(d)],1,max) > -500, ]

#Setting all values lower than 1e-10 to 1e-10
values <- d[,3:ncol(d)]
values[values < -500] <- -500
target <- rnorm(nrow(d))
values <- normalize.quantiles.use.target(as.matrix(values),target)
d[,3:ncol(d)] <- values

#==============================================
#Plotting Data
#==============================================
#Getting rid of gene descriptions for reformatting
# and plotting
plot_d <- d[-(1:2)]

plot_d <- melt(plot_d)
colnames(plot_d) <- c('cell_type','fpkm')
ggplot(plot_d, 
    aes(x=fpkm, fill=cell_type, group=cell_type)) + 
    geom_density(alpha=0.3) +
    facet_wrap(~ cell_type)

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

