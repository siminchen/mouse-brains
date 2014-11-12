library(ggplot2)
library(reshape2)

setwd('~/Development_Workspaces/mouse-brains/results/run1_1000genes/')
beta_weights <- read.csv('betas.txt',header=FALSE)
num_members <- read.csv('num_members.txt',header=FALSE)

sample_names <- c('Astrocytes1','Astrocytes2','Endothelial1','Endothelial2','Microglia1',
                  'Microglia2','MO1','MO2','Neuron1','Neuron2','NFO1','NFO2','OPC1','OPC2',
                  'WC1','WC2','WC3')

colnames(beta_weights) <- sample_names
colnames(num_members) <- sample_names

beta_weights$cluster <- as.factor(1:dim(beta_weights)[1])
num_members$cluster <- as.factor(1:dim(num_members)[1])

plot_beta <- melt(beta_weights)
colnames(plot_beta) <- c('cluster','cell_type','beta')
ggplot(plot_beta, aes(x=cluster,y=value, fill=cluster)) +
  geom_bar(stat='identity') +
  guides(fill=FALSE) +
  facet_wrap(~ variable) +
  ggtitle("Proportion of gene reads within each cluster (Beta)")

plot_mem <- melt(num_members)
colnames(plot_mem) <- c('cluster','cell_type','num_members')
ggplot(plot_mem, aes(x=cluster,y=num_members, fill=cluster)) +
  geom_bar(stat='identity') +
  guides(fill=FALSE) +
  facet_wrap(~ cell_type)

