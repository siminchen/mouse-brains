
setwd('~/Development_Workspaces/mouse-brains/data/')
data <- read.table('genes.fpkm_table', header=T)
samples <- read.table('sample_list')

data <- cbind(data[1] , 'dummy_id'=1:length(data[,1]),data[2:length(data)])
colnames(data) <- c('gene_id','dummy_id',as.vector(samples$V1))


