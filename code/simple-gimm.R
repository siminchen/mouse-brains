library(gimmR)
.GimmPath <- "~/gimmR/doc"
.PosthocPath <- "~/gimmR/doc"

data = read.table("data/formated_data.table")
data = data[1:100, ]
galGimm = runGimmNPosthoc(data, M=17, T=100, nreplicates=1, 
                          nIter=100, burnIn=50, verbose=T)

#jpeg("clustering-result.jpg")
heatmap(data.matrix(galGimm$clustData[, -(1:2)]), 
        Rowv=as.dendrogram(galGimm$hGClustData), Colv=NA, 
        labCol=colnames(galGimm$clustData)[-(1:2)], 
        labRow=galGimm$clustData[, 2], scale="none")
#dev.off()