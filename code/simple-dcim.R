library(gimmR)
.GimmPath <- "~/gimmR/doc"
.PosthocPath <- "~/gimmR/doc"

data = read.table("data/formated_data.table")

dceGimm = runGimmNPosthoc(data, M=17, T=dim(data)[1], nIter=100,
                          burnIn=50, estimate_contexts="y", verbose=T
                          intFiles=T)

drawHeatmap(dceGimm, scale="none")