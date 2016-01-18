source("autoBlock.R")
load(file.path("data", "model_SSMindependent.RData"))
dfSSMindependent <- autoBlock(code, constants, data, inits, 50000, runList)$summary
save(dfSSMindependent, file = file.path("results_hclust_average", "results_SSMindependent.RData"))

