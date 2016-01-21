source("autoBlock.R")
load(file.path("data", "model_SSMindependent.RData"))
dfSSMindependent <- autoBlock(code, constants, data, inits, 50000, runList)$summary
save(dfSSMindependent, file = file.path("results_hclust_complete2", "results_SSMindependent.RData"))

