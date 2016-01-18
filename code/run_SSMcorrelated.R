source("autoBlock.R")
load(file.path("data", "model_SSMcorrelated.RData"))
dfSSMcorrelated <- autoBlock(code, constants, data, inits, 50000, runList)$summary
save(dfSSMcorrelated, file = file.path("results_hclust_average", "results_SSMcorrelated.RData"))

