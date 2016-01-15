source("autoBlock.R")
load(file.path("data", "model_SSMcorrelated.RData"))
dfSSMcorrelated <- autoBlock(code, constants, data, inits, 5000, runList)$summary
save(dfSSMcorrelated, file = file.path("results_hclust_single", "results_SSMcorrelated.RData"))

