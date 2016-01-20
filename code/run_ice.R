source("autoBlock.R")
load(file.path("data", "model_ice.RData"))
dfice <- autoBlock(code, constants, data, inits, 50000, runList)$summary
save(dfice, file = file.path("results_hclust_wardd", "results_ice.RData"))

