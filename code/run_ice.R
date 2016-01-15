source("autoBlock.R")
load(file.path("data", "model_ice.RData"))
dfice <- autoBlock(code, constants, data, inits, 5000, runList)$summary
save(dfice, file = file.path("results_hclust_single", "results_ice.RData"))

