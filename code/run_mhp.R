source("autoBlock.R")
load(file.path("data", "model_mhp.RData"))
dfmhp <- autoBlock(code, constants, data, inits, 50000, runList)$summary
save(dfmhp, file = file.path("results_hclust_average", "results_mhp.RData"))

