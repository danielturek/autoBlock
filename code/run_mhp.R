source("autoBlock.R")
load(file.path("data", "model_mhp.RData"))
dfmhp <- autoBlock(code, constants, data, inits, 2e+05, runList)$summary
save(dfmhp, file = file.path("results_hclust_single", "results_mhp.RData"))

