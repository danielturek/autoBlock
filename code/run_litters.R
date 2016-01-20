source("autoBlock.R")
load(file.path("data", "model_litters.RData"))
dflitters <- autoBlock(code, constants, data, inits, 50000, runList)$summary
save(dflitters, file = file.path("results_hclust_median", "results_litters.RData"))

