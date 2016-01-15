source("autoBlock.R")
load(file.path("data", "model_spatial.RData"))
dfspatial <- autoBlock(code, constants, data, inits, 5000, runList)$summary
save(dfspatial, file = file.path("results_hclust_single", "results_spatial.RData"))

