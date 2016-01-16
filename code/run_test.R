source("autoBlock.R")
load(file.path("data", "model_test.RData"))
dftest <- autoBlock(code, constants, data, inits, 50000, runList)$summary
save(dftest, file = file.path("results_hclust_single", "results_test.RData"))

