source("autoBlock.R")
load(file.path("data", "model_spatial.RData"))
saveSamples <- TRUE
niter <- 50000
ab <- autoBlock(code, constants, data, inits, niter, runList, saveSamples = saveSamples)
dfspatial <- ab$summary
save(dfspatial, file = file.path("results_samples", "results_spatial.RData"))
if (saveSamples) {
    burnedSamplesList <- ab$samples
    for (i in 1:length(burnedSamplesList)) burnedSamplesList[[i]] <- burnedSamplesList[[i]][(floor(niter/2) + 1):niter, ]
    save(burnedSamplesList, niter, file = file.path("results_samples", "results_spatial_samples.RData"))
}

