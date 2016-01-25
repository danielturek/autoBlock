source("autoBlock.R")
load(file.path("data", "model_mhp.RData"))
saveSamples <- TRUE
niter <- 50000
ab <- autoBlock(code, constants, data, inits, niter, runList, saveSamples = saveSamples)
dfmhp <- ab$summary
save(dfmhp, file = file.path("results_samples", "results_mhp.RData"))
if (saveSamples) {
    burnedSamplesList <- ab$samples
    for (i in 1:length(burnedSamplesList)) burnedSamplesList[[i]] <- burnedSamplesList[[i]][(floor(niter/2) + 1):niter, ]
    save(burnedSamplesList, niter, file = file.path("results_samples", "results_mhp_samples.RData"))
}

