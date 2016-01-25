source("autoBlock.R")
load(file.path("data", "model_SSMcorrelated.RData"))
saveSamples <- TRUE
niter <- 50000
ab <- autoBlock(code, constants, data, inits, niter, runList, saveSamples = saveSamples)
dfSSMcorrelated <- ab$summary
save(dfSSMcorrelated, file = file.path("results_samples", "results_SSMcorrelated.RData"))
if (saveSamples) {
    burnedSamplesList <- ab$samples
    for (i in 1:length(burnedSamplesList)) burnedSamplesList[[i]] <- burnedSamplesList[[i]][(floor(niter/2) + 1):niter, ]
    save(burnedSamplesList, niter, file = file.path("results_samples", "results_SSMcorrelated_samples.RData"))
}

