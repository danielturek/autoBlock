

require(mcmcplots)
inputDir <- 'results_samples'
outputDir <- 'html_files'
fast <- TRUE


getSamplesFiles <- function(inputDir) {
    samplesFiles <- list.files(inputDir)
    samplesFiles <- samplesFiles[grepl('_samples.RData', samplesFiles)]
    samplesFiles <- file.path(inputDir, samplesFiles)
    return(samplesFiles)
}

## rename auto0, all, etc, to All Scalar, All Blocked, Auto Blocking
renameMCMCs <- function(burnedSamplesList) {
    finalAutoInd <- max(grep('^auto', names(burnedSamplesList)))
    finalAutoItName <- names(burnedSamplesList)[finalAutoInd]
    sampList <- list()
    sampList[['All Scalar']] <- burnedSamplesList[['auto0']]
    sampList[['All Blocked']] <- burnedSamplesList[['all']]
    sampList[['AutoBlock']] <- burnedSamplesList[[finalAutoItName]]
    return(sampList)
}


samplesFiles <- getSamplesFiles(inputDir) 
### TEMP REMOVE
print(samplesFiles)
samplesFiles <- samplesFiles[4:5]
file <- samplesFiles[1]
###

for(file in samplesFiles) {
    load(file)
    sampList <- renameMCMCs(burnedSamplesList)
    ## Optional: truncate all samples a lot!
    if(fast) {
        for(i in seq_along(sampList)) sampList[[i]] <- sampList[[i]][1:1000, ]
    }
    example_name <- gsub('_samples.RData', '', gsub(paste0(inputDir, '/', 'results_'), '', file))
    example_title <- switch(example_name,
                            test = 'Test model',
                            litters = 'Random Effects model'
                            )
    example_title_nospaces <- gsub(' ', '_', example_title)
    ## remove [],space from dimension names
    for(i in seq_along(sampList)) {
        dimnames(sampList[[i]])[[2]] <-
            gsub(perl=TRUE, '[\\] ]', '',
                 gsub(perl=TRUE, '[\\[,]', '_',
                      dimnames(sampList[[i]])[[2]]))
    }
    mcmcplot(sampList, style='plain', title=example_title, dir=outputDir, filename=example_title_nospaces)
    system(paste0('cd ', outputDir, '; pandoc -s -o ', example_title_nospaces, '.pdf ', example_title_nospaces, '.html'))
}




### TEMP REMOVE


