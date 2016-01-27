require(mcmcplots)
inputDir <- 'results_samples'
outputDir <- 'results_samples_plots'
outputFileName <- 'posterior_plots'
fast <- FALSE
verbose <- TRUE
openPDF <- TRUE

## read inputDir and get names of all *_samples.RData files
getSamplesFiles <- function(inputDir) {
    samplesFiles <- list.files(inputDir)
    samplesFiles <- samplesFiles[grepl('_samples.RData', samplesFiles)]
    samplesFiles <- file.path(inputDir, samplesFiles)
    return(samplesFiles)
}
makeAll <- function(samplesFiles, fast, verbose) {
    if(verbose) message('making list of all samples')
    ##if(fast) samplesFiles <- samplesFiles[3:4]
    all <- list()
    for(file in samplesFiles) {
        load(file)
        sampList <- renameMCMCs(burnedSamplesList)
        if(fast) for(i in seq_along(sampList)) sampList[[i]] <- sampList[[i]][1:1000, ]
        name <- getName(file)
        nameNS <- gsub(perl=TRUE, '[\\(\\) -]', '', name)
        all[[nameNS]] <- list(sampList=sampList, name=name, nameNS=nameNS)
    }
    return(all)
}
## change names auto0, all, etc, to All Scalar, All Blocked, Auto Blocking
renameMCMCs <- function(burnedSamplesList) {
    finalAutoInd <- max(grep('^auto', names(burnedSamplesList)))
    finalAutoItName <- names(burnedSamplesList)[finalAutoInd]
    sampList <- list()
    sampList[['All Scalar']] <- burnedSamplesList[['auto0']]
    sampList[['All Blocked']] <- burnedSamplesList[['all']]
    sampList[['AutoBlock']] <- burnedSamplesList[[finalAutoItName]]
    return(sampList)
}
## get the raw example name from filename
getName <- function(file) {
    name <- gsub(perl=TRUE, '.*/results_', '',
                 gsub(perl=TRUE, '_samples.RData', '', file))
    name <- switch(name, litters = 'Random Effects model', ice = 'Auto-Regressive model', SSMindependent = 'State Space model (independent)', SSMcorrelated = 'State Space model (correlated)', spatial = 'Spatial model', mhp = 'GLMM', test = 'Test model')
    return(name)
}
makeExamplePlots <- function(ex, outputDir, verbose) {
    if(verbose) message('making plot for ', ex$name)
    `All Scalar` <- ex$sampList[['All Scalar']]
    AutoBlock <- ex$sampList[['AutoBlock']]
    random <- switch(ex$nameNS,
                     RandomEffectsmodel = 14,
                     AutoRegressivemodel = NULL,
                     StateSpacemodelindependent = 20,
                     StateSpacemodelcorrelated = 20,
                     Spatialmodel = 20,
                     GLMM = 6)
    denoverplot(`All Scalar`, AutoBlock, style='plain', random=random, plot.title=ex$name,
                col = c('#009E73', 'blue'))
    dev.copy2pdf(file = paste0(outputDir, '/', ex$nameNS, '.pdf'))
    dev.off()
}
makePDF <- function(outputDir, outputFileName, verbose, openPDF) {
    if(verbose) message('building PDF file')
    system(paste0('cd ', outputDir, '; pdflatex ', outputFileName, '.tex'))
    system(paste0('rm ', outputDir, '/', outputFileName, '.aux'))
    system(paste0('rm ', outputDir, '/', outputFileName, '.log'))
    if(openPDF) system(paste0('open ', outputDir, '/', outputFileName, '.pdf'))
}

samplesFiles <- getSamplesFiles(inputDir)
all <- makeAll(samplesFiles, fast, verbose)
for(ex in all) makeExamplePlots(ex, outputDir, verbose)
makePDF(outputDir, outputFileName, verbose, openPDF)

##args(denoverplot)
##?denoverplot

##ex <- all[[2]]

##makeOutputFile(all, outputDir, outputFileName)
##if(openPDF) system(paste0('open ', outputDir, '/', outputFileName, '.pdf'))



## makeOutputFile <- function(all, outputDir, outputFileName) {
##     maxFloats <- 5; cur <- 0
##     outName <- paste0(outputDir, '/', outputFileName, '.tex')
##     system(paste0('rm -f ', outName))
##     system(paste0('touch ', outName))
##     con <- file(outName, open = 'wt')
##     writeLines('\\documentclass[12pt]{article}', con)
##     writeLines('\\usepackage{graphicx}', con)
##     writeLines('\\setlength{\\textheight}{9in}', con)
##     writeLines('\\setlength{\\topmargin}{-.5in}', con)
##     writeLines('\\setlength{\\headheight}{0.25in}', con)
##     writeLines('\\setlength{\\headsep}{0.25in}', con)
##     writeLines('\\setlength{\\topskip}{0in}', con)
##     writeLines('\\begin{document}', con)
##     writeLines('{\\LARGE Automated Blocking Posterior Plots} \\\\', con)
##     for(ex in all) {
##         writeLines(paste0('\\section{', ex$name, '}'), con)
##         nodes <- dimnames(ex$sampList[[1]])[[2]]
##         for(node in nodes) {
##             nodeForFilename <- gsub(perl=TRUE, '[\\[\\] ,]', '-', node)
##             dev.new(height=2, width=5)
##             par(mfrow = c(1,1), mar=c(3,3,2,1), mgp=c(0,0.6,0), tcl=-0.3)
##             ##  "All Scalar"  "All Blocked" "AutoBlock" 
##             tempArray <- cbind(ex$sampList[['All Scalar']][, node], ex$sampList[['All Blocked']][, node], ex$sampList[['AutoBlock']][, node])
##             densityList <- apply(tempArray, 2, density)
##             xlim <- range(unlist(lapply(densityList, function(d) d$x)))
##             xlim <- mean(xlim) + (xlim-mean(xlim)) * 1.5
##             ymax <- max(unlist(lapply(densityList, function(d) d$y))) * 1.1
##             plot(-100, -100, xlim=xlim, ylim=c(0,ymax), main=paste0('posterior density: ', node), xlab='', ylab='', yaxt='n', bty='n')
##             legend(x='topleft', legend=c('All Scalar', 'All Blocked', 'AutoBlock'), lty=1, lwd=2, col=c('green', 'red', 'blue'), bty='n')
##             for(i in 1:3) polygon(densityList[[i]], border=c('green', 'red', 'blue')[i])
##             abline(h=0, col='white')
##             ##plot(1, 1, main = node)
##             plotName <- paste0(ex$nameNS, '-', nodeForFilename)
##             dev.copy2pdf(file = paste0(outputDir, '/', plotName, '.pdf'))
##             dev.off()
##             writeLines('\\begin{figure}[ht]', con)
##             writeLines(paste0('\\includegraphics{', outputDir, '/', plotName, '}'), con)
##             writeLines('\\end{figure}', con)
##             cur <- cur+1;  if(cur==maxFloats) {cur <- 0;  writeLines('\\clearpage', con) }
##         }
##         writeLines('\\clearpage', con)
##     }
##     writeLines('\\end{document}', con)
##     close(con)
##     system(paste0('pdflatex -output-directory ', outputDir, ' ', outputFileName, '.tex'))
##     system(paste0('rm ', outputDir, '/', outputFileName, '.aux'))
##     system(paste0('rm ', outputDir, '/', outputFileName, '.log'))
## }
    
## remove { [ ] , space } from nodes (dimension) names
## renameNodes <- function(sampList) {
##     for(i in seq_along(sampList)) {
##         dimnames(sampList[[i]])[[2]] <-
##             gsub(perl=TRUE, '[\\] ]', '',
##                  gsub(perl=TRUE, '[\\[,]', '_',
##                       dimnames(sampList[[i]])[[2]])) }
##     return(sampList)
## }
    
