

autoBlockModel <- setRefClass(
    Class = 'autoBlockModel',
    fields = list(
        code = 'ANY',
        constants = 'list',
        data = 'list',
        inits = 'list',
        md = 'ANY',
        Rmodel = 'ANY',
        scalarNodeVector = 'character',
        nodeGroupScalars = 'list',
        nodeGroupAllBlocked = 'list',
        nodeGroupStochNodes = 'list',
        monitorsVector = 'character'
        ),
    methods = list(
        initialize = function(code, constants, data, inits) {
            library(nimble)
            code <<- code
            constants <<- if(missing(constants)) list() else constants
            data <<- if(missing(data)) list() else data
            inits <<- if(missing(inits)) list() else inits
            md <<- nimbleModel(code=code, constants=constants, returnDef=TRUE)
            Rmodel <<- md$newModel(data=data, inits=inits)
            scalarNodeVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=TRUE)
            nodeGroupScalars <<- lapply(scalarNodeVector, function(x) x)
            nodeGroupAllBlocked <<- list(scalarNodeVector)
            stochNodeVector <- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE, returnScalarComponents=FALSE)
            nodeGroupStochNodes <<- lapply(stochNodeVector, function(x) x)
            monitorsVector <<- Rmodel$getNodeNames(stochOnly=TRUE, includeData=FALSE)
        },
        createGroups = function(listOfBlocks = list()) {
            listOfBlocks <- lapply(listOfBlocks, function(blk) Rmodel$expandNodeNames(blk, returnScalarComponents=TRUE))
            nodes <- scalarNodeVector
            nodes <- setdiff(nodes, unlist(listOfBlocks))
            nodeList <- lapply(nodes, function(x) x)
            for(ng in listOfBlocks) nodeList[[length(nodeList)+1]] <- ng
            return(nodeList)
        },
        newModel = function() {
            newRmodel <- md$newModel(data=data, inits=inits)
            return(newRmodel)
        }
        )
    )

createCodeAndConstants <- function(N, listOfBlockIndexes, rhoVector) {
    code <- quote({})
    constants <- list()
    if(length(listOfBlockIndexes) != length(rhoVector)) stop()
    for(i in seq_along(listOfBlockIndexes)) {
        blockIndexes <- listOfBlockIndexes[[i]]
        rho <- rhoVector[i]
        numNodes <- as.numeric(length(blockIndexes))
        if(numNodes == 1) {
            code[[length(code)+1]] <- substitute(x[IND] ~ dnorm(0, 1), list(IND=as.numeric(blockIndexes)))
        } else {
            muText <- paste0('mu', i)
            sigmaText <- paste0('Sigma', i)
            indMin <- as.numeric(min(blockIndexes))
            indMax <- as.numeric(max(blockIndexes))
            code[[length(code)+1]] <- substitute(x[MIN:MAX] ~ dmnorm(MU[1:NUM], cov = SIGMA[1:NUM,1:NUM]), list(MIN=indMin, MAX=indMax, NUM=numNodes, MU=as.name(muText), SIGMA=as.name(sigmaText)))
            constants[[muText]] <- rep(0, numNodes)
            constants[[sigmaText]] <- createCov(N=numNodes, rho=rho)
        }
    }
    allInd <- 1:N
    leftoverInd <- setdiff(allInd, unlist(listOfBlockIndexes))
    for(ind in leftoverInd) code[[length(code)+1]] <- substitute(x[IND] ~ dnorm(0, 1), list(IND=as.numeric(ind)))
    codeAndConstantsList <- list(code=code, constants=constants)
    return(codeAndConstantsList)
}

createCov <- function(N, indList=list(1:N), rho=0.8, indList2=list(), rho2=0.3, indList3=list(), rho3=0.5) {
    Sigma <- diag(N)
    for(gp in indList)  { for(i1 in gp) for(i2 in gp) Sigma[i1,i2] <- Sigma[i2,i1] <- rho  }
    for(gp in indList2) { for(i1 in gp) for(i2 in gp) Sigma[i1,i2] <- Sigma[i2,i1] <- rho2 }
    for(gp in indList3) { for(i1 in gp) for(i2 in gp) Sigma[i1,i2] <- Sigma[i2,i1] <- rho3 }
    diag(Sigma) <- 1
    Sigma
}


autoBlockParamDefaults <- function() {
    list(
        adaptIntervalBlock = 200,
        cutree_heights = seq(0, 1, by=0.1),
        makePlots = TRUE,
        niter = 200000,
        saveSamples = FALSE,
        setSeed0 = FALSE,
        verbose = TRUE
        )
}


autoBlock <- setRefClass(

    Class = 'autoBlock',

    fields = list(
        
        ## special
        abModel = 'ANY',
        it = 'numeric',

        ## overall control
        adaptIntervalBlock = 'numeric',
        cutree_heights = 'numeric',
        makePlots = 'logical',
        niter = 'numeric',
        saveSamples = 'logical',
        setSeed0 = 'logical',
        verbose = 'logical',

        ## persistant lists of historical data
        naming = 'list',
        candidateGroups = 'list',
        grouping = 'list',
        groupSizes = 'list',
        groupIDs = 'list',
        samplers = 'list',
        Cmcmcs = 'list',
        timing = 'list',
        samples = 'list',
        means = 'list',
        sds = 'list',
        ess = 'list',
        essPT = 'list',
        burnedSamples = 'list',
        empCov = 'list',
        empCor = 'list',
        distMatrix = 'list',
        hTree = 'list'
        ),

    methods = list(

        initialize = function(code, constants=list(), data=list(), inits=list(), control=list()) {
            library(lattice)
            library(coda)
            library(nimble)
            abModel <<- autoBlockModel(code=code, constants=constants, data=data, inits=inits)
            defaultsList <- autoBlockParamDefaults()
            for(i in seq_along(defaultsList)) if(is.null(control[[names(defaultsList)[i]]])) control[[names(defaultsList)[i]]] <- defaultsList[[i]]
            for(i in seq_along(control)) eval(substitute(verbose <<- VALUE, list(verbose=as.name(names(control)[i]), VALUE=control[[i]])))
            it <<- 0
        },

        run = function(runList) {
            if(!is.list(runList)) stop('runList argument should be a list')
            if(is.null(names(runList))) names(runList) <- rep('', length(runList))
            
            for(i in seq_along(runList)) {
                runListElement <- runList[[i]]
                runListName <- names(runList)[i]
                if(is.character(runListElement)) {
                    type <- runListElement
                } else if(is.list(runListElement)) {
                    type <- 'blocks'
                } else if(class(runListElement) == '{') {
                    type <- 'spec'
                } else stop('don\'t understand element in run list')
                Rmodel <- abModel$newModel()
                switch(type,
                       
                       none =    { specList <- list(createSpecFromGroups(Rmodel, abModel$nodeGroupScalars))
                                   runSpecListAndSaveBest(Rmodel, specList, 'none') },

                       all =     { specList <- list(createSpecFromGroups(Rmodel, abModel$nodeGroupAllBlocked))
                                   runSpecListAndSaveBest(Rmodel, specList, 'all') },

                       default = { specList <- list(createSpecFromGroups(Rmodel, abModel$nodeGroupStochNodes, conjOveride=TRUE))
                                   runSpecListAndSaveBest(Rmodel, specList, 'default') },

                       blocks =  { specList <- list(createSpecFromGroups(Rmodel, abModel$createGroups(runListElement)))
                                   name <- if(runListName == '') 'customBlocks' else runListName
                                   runSpecListAndSaveBest(Rmodel, specList, name) },

                       spec =    { specList <- list(eval(runListElement, envir=environment()))
                                   name <- if(runListName == '') 'customSpec' else runListName
                                   runSpecListAndSaveBest(Rmodel, specList, name) },

                       auto =    { autoIt <- 0
                                   while((autoIt < 2) || ((!groupingsEquiv(grouping[[it]], grouping[[it-1]])) && (min(essPT[[it]]) > min(essPT[[it-1]])))) {
                                       Rmodel <- abModel$newModel()
                                       candidateGroupsList <- if(autoIt==0) list(abModel$nodeGroupScalars)  else determineCandidateGroupsFromCurrentSample()
                                       specList <- lapply(candidateGroupsList, function(groups) createSpecFromGroups(Rmodel, groups))
                                       runSpecListAndSaveBest(Rmodel, specList, paste0('auto',autoIt), auto=TRUE)
                                       autoIt <- autoIt + 1
                                   }
                               },
                       stop('don\'t understand element in run list'))
            }

            names(candidateGroups) <<- naming
            names(grouping) <<- naming
            names(groupSizes) <<- naming
            names(groupIDs) <<- naming
            names(samplers) <<- naming
            names(Cmcmcs) <<- naming
            names(timing) <<- naming
            if(saveSamples) names(samples) <<- naming
            names(means) <<- naming
            names(sds) <<- naming
            names(ess) <<- naming
            names(essPT) <<- naming
        },

        determineCandidateGroupsFromCurrentSample = function() {
            cutreeList <- lapply(cutree_heights, function(height) cutree(hTree[[it]], h = height))
            names(cutreeList) <- paste0('cut', cutree_heights)
            uniqueCutreeList <- unique(cutreeList)
            for(i in seq_along(uniqueCutreeList)) { for(j in seq_along(cutreeList)) { if(all(uniqueCutreeList[[i]]==cutreeList[[j]])) { names(uniqueCutreeList)[i] <- names(cutreeList)[j]; break } } }
            candidateGroupsList <- lapply(uniqueCutreeList, function(ct) determineGroupsFromCutree(ct))
            return(candidateGroupsList)
        },
        
        determineGroupsFromCutree = function(ct) {
            return(lapply(unique(ct), function(x) names(ct)[ct==x]))
        },
        
        runSpecListAndSaveBest = function(Rmodel, specList, name, auto=FALSE) {
            RmcmcList <- timingList <- samplesList <- meansList <- sdsList <- essList <- essPTList <- essPTminList <- list()
            for(i in seq_along(specList)) {
                checkOverMCMCspec(specList[[i]])
                specList[[i]]$addMonitors(abModel$monitorsVector, print=FALSE)
                RmcmcList[[i]] <- buildMCMC(specList[[i]])
            }
            toCompileList <- c(list(Rmodel), RmcmcList)
            compiledList <- compileNimble(toCompileList)
            Cmodel <- compiledList[[1]]
            CmcmcList <- compiledList[-1]
            for(i in seq_along(CmcmcList)) {
                Cmodel$setInits(abModel$inits)
                if(setSeed0) set.seed(0)
                timingList[[i]] <- as.numeric(system.time(CmcmcList[[i]](niter))[3])
                samplesList[[i]] <- as.matrix(nfVar(CmcmcList[[i]], 'mvSamples'))
                meansList[[i]] <- apply(samplesList[[i]], 2, mean)
                sdsList[[i]]   <- apply(samplesList[[i]], 2, sd)
                essList[[i]]   <- apply(samplesList[[i]], 2, effectiveSize)
                
                if(!saveSamples) samplesList[[i]] <- NA
                
                essPTList[[i]] <- essList[[i]] / timingList[[i]]
                essPTminList[[i]] <- sort(essPTList[[i]])[1]
            }
            bestInd <- as.numeric(which(unlist(essPTminList) == max(unlist(essPTminList))))
            if(!is.null(names(specList))) name <- paste0(name, '-', names(specList)[bestInd])
            
            it <<- it + 1
            naming[[it]] <<- name
            candidateGroups[[it]] <<- lapply(specList, function(spec) determineGroupsFromSpec(spec))
            grouping[[it]] <<- candidateGroups[[it]][[bestInd]]
            groupSizes[[it]] <<- determineNodeGroupSizesFromGroups(grouping[[it]])
            groupIDs[[it]] <<- determineNodeGroupIDsFromGroups(grouping[[it]])
            samplers[[it]] <<- determineSamplersFromGroupsAndSpec(grouping[[it]], specList[[bestInd]])
            Cmcmcs[[it]] <<- CmcmcList[[bestInd]]
            timing[[it]] <<- timingList[[bestInd]]
            samples[[it]] <<- samplesList[[bestInd]]
            means[[it]] <<- meansList[[bestInd]]
            sds[[it]] <<- sdsList[[bestInd]]
            ess[[it]] <<- essList[[bestInd]]
            essPT[[it]] <<- sort(essPTList[[bestInd]])
            
            if(auto) {
                samp <- as.matrix(nfVar(CmcmcList[[bestInd]], 'mvSamples'))
                burnedSamples[[it]] <<- samp[(floor(niter/2)+1):niter, ]
                empCov[[it]] <<- cov(burnedSamples[[it]])
                empCor[[it]] <<- cov2cor(empCov[[it]])
                distMatrix[[it]] <<- as.dist(1 - abs(empCor[[it]]))
                hTree[[it]] <<- hclust(distMatrix[[it]])
            }
            
            if(!saveSamples) burnedSamples[[it]] <<- NA
            
            if(verbose) printCurrent(name, specList[[bestInd]])
            if(makePlots && auto) makeCurrentPlots(name)
        },

        determineGroupsFromSpec = function(spec) {
            groups <- list()
            for(ss in spec$samplerSpecs) {
                if(ss$type == 'RW_block') {
                    nodes <- ss$control$targetNodes
                } else if(ss$type == 'crossLevel') {
                    topNodes <- ss$control$topNodes
                    lowNodes <- spec$model$getDependencies(topNodes, self=FALSE, stochOnly=TRUE, includeData=FALSE)
                    nodes <- c(topNodes, lowNodes)
                } else if(!is.null(ss$control$targetNode)) {
                    nodes <- ss$control$targetNode
                } else stop('don\'t understand sampler type')
                groups[[length(groups)+1]] <- spec$model$expandNodeNames(nodes, returnScalarComponents=TRUE)
            }
            return(groups)
        },

        determineNodeGroupSizesFromGroups = function(groups) {
            groupSizeVector <- numeric(0)
            for(gp in groups) for(node in gp) groupSizeVector[[node]] <- length(gp)
            return(groupSizeVector)
        },

        determineNodeGroupIDsFromGroups = function(groups) {
            groupIDvector <- numeric(0)
            for(i in seq_along(groups)) for(node in groups[[i]]) groupIDvector[[node]] <- i
            return(groupIDvector)
        },

        determineSamplersFromGroupsAndSpec = function(groups, spec) {
            samplerSpecs <- spec$samplerSpecs
            if(length(groups) != length(samplerSpecs)) stop('something wrong')
            samplerVector <- character(0)
            for(i in seq_along(groups)) for(node in groups[[i]]) samplerVector[[node]] <- samplerSpecs[[i]]$type
            return(samplerVector)
        },
        
        createSpecFromGroups = function(Rmodel, groups, conjOveride=FALSE) {
            spec <- MCMCspec(Rmodel, nodes=NULL, monitors=character(0))
            for(nodeGroup in groups) addSamplerToSpec(Rmodel, spec, nodeGroup, conjOveride)
            return(spec)
        },
        
        addSamplerToSpec = function(Rmodel, spec, nodeGroup, conjOveride) {
            if(length(nodeGroup) > 1) {
                spec$addSampler(type = 'RW_block', control = list(targetNodes=nodeGroup, adaptInterval=adaptIntervalBlock), print = FALSE); return()
            }
            if(!(nodeGroup %in% Rmodel$getNodeNames())) {
                spec$addSampler(type = 'RW', control = list(targetNode=nodeGroup), print = FALSE); return()
            }
            ## if(nodeGroup %in% Rmodel$getMaps('nodeNamesEnd')) {
            ##     cat(paste0('warning: using \'end\' sampler for node ', nodeGroup, ' may lead to results we don\'t want\n\n'))
            ##     spec$addSampler(type = 'end', control = list(targetNode=nodeGroup), print = FALSE); return()
            ## }
            conjugacyResult <- Rmodel$checkConjugacy(nodeGroup)
            if((!is.null(conjugacyResult)) && conjOveride) {
                spec$addSampler(type = conjugacyResult$samplerType, control = conjugacyResult$control, print = FALSE); return()
            }
            if(Rmodel$getNodeInfo()[[nodeGroup]]$isDiscrete()) {
                spec$addSampler(type = 'slice', control = list(targetNode=nodeGroup), print = FALSE); return()
            }
            if(length(Rmodel$expandNodeNames(nodeGroup, returnScalarComponents = TRUE)) > 1) {
                spec$addSampler(type = 'RW_block', control = list(targetNodes=nodeGroup, adaptInterval=adaptIntervalBlock), print = FALSE); return()
            }
            spec$addSampler(type = 'RW', control = list(targetNode=nodeGroup), print = FALSE); return()
        },

        checkOverMCMCspec = function(spec) {
            for(ss in spec$samplerSpecs) {
                if(ss$type == 'end') {
                    msg <- 'using \'end\' sampler may lead to results we don\'t want'
                    cat(paste0('\nWARNING: ', msg, '\n')); warning(msg)
                }
                if(grepl('^conjugate_', ss$type) && nimble:::nimbleOptions$verifyConjugatePosterior) {
                    msg <- 'conjugate sampler running slow due to checking the posterior'
                    cat(paste0('\nWARNING: ', msg, '\n')); warning(msg)
                }
            }
        },

        printCurrent = function(name, spec) {
            cat(paste0('\n################################\nbegin iteration ', it, ': ', name, '\n################################\n'))
            if(length(candidateGroups[[it]]) > 1) { cat('\ncandidate groups:\n'); cg<-candidateGroups[[it]]; for(i in seq_along(cg)) { cat(paste0('\n',names(cg)[i],':\n')); printGrouping(cg[[i]]) } }
            cat('\ngroups:\n'); printGrouping(grouping[[it]])
            cat('\nsamplers:\n'); spec$getSamplers()
            cat(paste0('\nMCMC runtime: ', round(timing[[it]], 1), ' seconds\n'))
            cat('\nESS:\n'); print(round(ess[[it]], 0))
            cat('\nESS/time:\n'); print(round(essPT[[it]], 1))
            cat(paste0('\n################################\nend iteration ', it, ': ', name, '\n################################\n\n'))
        },

        makeCurrentPlots = function(name) {
            dev.new()
            if(inherits(try(plot(as.dendrogram(hTree[[it]]), ylim=c(0,1), main=name), silent=TRUE), 'try-error')) dev.off()
        },

        printGrouping = function(g) {
            for(i in seq_along(g)) cat(paste0('[', i, '] ', paste0(g[[i]], collapse=', '), '\n'))
        },

        groupingsEquiv = function(grouping1, grouping2) {
            grouping1 <- lapply(grouping1, sort)
            grouping2 <- lapply(grouping2, sort)
            while(length(grouping1) > 0) {
                grp1 <- grouping1[[1]]
                found <- FALSE
                for(i in seq_along(grouping2)) {
                    grp2 <- grouping2[[i]]
                    if(identical(grp1, grp2)) {
                        found <- TRUE
                        grouping1[1] <- grouping2[i] <- NULL
                        break
                    }
                }
                if(!found) return(FALSE)
            }
            if(length(grouping2) == 0) return(TRUE) else return(FALSE)
        }
        )
    )


createDFfromABlist <- function(lst) {
    df <- data.frame(model=character(), blocking=character(), timing=numeric(), node=character(), groupSize = numeric(), groupID = numeric(), sampler = character(), mean=numeric(), sd=numeric(), ess=numeric(), essPT=numeric(), stringsAsFactors=FALSE)
    for(iAB in seq_along(lst)) {
        ab <- lst[[iAB]]
        abName <- names(lst)[iAB]
        for(iBlock in seq_along(ab$naming)) {
            blocking <- ab$naming[[iBlock]]
            timing <- ab$timing[[iBlock]]
            ess <- ab$ess[[iBlock]]
            nodes <- names(ess)
            means <- ab$means[[iBlock]][nodes]            ## sort
            sds <- ab$sds[[iBlock]][nodes]                ##
            essPT <- ab$essPT[[iBlock]][nodes]            ##
            groupSizes <- ab$groupSizes[[iBlock]][nodes]  ##
            groupIDs <- ab$groupIDs[[iBlock]][nodes]      ##
            samplers <- ab$samplers[[iBlock]][nodes]      ##
            newIndDF <- (1:length(nodes)) + dim(df)[1]
            df[newIndDF,] <- NA
            df[newIndDF,]$model <- abName
            df[newIndDF,]$blocking <- blocking
            df[newIndDF,]$timing <- timing
            df[newIndDF,]$node <- nodes
            df[newIndDF,]$groupSize <- groupSizes
            df[newIndDF,]$groupID <- groupIDs
            df[newIndDF,]$sampler <- samplers
            df[newIndDF,]$mean <- means
            df[newIndDF,]$sd <- sds
            df[newIndDF,]$ess <- ess
            df[newIndDF,]$essPT <- essPT
        }
    }
    return(df)
}



plotABS <- function(df, xlimToMin=FALSE, together) {
    models <- unique(df$model)
    nModels <- length(models)
    if(missing(together)) together <- if(nModels <= 5) TRUE else FALSE
    nVertPlots <- if(together) nModels*2 else nModels
    xVarNames <- c('ess', 'essPT')
    parCmd <- quote(par(mfrow=c(nVertPlots,1),mar=c(1,0,1,0),tcl=-.1,mgp=c(3,0,0),cex.axis=.7))
    if(together) { quartz(); eval(parCmd) }
    for(xVarName in xVarNames) {
        if(!together) { quartz(); eval(parCmd) }
        maxMinXVar<-0; for(mod in models) {dfMod<-df[df$model==mod,]; blks<-unique(dfMod$blocking); for(blk in blks) {maxMinXVar<-max(maxMinXVar,min(dfMod[dfMod$blocking==blk,xVarName]))}}
        maxXVar <- if(xlimToMin) maxMinXVar else max(df[, xVarName])
        xlim <- c(maxXVar*-0.05, maxXVar)
        maxTiming <- max(df[, 'timing'])
        for(mod in models) {
            dfMod <- df[df$model==mod,]
            blockings <- unique(dfMod$blocking)
            nBlockings <- length(blockings)
##            bestBlk<-''; bestEssPT<-0; for(blk in blockings) { if(min(dfMod[dfMod$blocking==blk,'essPT'])>bestEssPT && ((blk=='all')||(blk=='default')||grepl('^auto',blk))) {bestEssPT<-min(dfMod[dfMod$blocking==blk,'essPT']); bestBlk<-blk} }
            bestBlk<-''; bestEssPT<-0; for(blk in blockings) { if(min(dfMod[dfMod$blocking==blk,'essPT'])>bestEssPT) {bestEssPT<-min(dfMod[dfMod$blocking==blk,'essPT']); bestBlk<-blk} }
            plot(-100,-100,xlim=xlim,ylim=c(0,nBlockings+1),xlab='',ylab='',main=paste0(xVarName, ' for model ', mod))
            for(iBlocking in 1:nBlockings) {
                blocking <- blockings[iBlocking]
                dfModBlock <- dfMod[dfMod$blocking==blocking,]
                xVarValues <- dfModBlock[,xVarName]
                groupSizes <- dfModBlock[,'groupSize']
                timing <- dfModBlock[,'timing'][1]   # only first element
                timingOnXaxis <- timing/maxTiming * xlim[2]
                yCoord <- nBlockings+1-iBlocking
                lines(x=c(0,timingOnXaxis), y=rep(yCoord,2), lty=1, lwd=2, col='lightgrey')
                col <- if(blocking == bestBlk) 'green' else 'black'
                text(x=xVarValues, y=yCoord, labels=groupSizes, cex=0.7, col=col)
                col <- if(blocking == bestBlk) 'green' else 'blue'
                text(x=xlim[1], y=yCoord, labels=blocking, col=col)
                if(timing==maxTiming) text(xlim[2], yCoord+1, paste0('t = ',round(timing,1)))
            }
        }
    }
}


printMinTimeABS <- function(df) {
    namesToRemove <- c('groupID', 'sampler', 'mean', 'sd')
    for(name in namesToRemove) { ind <- which(names(df)==name); df <- df[, -ind] }
    models <- unique(df$model)
    cat('\n')
    for(mod in models) {
        dfMod <- df[df$model == mod, ]
        blockings <- unique(dfMod$blocking)
        dfOut <- dfMod[numeric(0), ]
        for(blk in blockings) {
            dfModBlk <- dfMod[dfMod$blocking == blk, ]
            ind <- which(dfModBlk$essPT == min(dfModBlk$essPT))[1]
            dfOut[dim(dfOut)[1] + 1, ] <- dfModBlk[ind, ]
        }
        dfOut <- dfOut[sort(dfOut$essPT,index.return=TRUE)$ix, ]
        dimnames(dfOut)[[1]] <- 1:(dim(dfOut)[1])
        dfOut$timing <- round(dfOut$timing, 1)
        dfOut$ess    <- round(dfOut$ess, 0)
        dfOut$essPT  <- round(dfOut$essPT, 1)
        print(dfOut)
        cat('\n')
    }
}


codeToText <- function(code) {
    a <- deparse(code, width.cutoff=500L)
    a <- a[c(-1, -length(a))]
    a <- sub('^    ', '', a)
    a[length(a) + 1] <- ''
    a[length(a) + 1] <- ''
    a <- paste0(a, collapse='\n')
    return(a)
}





#######################################
########         MODELS        ########
#######################################

library(nimble)
if(!exists('control')) control <- list()

################
### tester
################

code_tester <- modelCode({
    x[1:3] ~ dmnorm(mu[1:3], Q[1:3,1:3])
    y[1:3] ~ dmnorm(x[1:3], Q[1:3, 1:3])
    z1 ~ dnorm(y[1], 100)
})

constants_tester <- list(mu=rep(0,3), Q=diag(3))
data_tester <- list(y = c(1,2,3))
inits_tester <- list(x=rep(0,3), z1=0)

##abtester <- autoBlock(code=code_tester, constants=constants_tester, data=data_tester, inits=inits_tester, control=control)

################
### litters
################

G <- 2
N <- 16
n <- array(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 11, 8, 10, 13, 10, 12, 9, 10, 9, 10, 5, 9, 9, 13, 7, 5, 10, 7, 6, 10, 10, 10, 7), dim = c(2, 16))
r <- array(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 10, 8, 9, 12, 9, 11, 8, 9, 8, 9, 4, 8, 7, 11, 4, 4, 5, 5, 3, 7, 3, 7, 0), dim = c(2, 16))
p <- array(0.5, dim = c(2, 16))
a <- c(1, 1)
b <- c(1, 1)

constants_litters <- list(G=G, N=N, n=n)
data_litters      <- list(r=r)
inits_litters     <- list(a=a, b=b, p=p)

code_litters <- modelCode({
#    a[1] ~ dunif(0, 80000)
#    b[1] ~ dunif(0, 10000)
    a[1] ~ dgamma(1, 0.001)   # works well
    b[1] ~ dgamma(1, 0.001)   # works well
    a[2] ~ dunif(0, 100)   # works well
    b[2] ~ dunif(0, 50)    # works well
     for (i in 1:G) {
#         a[i] ~ dgamma(1, 0.001)
#         b[i] ~ dgamma(1, 0.001)
         for (j in 1:N) {
             r[i,j] ~ dbin(p[i,j], n[i,j])
             p[i,j] ~ dbeta(a[i], b[i])
         }
     }
})

##ablitters <- autoBlock(code=code_litters, constants=constants_litters, data=data_litters, inits=inits_litters, control=control)

rm(list = c('G','N','n','r','p','a','b'))

################
### SSMmub
################

## better parameterization: mean and autocorrelation
code_SSMmub<- modelCode({
    mu ~ dnorm(0, sd = 1000)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(mu, sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    a <- 1-(b/mu)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})

t <- 30
constants_SSMmub <- list(t = t)
Rmodel <- nimbleModel(code_SSMmub, constants = constants_SSMmub)
Rmodel$mu <- 10/(1-.5)
Rmodel$b <- 10
Rmodel$sigPN <- .1
Rmodel$sigOE <- .1
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('mu','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data_SSMmub <- list(y = Rmodel$y)
inits_SSMmub <- list(mu = Rmodel$mu, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)

##abSSMmub <- autoBlock(code=code_SSMmub, constants=constants_SSMmub, data=data_SSMmub, inits=inits_SSMmub, control=control)

rm(list = c('t','Rmodel'))

################
### SSMab
################

## parameterization in terms of slope (autocorelation) and intercept,
## which are highly correlated in mixing
code_SSMab <- modelCode({
    a ~ dunif(-0.9999, 0.9999)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
    x[1] ~ dnorm(b/(1-a), sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})

t <- 30
constants_SSMab <- list(t = t)
Rmodel <- nimbleModel(code_SSMab, constants = constants_SSMab)
Rmodel$a <- .5
Rmodel$b <- 10
Rmodel$sigPN <- .1
Rmodel$sigOE <- .1
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('a','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data_SSMab <- list(y = Rmodel$y)
inits_SSMab <- list(a = Rmodel$a, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)

##abSSMab <- autoBlock(code=code_SSMab, constants=constants_SSMab, data=data_SSMab, inits=inits_SSMab, control=control)

rm(list = c('t','Rmodel'))


################
### scallops
################

library(Imap)
myscallops <- read.table("http://www.biostat.umn.edu/~brad/data/myscallops.txt", header = TRUE)
#####myscallops <- myscallops[1:25,]
N <- dim(myscallops)[1]
catch <- myscallops$tcatch
lat <- myscallops$lat
long <- myscallops$long
dist <- array(NA, c(N,N))
for(i in 1:N) for(j in 1:N) dist[i,j] <- gdist(long[i], lat[i], long[j], lat[j])

code_scallops<- modelCode({
    mu ~ dunif(0, 100)
    sigma ~ dunif(0, 100)
    rho ~ dunif(20, 100)
    muVec[1:N] <- mu * onesVector[1:N]
    Cov[1:N,1:N] <- sigma^2 * exp(-dist[1:N,1:N] / rho)
    g[1:N] ~ dmnorm(muVec[1:N], cov = Cov[1:N,1:N])
    for(i in 1:N) {
        y[i] ~ dpois(exp(g[i]))
    }
})

constants_scallops <- list(N=N, onesVector=rep(1,N), dist=dist)
data_scallops <- list(y=catch)
inits_scallops <- list(mu=5, sigma=1, rho=30, g=rep(0,N))

##abscallops <- autoBlock(code=code_scallops, constants=constants_scallops, data=data_scallops, inits=inits_scallops)

rm(list = c('myscallops', 'N', 'catch', 'lat', 'long', 'dist', 'i', 'j'))



################
### dipper
################
## load('dipperData.RData')
## ## optionally, truncate data
## ind <- 1:3;   nind<-length(ind);   first<-first[ind];   y<-y[ind,,drop=FALSE];   x_init<-x_init[ind,,drop=FALSE]
## code_dipper <- modelCode({
##     phi ~ dunif(0, 1)
##     p ~ dunif(0, 1)
##     for (i in 1:nind) {
##         x[i, first[i]] <- 1
##         for (t in (first[i] + 1):k) {
##             mu.x[i, t] <- phi * x[i, t-1]
##             mu.y[i, t] <- p * x[i, t]
##             x[i, t] ~ dbin(mu.x[i, t], 1)
##             y[i, t] ~ dbin(mu.y[i, t], 1)
##         }
##     }
## })
## 
## constants_dipper <- list(k=k, nind=nind, first=first)
## data_dipper      <- list(y=y)
## inits_dipper     <- list(phi=0.6, p=0.9, x=x_init)
## abdipper <- autoBlock(code=code_dipper, constants=constants_dipper, data=data_dipper, inits=inits_dipper, control=control)
## if(exists('ind')) rm(list = 'ind')
## rm(list = c('first', 'k', 'nind', 'x_init', 'y'))




## cutree_custom <- function(ht, maxHeight, maxGroupSize, maxHeightRelativeFromBase) {
##     labels <- ht$labels;     height <- ht$height;     merge <- ht$merge
##     nNodes <- length(labels)
##     nMerges <- dim(merge)[1]
##     if(nMerges+1 != nNodes) stop('something fishy')
##     df <- data.frame(baseHeight=numeric(nNodes), recentHeight = numeric(nNodes), num=numeric(nNodes))
##     for(i in 1:nMerges) { for(j in 1:2) {
##         if(merge[i,j] < 0) df$baseHeight[abs(merge[i,j])] <- df$recentHeight[abs(merge[i,j])] <- height[i] } }
##     df$num <- rep(1, nNodes)
##     nodes <- lapply(1:nNodes, function(x) x)
##     for(i in 1:nMerges) {
##         ind1 <- if(merge[i,1] < 0) abs(merge[i,1]) else (merge[i,1]+nNodes)
##         ind2 <- if(merge[i,2] < 0) abs(merge[i,2]) else (merge[i,2]+nNodes)
##         indout <- i + nNodes
##         df[indout,] <- NA
##         if((df$num[ind1]==0) || (df$num[ind2]==0) || is.na(df$num[ind1]) || is.na(df$num[ind2])) stop('something wrong')
##         canMerge <- TRUE
##         if(height[i] > maxHeight) canMerge <- FALSE
##         if(df$num[ind1] + df$num[ind2] > maxGroupSize) canMerge <- FALSE
##         bh <- min(df$baseHeight[c(ind1,ind2)])
##         if(height[i] > bh + maxHeightRelativeFromBase*(1-bh)) canMerge <- FALSE
##         if(canMerge) {
##             ## merge
##             df$baseHeight[indout] <- min(df$baseHeight[c(ind1,ind2)])
##             df$recentHeight[indout] <- height[i]
##             df$num[indout] <- df$num[ind1] + df$num[ind2]
##             nodes[[indout]] <- c(nodes[[ind2]], nodes[[ind1]])
##             if(df$num[indout] != length(nodes[[indout]])) stop('something wrong')
##             df[ind1,] <- NA;     df[ind2,] <- NA
##             df$num[c(ind1,ind2)] <- c(0,0)
##             nodes[[ind1]] <- NA;     nodes[[ind2]] <- NA
##         } else {
##             ## don't merge; lower branch becomes a group, higher branch propagates up
##             indhigher <- if(df$recentHeight[ind1] > df$recentHeight[ind2]) ind1 else ind2
##             df[indout,] <- df[indhigher,]
##             df[indhigher,] <- NA
##             df$num[indhigher] <- 0
##             nodes[[indout]] <- nodes[[indhigher]]
##             nodes[[indhigher]] <- NA
##         }
##     }
##     groupInd <- which(df$num > 0)
##     groupNodeNumbers <- nodes[groupInd]
##     groups <- lapply(groupNodeNumbers, function(gnn) labels[gnn])
##     return(groups)
## }



