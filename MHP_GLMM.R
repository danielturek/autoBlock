



##
## create constants from 'mhp' data frame
##
if(FALSE) {
    
    setwd('~/GitHub/autoBlock')
    rm(list=ls())
    load('mhp.RData')
    
    N <- dim(mhp)[1]
    count <- as.numeric(mhp$Count)
    subject <- as.numeric(mhp$Subject)
    event <- as.numeric(factor(mhp$Event, levels=1:2))
    period <- as.numeric(factor(mhp$Period, levels=1:4))
    eventperiod <- period
    eventperiod[event==2] <- 4
    Nsubject <- max(subject)
    subjectevent <- 2*subject-event+1
    Nsubjectevent <- max(subjectevent)
    subjectperiod <- 4*subject-period+1
    Nsubjectperiod <- max(subjectperiod)
    
    save(mhp, N, count, subject, event, period, eventperiod, Nsubject, subjectevent, Nsubjectevent, subjectperiod, Nsubjectperiod, file = 'mhp2.RData')
    
}

##
## create code, constants, data, and inits
##
if(FALSE) {
    
    setwd('~/GitHub/autoBlock')
    rm(list=ls())
    load('mhp2.RData')
    
    code_mhp <- quote({
        ## fixed effects
        mu ~ dnorm(0, 0.001)
        be[1] ~ dnorm(0, 0.001)
        be[2] <- 0
        for(i in 1:3) {
            bp[i] ~ dnorm(0, 0.001)
        }
        bp[4] <- 0
        for(i in 1:3) {
            bep[i] ~ dnorm(0, 0.001)
        }
        bep[4] <- 0
        ## random effects
        sds ~ dunif(0, 10)
        taus <- pow(sds, -2)
        for(i in 1:Nsubject) {
            us[i] ~ dnorm(0, taus)
        }
        sdse ~ dunif(0, 10)
        tause <- pow(sdse, -2)
        for(i in 1:Nsubjectevent) {
            use[i] ~ dnorm(0, tause)
        }
        sdsp ~ dunif(0, 10)
        tausp <- pow(sdsp, -2)
        for(i in 1:Nsubjectperiod) {
            usp[i] ~ dnorm(0, tausp)
        }
        ## likelihood
        for(i in 1:N) {
            eta[i] <- mu + be[event[i]] + bp[period[i]] + bep[eventperiod[i]] + us[subject[i]] + use[subjectevent[i]] + usp[subjectperiod[i]]
            lambda[i] <- exp(eta[i])
            count[i] ~ dpois(lambda[i])
        }
    })

    constants_mhp <- list(N=N, Nsubject=Nsubject, Nsubjectevent=Nsubjectevent, Nsubjectperiod=Nsubjectperiod, event=event, eventperiod=eventperiod, period=period, subject=subject, subjectevent=subjectevent, subjectperiod=subjectperiod)

    data_mhp <- list(count=count)

    inits_mhp <- list(mu=0, be=c(0,NA), bp=c(0,0,0,NA), bep=c(0,0,0,NA), sds=0.5, sdse=0.6, sdsp=0.6, us=rep(0,Nsubject), use=rep(0,Nsubjectevent), usp=rep(0,Nsubjectperiod))

    save(list = c(ls(), 'code_mhp', 'constants_mhp', 'data_mhp', 'inits_mhp'), file = 'mhp3.RData')
    
}


##
## run JAGS
##
if(FALSE) {
    
    setwd('~/GitHub/autoBlock')
    rm(list=ls())
    library(rjags)
    load('mhp3.RData')
    code <- code_mhp
    constants <- constants_mhp
    data <- data_mhp
    inits <- inits_mhp

    ## MCMC control
    monitorVars <- c('mu', 'be', 'bp', 'bep', 'sds', 'sdse', 'sdsp')
    niter <- 20000
    thin <- 1

    ## additional processing -- don't touch
    constsAndData <- c(constants, data)
    modelfile <- file.path(tempdir(), 'model.txt')
    writeLines(paste0('model\n', paste0(deparse(code), collapse='\n')), con=modelfile)

    jags_mod <- jags.model(file=modelfile, data=constsAndData, inits=inits, n.chains=1, quiet=FALSE)
    jags_out <- coda.samples(model=jags_mod, variable.names=monitorVars, n.iter=niter, thin=thin)
    
    dimnames(jags_out[[1]])
    means <- apply(jags_out[[1]][,], 2, mean)
    sds <- apply(jags_out[[1]][,], 2, sd)
    cbind(means, sds)  # compare to: http://glmm.wikidot.com/minnesota-health-plan
    
}


##
## run NIMBLE default MCMC
##
if(FALSE) {
    
    setwd('~/GitHub/autoBlock')
    rm(list=ls())
    library(nimble)
    load('mhp3.RData')

    Rmodel <- nimbleModel(code=code_mhp, constants=constants_mhp, data=data_mhp, inits=inits_mhp)

    spec <- configureMCMC(Rmodel)
    spec$getSamplers()
    spec$getMonitors()

    Rmcmc <- buildMCMC(spec)

    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

    Cmcmc$run(20000)

    samples <- as.matrix(Cmcmc$mvSamples)
    means <- apply(samples, 2, mean)
    sds <- apply(samples, 2, sd)
    cbind(means, sds)  # compare to: http://glmm.wikidot.com/minnesota-health-plan

}






