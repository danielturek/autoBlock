


##
## create constants from file: '2000_labeled_processed_race.dta', and 'st.dat'
##
if(FALSE) {
    
    setwd('~/GitHub/autoBlock')
    library(foreign)
    rm(list=ls())
    ## individual level data
    df <- read.dta('2000_labeled_processed_race.dta')
    ###df <- df[seq(1, dim(df)[1], 200), ]     ######################################################## TEMPOPRARY ###############################################
    income <- as.numeric(df$income) - 3   ## contains NAs
    state <- df$state
    y <- df$y
    naIndex <- is.na(income) | is.na(y)    ## get indexes of NAs
    income <- income[!naIndex]  ## remove thewe from income, state, y
    state <- state[!naIndex]    ##
    y <- y[!naIndex]            ##
    N <- length(y)
    Nstates <- max(state)
    ## state level data
    st <- read.table('st.dat')
    st2000 <- st[st$year == 2000, ]
    stateIncome <- st2000$z.st.inc
    save(N, income, state, Nstates, stateIncome, y, file = 'redblue2.RData')
    
}

##
## create code, constants, data, and inits
## BUGS model, data, inits for Red State Blue State example from Gelman book
##
if(FALSE) {
    
    setwd('~/GitHub/autoBlock')
    rm(list=ls())
    load('redblue2.RData')
    constants_redblue <- list(
        N = N,
        income = income,
        state = state,
        Nstates = Nstates,
        stateIncome = stateIncome
    )
    data_redblue <- list(y = y)
    inits_redblue <- list(
        gamma = matrix(rep(0, 4), nrow = 2),
        sigmaIntercept = 1,
        sigmaSlope = 1,
        rho = 0
    )
    code_redblue <- quote({
        ## top-level parameters
        for(i in 1:2) {     for(j in 1:2) {     gamma[i, j] ~ dnorm(0, 0.0001)     }     }
        sigmaIntercept ~ dunif(0, 100)
        sigmaSlope ~ dunif(0, 100)
        rho ~ dunif(-1, 1)
        ## Sigma matrix
        Sigma[1, 1] <- sigmaIntercept^2
        Sigma[2, 2] <- sigmaSlope^2
        Sigma[1, 2] <- rho * sigmaIntercept * sigmaSlope
        Sigma[2, 1] <- rho * sigmaIntercept * sigmaSlope
        ## latent "states"
        for(i in 1:Nstates) {
            stateBetaMeans[i, 1] <- gamma[1, 1] + gamma[1, 2] * stateIncome[i]
            stateBetaMeans[i, 2] <- gamma[2, 1] + gamma[2, 2] * stateIncome[i]
            stateBetas[i, 1:2] ~ dmnorm(mean = stateBetaMeans[i, 1:2], cov = Sigma[1:2, 1:2])
        }
        ## likelihood
        for(i in 1:N) {
            eta[i] <- stateBetas[ state[i], 1 ] + stateBetas[ state[i], 2 ] * income[i]
            p[i] <- expit( eta[i] )
            y[i] ~ dbern( p[i] )
        }
    })
    code_redblueJAGS <- quote({
        ## top-level parameters
        for(i in 1:2) {     for(j in 1:2) {     gamma[i, j] ~ dnorm(0, 0.0001)     }     }
        sigmaIntercept ~ dunif(0, 100)
        sigmaSlope ~ dunif(0, 100)
        rho ~ dunif(-1, 1)
        ## Sigma matrix
        Prec[1, 1] <- 1 / (sigmaIntercept^2 * (1-rho^2))
        Prec[2, 2] <- 1 / (sigmaSlope^2     * (1-rho^2))
        Prec[1, 2] <- -1 * rho / (sigmaIntercept * sigmaSlope * (1-rho^2))
        Prec[2, 1] <- -1 * rho / (sigmaIntercept * sigmaSlope * (1-rho^2))
        ## latent "states"
        for(i in 1:Nstates) {
            stateBetaMeans[i, 1] <- gamma[1, 1] + gamma[1, 2] * stateIncome[i]
            stateBetaMeans[i, 2] <- gamma[2, 1] + gamma[2, 2] * stateIncome[i]
            stateBetas[i, 1:2] ~ dmnorm(stateBetaMeans[i, 1:2], Prec[1:2, 1:2])
        }
        ## likelihood
        for(i in 1:N) {
            eta[i] <- stateBetas[ state[i], 1 ] + stateBetas[ state[i], 2 ] * income[i]
            logit(p[i]) <- eta[i]
            y[i] ~ dbern( p[i] )
        }
    })
    save(list = c(ls(), 'code_redblue', 'code_redblueJAGS', 'constants_redblue', 'data_redblue', 'inits_redblue'), file = 'redblue3.RData')
    
}


##
## run NIMBLE default MCMC
##
if(FALSE) {
    
    setwd('~/GitHub/autoBlock')
    rm(list=ls())
    library(nimble)
    load('redblue3.RData')
    Rmodel <- nimbleModel(code=code_redblue, constants=constants_redblue, data=data_redblue, inits=inits_redblue)
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
    cbind(means, sds)

}













