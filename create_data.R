

################
### litters (aka 'random effects')
################

rm(list = ls())
modelName <- 'litters'
G <- 2
N <- 16
n <- array(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 11, 8, 10, 13, 10, 12, 9, 10, 9, 10, 5, 9, 9, 13, 7, 5, 10, 7, 6, 10, 10, 10, 7), dim = c(2, 16))
r <- array(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 10, 8, 9, 12, 9, 11, 8, 9, 8, 9, 4, 8, 7, 11, 4, 4, 5, 5, 3, 7, 3, 7, 0), dim = c(2, 16))
p <- array(0.5, dim = c(2, 16))
a <- c(1, 1)
b <- c(1, 1)
constants <- list(G=G, N=N, n=n)
data      <- list(r=r)
inits     <- list(a=a, b=b, p=p)
code <- quote({
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
runList <- list('all',
                blockAB = list(c('a[1]','b[1]'), c('a[2]','b[2]')),
                crossLevel = quote({
                    spec <- configureMCMC(oldSpec = abModel$initialMCMCspec)  ## new version
                    spec$setSamplers()  ## new version -- removes all the samplers from initalMCMCspec
                    spec$addSampler('crossLevel', c('a[1]', 'b[1]'), print=FALSE)
                    spec$addSampler('crossLevel', c('a[2]', 'b[2]'), print=FALSE)
                    spec
                }),
                'default',
                'auto')
modelfileName <- paste0('~/GitHub/legacy/autoBlock/data/model_', modelName, '.RData')
save(code, constants, data, inits, runList, file = modelfileName)



################
## ice (aka 'auto-regressive')
################

rm(list = ls())
modelName <- 'ice'
code <- quote({
    for (i in 1:I) {
        eta[i] <- log(pyr[i]) + alpha[age[i]] + beta[year[i]]
        lambda[i] <- exp(eta[i])
        cases[i] ~ dpois(lambda[i])
    }
    betamean[1] <- 0
    betaprec[1] <- tau * 0.00001
    betamean[2] <- 0.0
    betaprec[2] <- tau * 0.00001
    for (k in 3:K) {
        betamean[k] <- 2*beta[k-1] - beta[k-2]
        betaprec[k] <- tau
    } 
    for (k in 1:K) {
        beta[k] ~ dnorm(betamean[k], betaprec[k])
        logRR[k] <- beta[k] - beta[5]
    }
    alpha[1] <- 0
    for (j in 2:Nage) {
        alpha[j] ~ dnorm(0, 0.00001)
    }
    sigma ~ dunif(0,1)
    tau <- 1/(sigma*sigma)
})
constants <- list(
    I = 77,
    K = 11,
    Nage = 13,
    pyr = c(41380, 43650, 49810, 58105, 57105, 76380, 39615, 
        42205, 48315, 56785, 55965, 33955, 29150, 38460, 40810, 47490, 
        55720, 55145, 27950, 37375, 39935, 46895, 54980, 27810, 25055, 
        27040, 36400, 39355, 46280, 54350, 24040, 26290, 35480, 38725, 
        45595, 25710, 22890, 23095, 25410, 34420, 37725, 44740, 21415, 
        21870, 24240, 33175, 36345, 21320, 17450, 19765, 20255, 22760, 
        31695, 34705, 15350, 17720, 18280, 20850, 29600, 15635, 9965, 
        12850, 15015, 15725, 18345, 26400, 8175, 11020, 13095, 14050, 
        16480, 10885, 7425, 10810, 12260, 14780, 13600),
    age = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 
        3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 
        7, 7, 7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 
        10, 10, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 
        13, 13, 13, 13),
    year = c(6, 7, 8, 9, 10, 11, 6, 7, 8, 9, 10, 11, 5, 6, 7, 
        8, 9, 10, 5, 6, 7, 8, 9, 10, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 
        9, 3, 4, 5, 6, 7, 8, 3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 2, 3, 
        4, 5, 6, 7, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5)
)
data <- list(
    cases = c(2, 0, 1, 1, 1, 2, 0, 2, 1, 1, 5, 5, 1, 1, 3, 7, 
        12, 10, 6, 11, 9, 14, 20, 14, 7, 14, 22, 25, 29, 37, 21, 11, 
        29, 33, 57, 24, 15, 8, 22, 27, 38, 52, 10, 15, 22, 26, 47, 31, 
        8, 11, 17, 23, 31, 38, 8, 10, 24, 30, 53, 26, 5, 3, 10, 18, 22, 
        30, 1, 7, 11, 26, 32, 17, 5, 8, 17, 32, 31)
)
inits <- list(
    tau = 1,
    alpha = c(NA,0,0,0,0,0,0,0,0,0,0,0,0),
    beta = c(0.05,0.1,0,0,0,0,0,0,0,0,0)
)
runList <- list('all',
                'default',
                'auto')
modelfileName <- paste0('~/GitHub/legacy/autoBlock/data/model_', modelName, '.RData')
save(code, constants, data, inits, runList, file = modelfileName)



################
### SSMindependent
################

rm(list = ls())
modelName <- 'SSMindependent'
library(nimble)
## better parameterization: mean and autocorrelation
code <- quote({
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
t <- 100
constants <- list(t = t)
Rmodel <- nimbleModel(code = code, constants = constants)
## Rmodel$mu <- 10/(1-.5)   ## original SSM
Rmodel$mu <- 1/(1-.95)   ## next SSM attempt (v2)
Rmodel$b <- 1
## Rmodel$sigPN <- .1  ## original SSM
Rmodel$sigPN <- .2  ## next SSM attempt (v2)
Rmodel$sigOE <- .05
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('mu','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data <- list(y = Rmodel$y)
inits <- list(mu = Rmodel$mu, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)
runList <- list('all',
                'default',
                'auto')
modelfileName <- paste0('~/GitHub/legacy/autoBlock/data/model_', modelName, '.RData')
save(code, constants, data, inits, runList, file = modelfileName)



################
### SSMcorrelated
################

rm(list = ls())
modelName <- 'SSMcorrelated'
library(nimble)
## parameterization in terms of slope (autocorelation) and intercept,
## which are highly correlated in mixing
code <- quote({
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
t <- 100
constants <- list(t = t)
Rmodel <- nimbleModel(code = code, constants = constants)
## Rmodel$a <- .5  ## original SSM
Rmodel$a <- .95  ## next SSM attempt (v2)
Rmodel$b <- 1
## Rmodel$sigPN <- .1  ## original SSM
Rmodel$sigPN <- .2  ## next SSM attempt (v2)
Rmodel$sigOE <- .05
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('a','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data <- list(y = Rmodel$y)
inits <- list(a = Rmodel$a, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)
runList <- list('all',
                blockAB = list(c('a', 'b')),
                'default',
                'auto')
modelfileName <- paste0('~/GitHub/legacy/autoBlock/data/model_', modelName, '.RData')
save(code, constants, data, inits, runList, file = modelfileName)



################
### spatial (scallop abundance)
################

rm(list = ls())
modelName <- 'spatial'
library(Imap)
myscallops <- read.table('http://www.biostat.umn.edu/~brad/data/myscallops.txt', header = TRUE)
#####myscallops <- myscallops[1:25,]
N <- dim(myscallops)[1]
catch <- myscallops$tcatch
lat <- myscallops$lat
long <- myscallops$long
dist <- array(NA, c(N,N))
for(i in 1:N) for(j in 1:N) dist[i,j] <- gdist(long[i], lat[i], long[j], lat[j])
code <- nimbleCode({
    mu ~ dunif(-100, 100)
    sigma ~ dunif(0, 100)
    rho ~ dunif(20, 100)
    muVec[1:N] <- mu * onesVector[1:N]
    Cov[1:N,1:N] <- sigma^2 * exp(-dist[1:N,1:N] / rho)
    g[1:N] ~ dmnorm(muVec[1:N], cov = Cov[1:N,1:N])
    for(i in 1:N) {
        y[i] ~ dpois(exp(g[i]))
    }
})
constants <- list(N=N, onesVector=rep(1,N), dist=dist)
data <- list(y=catch)
inits <- list(mu=0, sigma=5, rho=60, g=rep(0,N))
runList <- list('all',
                'default',
                'auto')
modelfileName <- paste0('~/GitHub/legacy/autoBlock/data/model_', modelName, '.RData')
save(code, constants, data, inits, runList, file = modelfileName)



################
### mhp (aka 'GLMM'; Minnesota Health Plan, http://glmm.wikidot.com/minnesota-health-plan)
################

rm(list = ls())
modelName <- 'mhp'
## create constants from 'mhp' data frame
load('~/GitHub/legacy/autoBlock/data/mhp.RData')
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
## create code, constants, data, and inits
code <- quote({
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
constants <- list(N=N, Nsubject=Nsubject, Nsubjectevent=Nsubjectevent, Nsubjectperiod=Nsubjectperiod, event=event, eventperiod=eventperiod, period=period, subject=subject, subjectevent=subjectevent, subjectperiod=subjectperiod)
data <- list(count=count)
inits <- list(mu=0, be=c(0,NA), bp=c(0,0,0,NA), bep=c(0,0,0,NA), sds=0.5, sdse=0.6, sdsp=0.6, us=rep(0,Nsubject), use=rep(0,Nsubjectevent), usp=rep(0,Nsubjectperiod))
runList <- list('all',
                'default',
                'auto')
modelfileName <- paste0('~/GitHub/legacy/autoBlock/data/model_', modelName, '.RData')
save(code, constants, data, inits, runList, file = modelfileName)



################
### redblue (big GLMM from Gelman book)
################

rm(list = ls())
modelName <- 'redblue'
library(foreign)
## create constants from file: '2000_labeled_processed_race.dta', and 'st.dat'
## individual level data
df <- read.dta('~/GitHub/legacy/autoBlock/data/2000_labeled_processed_race.dta')
###df <- df[seq(1, dim(df)[1], 200), ]     ########## TEMPOPRARY #############
income <- as.numeric(df$income) - 3   ## contains NAs
state <- df$state
y <- df$y
naIndex <- is.na(income) | is.na(y)    ## get indexes of NAs
income <- income[!naIndex]  ## remove these from income, state, y
state <- state[!naIndex]    ##
y <- y[!naIndex]            ##
N <- length(y)
Nstates <- max(state)
## state level data
st <- read.table('~/GitHub/legacy/autoBlock/data/st.dat')
st2000 <- st[st$year == 2000, ]
stateIncome <- st2000$z.st.inc
## create code, constants, data, and inits
constants <- list(
    N = N,
    income = income,
    state = state,
    Nstates = Nstates,
    stateIncome = stateIncome
)
data <- list(y = y)
inits <- list(
    gamma = matrix(rep(0, 4), nrow = 2),
    sigmaIntercept = 1,
    sigmaSlope = 1,
    rho = 0
)
code <- quote({
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
runList <- list('all',
                'default',
                'auto')
modelfileName <- paste0('~/GitHub/legacy/autoBlock/data/model_', modelName, '.RData')
save(code, constants, data, inits, runList, file = modelfileName)




################
### test
################

rm(list = ls())
modelName <- 'test'

code <- quote({
    a ~ dnorm(0, 1)
    b ~ dnorm(a, 1)
    c ~ dnorm(b, 1)
})

constants <- list()
data <- list(c = 1)
inits <- list(a = 0, b = 0)

runList <- list('all',
                'default',
                'auto',
                blockAB = list(c('a', 'b')))

modelfileName <- paste0('~/GitHub/legacy/autoBlock/data/model_', modelName, '.RData')
save(code, constants, data, inits, runList, file = modelfileName)


