library(nimble)

## parameterization in terms of slope (autocorelation) and intercept,
## which are highly correlated in mixing
## code <- modelCode({
##     a ~ dunif(-0.9999, 0.9999)
##     b ~ dnorm(0, sd = 1000)
##     sigPN ~ dunif(0.0001, 1)
##     sigOE ~ dunif(0.0001, 1)
##     x[1] ~ dnorm(b/(1-a), sd = sqrt(sigPN^2 + sigOE^2))
##     y[1] ~ dnorm(x[1], sd = sigOE)
##     for(i in 2:t){
##         x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
##         y[i] ~ dnorm(x[i], sd = sigOE)
##     }
## })



## better parameterization: mean and autocorrelation
code <- modelCode({
    mu ~ dnorm(0, sd = 1000)
    b ~ dnorm(0, sd = 1000)
    sigPN ~ dunif(0.0001, 1)
    sigOE ~ dunif(0.0001, 1)
##    x[1] ~ dnorm(mu, sd = sqrt((sigPN^2)/abs(1-a^2)))
    x[1] ~ dnorm(mu, sd = sqrt(sigPN^2 + sigOE^2))
    y[1] ~ dnorm(x[1], sd = sigOE)
    a <- 1-(b/mu)
    for(i in 2:t){
        x[i] ~ dnorm(x[i-1] * a + b, sd = sigPN)
        y[i] ~ dnorm(x[i], sd = sigOE)
    }
})



t <- 30
constants <- list(t = t)

Rmodel <- nimbleModel(code, constants = constants)
Rmodel$mu <- 10/(1-.5)
Rmodel$b <- 10
Rmodel$sigPN <- .1
Rmodel$sigOE <- .1

set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('mu','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))

data <- list(y = Rmodel$y)

inits <- list(mu = Rmodel$mu, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)








