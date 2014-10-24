library(nimble)


################
### tester
################

code <- modelCode({
    x[1:5] ~ dmnorm(mu[1:5], cov=Sigma[1:5,1:5])
    y[1:5] ~ dmnorm(x[1:5], cov=Ident[1:5,1:5])
})

mu <- rep(0,5)

Chol <- matrix(0, nrow=5, ncol=5)
Chol[1,1:5] <- c(1,   .99, .9, -.8, 0)
Chol[2,1:5] <- c(0, 1,   .8, 0,   0)
Chol[3,1:5] <- c(0,  0,  1,  -.9, 0)
Chol[4,1:5] <- c(0, 0,  0,  1,   0)
Chol[5,1:5] <- c(0,   0,   0,  0,   1)
Sigma <- t(Chol) %*% Chol
Prec <- solve(Sigma)
cov2cor(Sigma)
Ident <- diag(5)
y <- c(3, 3, 2, 0, -4)
postPrec <- Prec + Ident
postSigma <- solve(postPrec)
postCor <- cov2cor(postSigma)
postMean <- postSigma %*% t(t(y))

constants <- list(mu=mu, Sigma=Sigma, Ident=Ident)
data <- list(y = y)
inits <- list()

abtester <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)


################
### litters
################

G <- 2
N <- 16
n <- array(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 11, 8, 10, 13, 10, 12, 9, 10, 9, 10, 5, 9, 9, 13, 7, 5, 10, 7, 6, 10, 10, 10, 7), dim = c(2, 16))
r <- array(c(13, 12, 12, 11, 9, 10, 9, 9, 8, 10, 8, 9, 12, 9, 11, 8, 9, 8, 9, 4, 8, 7, 11, 4, 4, 5, 5, 3, 7, 3, 7, 0), dim = c(2, 16))
p <- array(0.5, dim = c(2, 16))

constants <- list(G=G, N=N, n=n)
data      <- list(r=r)
inits     <- list(p=p)

code <- modelCode({
     for (i in 1:G) {
         a[i] ~ dgamma(1, 0.001)
         b[i] ~ dgamma(1, 0.001)
         # mu[i] <- a[i] / (a[i] + b[i])
         # theta[i] <- 1 / (a[i] + b[i])
         for (j in 1:N) {
             r[i,j] ~ dbin(p[i,j], n[i,j])
             p[i,j] ~ dbeta(a[i], b[i])
         }
     }
})

ablitters <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)



################
### SSMmub
################

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

abSSMmub <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)


################
### SSMab
################

## parameterization in terms of slope (autocorelation) and intercept,
## which are highly correlated in mixing
code <- modelCode({
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
constants <- list(t = t)
Rmodel <- nimbleModel(code, constants = constants)
Rmodel$a <- .5
Rmodel$b <- 10
Rmodel$sigPN <- .1
Rmodel$sigOE <- .1
set.seed(0)
calculate(Rmodel, Rmodel$getDependencies(c('a','b','sigPN','sigOE'), determOnly = TRUE))
simulate(Rmodel, Rmodel$getDependencies(c('x', 'y')))
data <- list(y = Rmodel$y)
inits <- list(a = Rmodel$a, b = Rmodel$b, sigPN = Rmodel$sigPN, sigOE = Rmodel$sigOE, x = Rmodel$x)

abSSMab <- autoBlock(code=code, constants=constants, data=data, inits=inits, control=control)

