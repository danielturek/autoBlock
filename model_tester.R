library(nimble)


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






