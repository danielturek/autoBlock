library(nimble)
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

