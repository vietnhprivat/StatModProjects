##################################################
## Example: Multinomial
rm(list=ls())
y <- c(5, 5, 5)
y <-y*10

nll <- function(theta, y){
    p <- c(theta, 1 - sum(theta))
    - sum(y * log(p))
}

theta <- c(0.1, 0.1)
nll(theta, y)


?nlminb
opt <- nlminb(c(0.25, 0.25), nll, y = y, 
              lower = c(0, 0), upper = c(1, 1))
opt
y[1:2] / sum(y)

## Plotting
p1 <- seq(0,1, by = 0.01)
p2 <- seq(0,1, by = 0.01)
logL <- matrix(ncol = length(p2), nrow = length(p1))
for(i in 1:length(p1)){
    for(j in 1:length(p2)){
        logL[i, j] <- -nll(c(p1[i], p2[j]), y)
    }
}

warnings()

par(mfrow=c(1,2))
xlim <- c(0,1); ylim <- xlim
## xlim <- c(0.2,0.5); ylim <- xlim

contour(exp(logL + opt$objective),
        levels=c(0.01,0.05,0.1,0.25,0.5,0.75),xlim=xlim,ylim=ylim)
lines(c(0,1),c(1,0),col=2,lwd=2)
lines(c(0,0),c(1,0),col=2,lwd=2)
lines(c(0,1),c(0,0),col=2,lwd=2)

##################################################
## quadratic app.
##################################################
library(numDeriv)
?hessian
I <- hessian(nll, opt$par, y = y)
I

## Compare numerical results with theoretical result
n <- sum(y)
n^2*(1 / y[1] + 1/y[3])
n^2 / y[3]

## Plotting
logLapp <- matrix(ncol = length(p2), nrow = length(p1))
for(i in 1:length(p1)){
    for(j in 1:length(p2)){
        theta <- c(p1[i], p2[j])
        logLapp[i,j] <- -0.5 * as.numeric(
           (theta - opt$par) %*% I %*% (theta - opt$par))
    }
}

contour(exp(logLapp), levels = c(0.01,0.05,0.1,0.25,0.5,0.75),
        xlim = xlim, ylim = ylim)
lines(c(0, 1), c(1, 0), col = 2, lwd = 2)
lines(c(0, 0), c(1, 0), col = 2, lwd = 2)
lines(c(0, 1), c(0, 0), col = 2, lwd = 2)


#################################################################
## Profile likelihood
#################################################################
lp.p1 <- function(p1, y){
    fun.tmp <- function(p2, p1, y){
        nll(c(p1,p2), y = y)
    }
    optimise(fun.tmp,c(0,1 - p1), p1 = p1, y = y)$objective
}

lp.p1(0.1, y)
p1 <- seq(0.01, 0.99, by = 0.01)
logLp1 <- sapply(p1, lp.p1, y = y) ## note sapply!
logLp1 <- logLp1 - min(logLp1) ## normalization

par(mfrow = c(1,1))
plot(p1, exp(-logLp1), type = "l")
lines(c(0,1), exp(-qchisq(0.95,df=1)/2) * c(1,1), col = 2)


##################################################
## Curvature of profile likelihood (side note p. 63)
I
hessian(lp.p1,opt$par[1],y=y)
1/solve(I)[1,1]

##################################################
## Compare Profile likelihood and Wald CI
rug(opt$par[1] + c(-1, 1) * qnorm(0.975) * sqrt(solve(I)[1,1]),
    lwd=2,col="red")


##################################################
# End.......
##################################################
