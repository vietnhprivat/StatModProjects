## Prior
p <- seq(-0.01,1.01,length=1000)
plot(p,dbeta(p,shape1=1,shape2=1),type="l",ylim=c(0,15))

## Agla
n1 <- 252
x1 <- 194
(p.hat <- 194/252)

## Posterior
alpha <- x1+1
beta <- n1 - x1 + 1
lines(p,dbeta(p,shape1 = alpha, shape2 = beta),
      type = "l", ylim = c(0, 2), col = 2)
rug(p.hat,col=2,lwd=3)

##################################################
## Confidence intervals
##################################################

## Baysian CI
qbeta(c(0.025,0.975),shape1=alpha,shape2=beta)
## Wald CI
p.hat + c(-1, 1) * qnorm(0.975) * sqrt(p.hat*(1-p.hat)/n1)


##################################################
## Combine experiments
##################################################

## Our experiment (from lecture 1)
n2 <- 8
x2 <- 7

## Combine experiment
alpha <- x1 + x2 + 1
beta <- n1 + n2 - x1 - x2 +1

(p.hat <- (x1+x2)/(n1+n2))

lines(p,dbeta(p,shape1 = alpha, shape2 = beta),
      type = "l", ylim = c(0, 2), col = 4)
rug(p.hat,col=4,lwd=3)

## Baysian CI
qbeta(c(0.025,0.975),shape1=alpha,shape2=beta)
## Wald CI
p.hat + c(-1, 1) * qnorm(0.975) * sqrt(p.hat*(1-p.hat)/(n1+n2))


##################################################
## Baysian prediction
##################################################

## Beta binomial model 
dbeta.binom <- function(n, k, a, b){
    choose(n,k) * beta(a + x, b + n - x) / beta(a, b)
}


n <- 1
x <- 0:1

## Compare beta-binomial distribution and binomial
dbeta.binom(n, x, alpha, beta)
dbinom(x, size = n, prob = p.hat)


n <- 10
x <- 0:10

plot(x, dbinom(x, size = n, prob = p.hat), type="h")
lines(x+0.1, dbeta.binom(n, x, alpha, beta), type="h",col=2)

##################################################
## Example to play with
##################################################
n1 <- 10
x1 <- 7
alpha <- x1 + 1
beta <- n1 - x1 + 1
p.hat <- x1/n1

plot(x, dbinom(x, size = n, prob = p.hat), type="h")
lines(x+0.1, dbeta.binom(n, x, alpha, beta), type="h",col=2)

##################################################
## End
##################################################
