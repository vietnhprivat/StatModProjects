rm(list=ls())
setwd("~/Kurser/02418/2024Fall/lect10/")
##################################################
## AR(1) model

## Simulation:

## Parameters 
x0 <- rnorm(1)
phi <- 0.8
sigma <- 1
n <- 100
set.seed(2454)
## Simulation
x <- numeric(n+1)
x[1] <- x0
for(i in 1:n){
    x[i+1]<-phi*x[i]+rnorm(1,sd=sigma)
}

plot(x,type="l")
acf(x)

## Moment estimates
(phi.hat <- acf(x,plot=FALSE)$acf[2]) ## ACF in lag 1
cor(x[-1],x[-length(x)])
(sigma <- var(x)*(1-phi.hat^2))


## Likelihood estimation
nll <- function(theta,y){
    n <- length(y) - 1
    n/2 * log(theta[1]) + 1/(2*theta[1]) * sum((y[-1]-theta[2]*y[-(n+1)])^2)
}



## MLE
(opt <- nlminb(c("sigma"=1,"phi"=1/2),nll,y=x,lower=c(0.001,-0.999),upper=c(Inf,0.999)))

## Check MLEs
(phi <- sum(x[-1]*x[-length(x)])/sum(x[-length(x)]^2))
(sigma.sq <- 1/(length(x)-1) * sum((x[-1]-phi * x[-length(x)])^2))

## standard errors
library(numDeriv)
V <- solve(hessian(nll,opt$par,y=x))
(se <- sqrt(diag(V)))

## Profile likelihood for phi
llp.phi <- function(phi,x,x0){
    n <- length(x) - 1
    -n/2 * log(sum((x[-1]-phi*x[-(n+1)])^2))
}

## Plot profile likelihood
phi <- seq(max(opt$par[2]-3*se[2],-1),
           min(opt$par[2]+3*se[2],1),
           length=200)

llp <- sapply(phi,llp.phi,x=x[-1],x0)

plot(phi,exp(llp-max(llp)),type="l")
lines(range(phi),
      exp(-c(1,1) * qchisq(0.95,df=1)/2),
          col=2,lty=2)

## Directly in R
arima(x,order=c(1,0,0),include.mean=FALSE)
opt$par ## rather close. 

## Did it help?
e <- x[-1]-opt$par[2]*x[-n]
acf(e)

##################################################
## Mixture models
##################################################

## Data (including year 2007-2016)
dat <- read.table("earthquakes.txt",header=FALSE)
names(dat) <- c("year","eq")

## Plot data
par(mfrow=c(1,2))
plot(dat,type="b",xlab="Year",ylab="Count",
     main="No. of major earthquakes")
tab <- table(dat$eq)
plot(tab,xlim=c(0,45),axes=FALSE,xlab="Count",ylab="Freq")
points(0:45,dpois(0:45,lambda=mean(dat$eq))*dim(dat)[1],pch=19)
mean(dat$eq)
var(dat$eq)
## hence overdispersion

par(mfrow=c(1,1))
acf(dat$eq)
## hence autocorreltaion

## looking into overdispersion
summary(fit1 <- glm(eq~1,data=dat,family=poisson))
summary(glm(eq~1,data=dat,family=quasipoisson))

summary(fit1)
## Goodness of fit test
1-pchisq(285.8,df=116)
## Hence we need either overdispersion or other model.

## Poisson mixture: transform
## natural to working parameters
pois.mix.pn2pw <- function(n.dist,lambda,delta){
    if(sum(delta) >= 1){print("sum(delta) should be < 1")
    return()}
    if(length(lambda) != n.dist){
        print("length(lambda) should be n.dist")
        return()
    }
    if(length(delta) != (n.dist-1)){
        print("length(delta) should be n.dist-1")
        return()
    }
    eta <- log(lambda)
    tau <- log(delta/(1-sum(delta)))
    return(list(eta=eta,tau=tau))
}

## Poisson mixture: transform
## working to natural parameters
pois.mix.pw2pn <- function(n.dist,eta,tau){
    if(n.dist==1){return(exp(eta))}
    if(length(eta) != n.dist){
        print("length(lambda) should be n.dist")
    return()}
    if(length(tau) != (n.dist-1)){
        print("length(delta) should be n.dist-1")
        return()
    }
    lambda <- exp(eta)
    delta <- exp(tau)/(1+sum(exp(tau)))
    delta <- c(1-sum(delta),delta)
    return(list(lambda=lambda,delta=delta))
}


## negative log likelihood function
nll <- function(theta,n.dist,y){
    if(n.dist==1){
        return(-sum(dpois(y,lambda=exp(theta),log=TRUE)))
    }
    eta <- theta[1:n.dist]
    tau <- theta[(n.dist+1):(2*n.dist-1)]
    n.pars <- pois.mix.pw2pn(n.dist,eta,tau)
    n <- length(y)
    nll <- 0
    for(i in 1:n){
        nll <- nll - log(sum(
                         n.pars$delta *
                         dpois(y[i],
                               lambda=n.pars$lambda)))
    }
    return(nll)
}

## Estimation with one distribution
n.dist <- 1; ## No of states
## Initial values
lambda <- mean(dat$eq); 
delta <- c()
## Working parameters
wpars <- pois.mix.pn2pw(n.dist,lambda,delta)
theta <- c(wpars$eta,wpars$tau)
## MLE
(opt1 <- nlminb(theta,nll,n.dist=n.dist,y=dat$eq))
## Natural parameters
(npars1 <- pois.mix.pw2pn(n.dist, opt1$par[1:n.dist],
                          opt1$par[(n.dist+1) :
                          (2*n.dist-1)]))

## Estimation with two distributions
n.dist <- 2; ## No of states
## Initial values
lambda <- mean(dat$eq)*c(1/2,3/2);
delta <- c(0.5)
## Working parameters
wpars <- pois.mix.pn2pw(n.dist,lambda,delta)
theta <- c(wpars$eta,wpars$tau)
## MLE
(opt2 <- nlminb(theta,nll,n.dist=n.dist,y=dat$eq))
## Natural parameters
(npars2 <-
     pois.mix.pw2pn(n.dist, opt2$par[1:n.dist],
                    opt2$par[(n.dist+1):(2*n.dist-1)]))

## Estimation with 3 distributions
n.dist <- 3;
## Initial values
lambda <- mean(dat$eq)*c(1/2,1,3/2);
delta <- c(1/3,1/3)
wpars <- pois.mix.pn2pw(n.dist,lambda,delta)
theta <- c(wpars$eta,wpars$tau)
## MLE
(opt3 <- nlminb(theta,nll,n.dist=n.dist,y=dat$eq))
## Natural parameters
(npars3 <-
     pois.mix.pw2pn(n.dist, opt3$par[1:n.dist],
                    opt3$par[(n.dist+1):(2*n.dist-1)]))

## Estimation with 4 distributions
n.dist <- 4;
## Initial values
lambda <- mean(dat$eq)*c(1/4,3/4,5/4,7/4);
delta <- c(1,1,1)/4
wpars <- pois.mix.pn2pw(n.dist,lambda,delta)
theta <- c(wpars$eta,wpars$tau)
## MLE
(opt4 <- nlminb(theta,nll,n.dist=n.dist,y=dat$eq))
## natural parameters
(npars4 <-
    pois.mix.pw2pn(n.dist, opt4$par[1:n.dist],
                   opt4$par[(n.dist+1):(2*n.dist-1)]))

## Plot the result
plot(tab/sum(tab),xlim=c(0,45),axes=FALSE,
     xlab="Count",ylab="Density")
axis(1);axis(2)
x <- 0:45

## Small function to calculate mixture dist
mix.dist <- function(x,npars){
    sum(npars$delta*dpois(x,lambda=npars$lambda))
}
##################################################

lines(x,dpois(x,lambda=npars1),type="b",col=1,pch=19)
lines(x,sapply(x,mix.dist,npars=npars2),type="b",col=2,pch=19)
lines(x,sapply(x,mix.dist,npars=npars3),type="b",col=3,pch=19)
lines(x,sapply(x,mix.dist,npars=npars4),type="b",col=4,pch=19)

## Choose model
(AIC <- 2 * c(opt1$objective, opt2$objective,
              opt3$objective, opt4$objective) +
     2 * c(length(opt1$par), length(opt2$par),
           length(opt3$par), length(opt4$par)))
## Hence we choose 3-state model
AIC.pois.mix <- AIC

##################################################
## likelihood ratio test 
1-pchisq(2*(opt1$objective-opt2$objective),
       df=length(opt2$par)-length(opt1$par))

1-pchisq(2*(opt2$objective-opt3$objective),
       df=length(opt3$par)-length(opt2$par))

1-pchisq(2*(opt3$objective-opt4$objective),
       df=length(opt4$par)-length(opt3$par))
## Hence same conclusion as for AIC

##################################################
## CI for working parameters
library(numDeriv)

## Present the result
H <- hessian(nll,opt3$par,n.dist=3,y=dat$eq)
V.pars <- solve(H)
round(V.pars,digits=3)
se <- sqrt(diag(V.pars))
tab <- round(cbind(opt3$par,se),digits=2)
rownames(tab) <- c("eta1","eta2","eta3",
                   "tau1","tau2")
tab

## Confidence intervals for natural pars
## by simulation from normal dist.
library(mvtnorm)
k <- 100000
PARS <- rmvnorm(k,mean=opt3$par,sigma=V.pars)
dim(PARS)
## Simulated (lambda1)
(CIlambda1 <- quantile(exp(PARS[ ,1]),probs=c(0.025,0.975)))
## Wald based
exp(opt3$par[1]+qnorm(c(0.025,0.975))*se[1])

## Simulated (lambda2)
CIlambda1 <- rbind(CIlambda1,quantile(exp(PARS[ ,2]),probs=c(0.025,0.975)))
CIlambda1[2, ]

## Wald based
exp(opt3$par[2]+qnorm(c(0.025,0.975))*se[2])

## Simulated (lambda3)
CIlambda1 <- rbind(CIlambda1,quantile(exp(PARS[ ,3]),probs=c(0.025,0.975)))
CIlambda1[3, ]
## Wald based
exp(opt3$par[3]+qnorm(c(0.025,0.975))*se[3])

## Simulated delta1-delta3
delta2 <- exp(PARS[ ,4])/(1+rowSums(exp(PARS[ ,4:5])))
delta3 <- exp(PARS[ ,5])/(1+rowSums(exp(PARS[ ,4:5])))
delta1 <- 1-delta2-delta3

## Estimated values of delta 1-3
delta2.hat <- exp(opt3$par[4]) /
    (1+sum(exp(opt3$par[4:5])))
delta3.hat <- exp(opt3$par[5]) /
    (1+sum(exp(opt3$par[4:5])))
delta1.hat <- 1-delta2.hat-delta3.hat

## CI for delta
CIdelta1 <- c(delta1.hat,quantile(delta1,probs=c(0.025,0.5,0.975)))
CIdelta1 <- rbind(CIdelta1,
                  c(delta2.hat,quantile(delta2,
                                        probs=c(0.025,0.5,0.975))))
CIdelta1 <- rbind(CIdelta1,
                  c(delta3.hat,quantile(delta3,
                                        probs=c(0.025,0.5,0.975))))
round(CIdelta1, digits = 2)

## Parametric boot strap
delta.hat <- c(delta1.hat,delta2.hat,delta3.hat)
lambda <- exp(opt3$par[1:3])

## Sample from distribution
rmix.pois <- function(k,delta,lambda){
  states <- sample(1:length(delta),size=k,replace=TRUE,
                   prob=delta)
  rpois(k,lambda[states])
}

## Tjek that is works
k <- 100000
y <- rmix.pois(k,npars3$delta,npars3$lambda)
plot(table(y)/k,type="h")
x <- seq(0,60)
lines(x, sapply(x, mix.dist, npars = npars3),
      type = "b", col = 3, pch = 19)


## Bootstrap
lambda <- mean(dat$eq)*c(1/2,1,3/2);
delta <- c(1/3,1/3)
n.dist <- 3
wpars <- pois.mix.pn2pw(n.dist,lambda,delta)
theta <- c(wpars$eta,wpars$tau)
n <- length(dat$eq)
PARS2 <- c();j<-1
for(i in 1:200){ ## Should choose a larger number
  y <- rmix.pois(n,npars3$delta,npars3$lambda)
    pars <- nlminb(theta,nll,n.dist=n.dist,y=y)$par
    n.pars <- pois.mix.pw2pn(3,pars[1:3],pars[4:5])
    I <- sort(n.pars$lambda,index.return=TRUE)$ix
    n.pars$lambda <- n.pars$lambda[I]
    n.pars$delta <- n.pars$delta[I]
    PARS2 <-
        rbind(PARS2, c(n.pars$lambda,n.pars$delta))
}

round(CIlambda1, digits = 2)
## Boot strap CI for lambda
round(rbind(quantile(PARS2[ ,1],prob=c(0.025,0.975)),
            quantile(PARS2[ ,2],prob=c(0.025,0.975)),
            quantile(PARS2[ ,3],prob=c(0.025,0.975))),
      digits = 2)

## Boot strap CI for delta
round(CIdelta1, digits = 2)
round(rbind(quantile(PARS2[,4],prob=c(0.025,0.975)),
            quantile(PARS2[,5],prob=c(0.025,0.975)),
            quantile(PARS2[,6],prob=c(0.025,0.975))),
      digits = 2)



##################################################
## Markov models
##################################################
## data
  y <- c(2,3,3,2,1,1,1,1,1,2,
       3,1,3,2,3,3,2,1,2,2,
       3,2,3,2,3,3,2,2,2,2,
       3,1,3,2,3,3,2,2,1,2,
       3,2,3,2,1,3,2,2,3,2,
       3,1,3,2,3,3,2,2,2,3,
       3,2,3,2,3,3,1,2,3,2,
       3,2,3,2,3,3,1,2,2,2,
       3,2,3,2,1,3,2,1,2,3,
       3,1,3,2,3,3,2,1,2,1)

## Plot data
plot(y,type="s")
points(1:100,y,pch=19,cex=0.5)

## Observed transtions
f <- matrix(ncol=3,nrow=3)
f
n <- length(y)
f[1,1] <- sum(y[-n]==1 & (y[-1]-y[-n])==0)
f[1,2] <- sum(y[-n]==1 & (y[-1]-y[-n])==1)
f[1,3] <- sum(y[-n]==1 & (y[-1]-y[-n])==2)

f[2,1] <- sum(y[-n]==2 & (y[-1]-y[-n])==-1)
f[2,2] <- sum(y[-n]==2 & (y[-1]-y[-n])==0)
f[2,3] <- sum(y[-n]==2 & (y[-1]-y[-n])==1)

f[3,1] <- sum(y[-n]==3 & (y[-1]-y[-n])==-2)
f[3,2] <- sum(y[-n]==3 & (y[-1]-y[-n])==-1)
f[3,3] <- sum(y[-n]==3 & (y[-1]-y[-n])==0)
f

## MLE
Gamma <- f
Gamma[1, ] <- Gamma[1, ]/sum(y[-n]==1)
Gamma[2, ] <- Gamma[2, ]/sum(y[-n]==2)
Gamma[3, ] <- Gamma[3, ]/sum(y[-n]==3)
round(Gamma, digits = 2)



## stationary dist
U <- matrix(1, ncol = 3, nrow = 3)
I <- diag(3)
(delta <- c(1, 1, 1) %*% solve(I - Gamma + U))

delta %*% Gamma-delta ## It is a stationary distribution

## Log likelihood
(l <- sum(f*log(Gamma)))

## Testing H0: independence
Gamma0 <- matrix(ncol=3,nrow=3)
Gamma0[ ,1] <- sum(y[-1]==1)/(length(y)-1)
Gamma0[ ,2] <- sum(y[-1]==2)/(length(y)-1)
Gamma0[ ,3] <- sum(y[-1]==3)/(length(y)-1)
Gamma0

l0 <- sum(f*log(Gamma0))
(Q<-2*(l-l0))
1-pchisq(Q,df=6-2)
## So not independence 

## test using contingency tables
f0 <- Gamma0 
f0[1, ] <- f0[1, ] * sum(y[-1]==1)
f0[2, ] <- f0[2, ] * sum(y[-1]==2)
f0[3, ] <- f0[3, ] * sum(y[-1]==3)

Q <- sum((f-f0)^2/f0)
1-pchisq(Q,df=2*2)
## So same conclusion
##################################################
## End
##################################################
