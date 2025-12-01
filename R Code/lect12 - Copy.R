rm(list=ls())
##################################################
## Ex 10.1 (my own simulation)
set.seed(2134)
N <- 20 ## sample size
(mu <- rnorm(N,sd=5)) ## different mu


## observations:
y1 <- rnorm(N,mu)  
y2 <- rnorm(N,mu)

## Averages
(ybar <- (y1+y2)/2)
(RSS <- sum((y1-ybar)^2)+sum((y2-ybar)^2))

## Profile likelihood for sigma^2
l.sig <- function(sig, y1, y2, ybar){
    RSS <- sum((y1-ybar)^2)+sum((y2-ybar)^2)
    N <- length(y1)
    -N*log(sig)-1/(2*sig)*RSS
}

## Plot the profile likelihood
sig <- seq(0.1,2,length=200)
lp.sig <- sapply(sig,l.sig,y1=y1,y2=y2,ybar=ybar) 
plot(sig,exp(lp.sig-max(lp.sig)),type="l")
## We see a strong bias (sigma^2 is 1, but we estimate about 1/2)

## "true" likelihood (in real applications we do not have
## this) 
lp.true <- sapply(sig,l.sig,y1=y1,y2=y2,ybar=mu) ## We insert the true mean
lines(sig,exp(lp.true-max(lp.true)),col=2)
lines(range(sig),c(1,1)*exp(-qchisq(0.95,df=1)/2),lty=2,lwd=2,col=3)

## Estimates
RSS/(2*N) ## Biased 
RSS/(N) ## unbiased

## Anova... as we know it
y <- c(y1,y2)
b <- rep(1:N,2)
anova(lm(y~as.factor(b))) ## The unbiased estimate
##################################################


##################################################
## Ex 10.2

## Simulate data
set.seed(132)
N <- 60
beta0 <- -1
beta1 <- 1
tau <- rep(c(0,1),each=N) ## Treatment effect
(S <- rnorm(N,sd=3))  ## Individual (random) effect
S <- rep(S,2)

## logit <- beta0+S+tau
p <- exp(beta0+S+tau*beta1)/(1+exp(beta0+S+tau*beta1))
y <- rbinom(2*N,size=1,prob=p)

## summarize data
tab <- cbind(c(0,1),c(sum(y[tau==0]),sum(y[tau==1])),N)
colnames(tab) <- c("tau", "y","N")
tab

## Set up model 
b <- rep(1:N,2)
b
tau

## Estimate
fit <- glm(y~factor(tau)+factor(b),family="binomial")
dim(summary(fit)$coef)
summary(fit)$coef[1:2, ]
## True values beta0 = -1 beta1 = 1

## Profile likelihood confint.
confint(fit)[1:2, ] ## We see that the
## estimate highly biased (true tau is 1)
c(beta0,beta1)

plot(coef(fit))
##################################################


##################################################
## Example 10.4 cont.
## Simulation
set.seed(2354)
n <- 10
x <- rnorm(n)


##  likelihood
lp.sig <- function(sig,x){
    n <- length(x)
    -(n)/2 * log(sig) - 1/(2*sig) * (n-1) * var(x)
}

## plot profile log-likelihood 
sig <- seq(0.1,5,length=100)
lp <- sapply(sig,lp.sig,x=x)
plot(sig,lp-max(lp),type="l",ylim=c(-4,0))


## orthogonal likelihood 
l.orth.sig <- function(sig,s2,n){
    -(n-1)/2 * log(sig) - (n-1) * s2 /(2*sig)
}

## plot result
l.orth <- sapply(sig,l.orth.sig,s2=var(x),n=n)
lines(sig,l.orth-max(l.orth),col=2)
lines(range(sig),-c(1,1)*qchisq(0.95,df=1)/2,lty=2,lwd=2,col=3)

## The usual CI for variances
rug((n-1)*var(x)*c(1/qchisq(0.975,df=n-1),1/qchisq(0.025,df=n-1)),
    lwd=3,col="blue")

## plot likelihood
plot(sig,exp(lp-max(lp)),type="l",ylim=c(0,1))

l.orth <- sapply(sig,l.orth.sig,s2=var(x),n=n)
lines(sig,exp(l.orth-max(l.orth)),col=2)
lines(range(sig),exp(-c(1,1)*qchisq(0.95,df=1)/2),lty=2,lwd=2,col=3)
rug((n-1)*var(x)*c(1/qchisq(0.975,df=n-1),1/qchisq(0.025,df=n-1)),
    lwd=2,col=4)
##################################################


##################################################
## Ex 10.1 (my own simulation), cont..

## Just same simulation as before
set.seed(2134)
N <- 20
mu <- rnorm(N,sd=5)
y1 <- rnorm(N,mu)
y2 <- rnorm(N,mu)
ybar <- (y1+y2)/2
RSS <- sum((y1-ybar)^2)+sum((y2-ybar)^2)
########################################


## profile likelihood (as above)
l.sig <- function(sig, y1, y2, ybar){
    RSS <- sum((y1-ybar)^2)+sum((y2-ybar)^2)
    N <- length(y1)
    -N*log(sig)-1/(2*sig)*RSS
}

sig <- seq(0.1,2,length=200)
lp.sig <- sapply(sig,l.sig,y1=y1,y2=y2,ybar=ybar) 
plot(sig,exp(lp.sig-max(lp.sig)),type="l")
## "true" likelihood
lp.true <- sapply(sig,l.sig,y1=y1,y2=y2,ybar=mu)
lines(sig,exp(lp.true-max(lp.true)),col=2)
lines(range(sig),c(1,1)*exp(-qchisq(0.95,df=1)/2),lty=2,lwd=2,col=3)

## Marginal likelihood
l.sig.marg <- function(sig, v){
    N <- length(v)
    -N/2*log(sig)-1/(2*sig)*sum(v^2)
}

lp.marg <- sapply(sig,l.sig.marg,v=(y1-y2)/sqrt(2))
lines(sig,exp(lp.marg-max(lp.marg)),col=4)



## Anova...
y <- c(y1,y2)
b <- rep(1:N,2)
anova(lm(y~as.factor(b)))
rug(anova(lm(y~as.factor(b)))[2,3],lwd=2,col="blue")
##################################################


##################################################
## Example 4.4 (section 10.5)
y1 <- 5; n1 <- 15
y2 <- 1; n2 <- 10

(p1 <- y1 / n1); (p2 <- y2 /n2)



## Likelihood for profile likelihood
ll1 <- function(p2,theta,n1,n2,y1,y2){
    p1 <- exp(theta) * p2 / (1 + p2 * (exp(theta) - 1))
    log(choose(n1,y1)) + log(choose(n2,y2)) + y1*log(p1) +
        (n1-y1)*log(1-p1) +
        y2*log(p2) + (n2 - y2) * log(1-p2)
}

## Profile likelihood
lp <- function(theta,n1,n2,y1,y2){
    ll <- optimise(ll1, c(0,1),theta=theta, n1=n1, n2=n2, y1=y1,
                   y2=y2, maximum=TRUE)
    ll$objective
}

## Plot the profilelikelihood
theta <- seq(-1,7,by=0.1)
lP<- theta
for(i in 1:length(theta)){
    lP[i]<-lp(theta[i],n1,n2,y1,y2)
}

plot(theta,exp(lP-max(lP)),type="l")

## likelihood for conditional likelihood
l <- function(theta,n1,n2,y1,y2){
    s <- seq(0,y1+y2)
    choose(n1,y1)*choose(n2,y2)*exp(theta*y1)/
        sum(choose(n1,s)*choose(n2,y1+y2-s)*exp(theta*s))
}

## plot the result
ll<-theta
for(i in 1:length(theta)){
    ll[i] <- l(theta[i],n1,n2,y1,y2)
}


lines(theta,ll/max(ll),type="l",col=2)
lines(range(theta),exp(-qchisq(0.95,df=1)/2)*c(1,1),
      col=4,lty=2,lwd=2)
lines(c(0,0),c(0,1))
##################################################

##################################################
## Ex 10.2 (our own data)

## simulation
set.seed(132)
N <- 60
beta0 <- -1
tau <- rep(c(0,1),each=N)
S <- rnorm(N,sd=3)
S <- rep(S,2)

logit <- beta0+S+tau

p <- exp(beta0+S+tau)/(1+exp(beta0+S+tau))
y <- rbinom(2*N,size=1,prob=p)
y

## set up the model
b <- rep(1:N,2)
b
tau

## fit the model 
fit <- glm(y~factor(tau)+factor(b),family="binomial")
summary(fit)$coef[1:2, ]
confint(fit)[1:2, ]

## set up data for conditional likelihood 
y1 <- y[1:N]
y2 <- y[(N+1):(2*N)]

n1 <- sum((y1-y2)==1) ## y1==1 and y2==0
n2 <- sum((y1-y2)==-1)  ## y1==0 and y2==1
n1
n2

## conditional likliehood
tau <- seq(0,5,length=200)
l.cond <- n1*log(1/(1+exp(tau)))+n2*log(1/(1+exp(-tau)))

## plot likelihood
plot(tau,exp(l.cond-max(l.cond)),type="l")
## cut off
lines(range(tau),exp(-qchisq(0.95,df=1)/2)*c(1,1),lty=2,lwd=2,col=2)
lines(c(1,1),c(0,1)) ## True value
## The estimate
lines(log(n2/n1)*c(1,1),c(0,1))

## ignoring the strata effect
tau <- rep(c(0,1),each=N)
fit2 <- glm(y~factor(tau),family="binomial")
summary(fit2)
confint(fit2)
## we see a quite different CI (shifted to the left), and too confident

##################################################
## The end.... 
##################################################
