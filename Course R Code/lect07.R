setwd("~/Kurser/02418/2021Fall/lect07")
##################################################
## rats (example 11.5)
rm(list=ls())
xdat<- scan('rat.dat',skip=1,
            what=list(group=0,surv=0,status=0))
xdat


## Non parametric analysis (Kaplan Meier)
library(survival)
dat <- data.frame(t=xdat$surv,dead=xdat$status,group=xdat$group)
Surv.Ex <- survfit(Surv(t, dead) ~ group, conf.type = "log-log",
                   data = dat)

## Kaplan-Meier plot
par(mfrow=c(1,1))
plot(Surv.Ex,mark.time=TRUE,col=1:2) ## Note mark time!

##################################################
## Exponential regression model
## likelihood function
nll.exp <- function(theta,t,d,X){
    mu <- exp(X %*% theta)
    - sum(d * dexp(t, rate = 1/mu, log=TRUE) + (1-d) *
          pexp(t, rate = 1/mu, lower.tail = FALSE, log.p=TRUE))
}

## Design matrix
X <- cbind(1,dat$group==2)
head(X)
t <- dat$t
d <- dat$dead

## test that it works
nll.exp(c(1,1),t=t,d=dat$dead,X)

## Find MLE
(opt.exp <- nlminb(c(0,0), nll.exp, t=t, d = d, X=X))


## Parameter unceartainty
library(numDeriv)
Hes <- hessian(nll.exp,opt.exp$par,t=t,d=d,X=X)
V <- solve(Hes)

## Confidence interval for parameters
CI <- cbind(opt.exp$par + qnorm(0.025) * sqrt(diag(V)),
      opt.exp$par + qnorm(0.975) * sqrt(diag(V)))

## Regression model in R
(mod1 <- survreg(Surv(t, dead) ~ factor(group), data = dat,
                 dist = "exponential"))
summary(mod1)
confint(mod1)


## Testing the effect

## LRT
nll.exp(1,t,d,X[ ,1]) ## does not work!
## forcing X[ ,1] to be a 1 column matrix
nll.exp(1,t,d,t(t(X[ ,1])))

## Optimise smaller model
opt0 <- nlminb(1, nll.exp, t=t, d = d, X=t(t(X[ ,1])))
chi2 <- 2*(opt0$objective-opt.exp$objective)
c(chi2, 1-pchisq(chi2,df=1)) ## Likelihood ratio test

## Profile likelihood
pnll.beta1 <- function(beta1,t,d,X){
    f.tmp <- function(beta0,beta1,t,d,X){
        nll.exp(c(beta0,beta1),t,d,X)
    }
    nlminb(1, f.tmp, beta1=beta1, t=t, d = d, X=X)$objective
}

se.beta1 <-  sqrt(diag(V))[2] ## Standard error for beta1
k <- 100
beta1 <- seq(opt.exp$par[2] + qnorm(0.01) * se.beta1,
             opt.exp$par[2] + qnorm(0.99) * se.beta1,
             length=k) ## The values to look at
pnll <- numeric(k) ## initalize vector
for(i in 1:length(beta1)){
     pnll[i]<-pnll.beta1(beta1[i], t = t, d = d, X = X)
}


## Profilelikelihood plot
plot(beta1, -(pnll-min(pnll)),type = "l")
lines(range(beta1), -qchisq(0.95,df=1)/2*c(1,1),col=2,lwd=2,lty=2)
confint(mod1)
## compare with wald interval
lines(confint(mod1)[2,1]*c(1,1),c(-3,0))
lines(confint(mod1)[2,2]*c(1,1),c(-3,0))
## so very close to the wald CI


##################################################
## Interpretating the results

## Time ratio
exp(c(coef(mod1)[2],confint(mod1)[2, ]))

## Hazard ratio
exp(-c(coef(mod1)[2],confint(mod1)[2,2:1 ]))

## So none of them are significiantly different from 1
## (on a 5% level)

##################################################
## Residual analysis
dat$CS.exp <- dat$t*exp(-mod1$linear.predictors) ## r_i
surv.exp <- survfit(Surv(CS.exp,dead==1)~1, data = dat) 
plot(surv.exp$time,-log(surv.exp$surv))
abline(a=0,b=1)
## so a rather bad fit!

##################################################
## Compare expected survival (under model) and KM
par(mfrow=c(1,2))
t <- seq(0,320,by=1)

h0.exp <- exp(-coef(mod1)[1])
h1.exp <- exp(-coef(mod1)[1]-coef(mod1)[2])
H0.exp <- t*h0.exp
H1.exp <- t*h1.exp


## Hazard
n <- length(t)
matplot(t,cbind(rep(h0.exp,n),rep(h1.exp,n)),col=1:2,
        type="l",lwd=2)

## Survival
plot(Surv.Ex,mark.time=TRUE,col=1:2)
matlines(t,exp(-cbind(H0.exp,H1.exp)),lty=2,lwd=2)
## So clearly insufficient

##################################################
## Weibull
##################################################


## Likelihood estimation
log.h.w <- function(theta, sigma, t, X){
    - log(sigma) + (1/sigma-1) * log(t) - X %*% theta /sigma
}

H.w <- function(theta,sigma,t,X){
    t^(1/sigma) * exp(- X %*% theta/sigma)
}

nll.w <- function(theta,t,d,X){
    p <- length(theta) - 1
    sigma <- theta[p+1]
    theta <- theta[1:p]
    - sum(d * log.h.w(theta, sigma, t, X) - H.w(theta, sigma, t, X))
}

    
nll.w(c(1,1,1),t=dat$t,d=dat$dead,X)

## MLE
opt.wei <- nlminb(c(1,1,1), nll.w, t=dat$t, d = dat$dead, X=X,
              lower=c(-Inf,-Inf,0.01))

opt.wei
opt.exp$objective

## Unceartainty estimation
Hes <- hessian(nll.w,opt.wei$par,t=dat$t,d=dat$d,X=X)
V <- solve(Hes)
(sd.par <- sqrt(diag(V)))

## profile for sigma
pnll.sigma <- function(sigma,t,d,X){
    f.tmp <- function(theta,sigma,t,d,X){
        nll.w(c(theta,sigma),t,d,X)
    }
    nlminb(c(1,1), f.tmp, sigma=sigma, t=dat$t, d = dat$d, X=X)$objective
}


sigma <- seq(opt.wei$par[3] - 3 * sd.par[3],
             opt.wei$par[3] + 3 * sd.par[3],
             length=100)
pnll <- sigma
for(i in 1:length(sigma)){
     pnll[i]<-pnll.sigma(sigma[i], t = t, d = d, X = X)
}

par(mfrow=c(1,1))
plot(sigma, -(pnll-min(pnll)),type = "l")
lines(range(sigma), -qchisq(0.95,df=1)/2*c(1,1),
      col=2,lwd=2,lty=2)

## Test of the weibull against the exponential
(chi2<--2*(opt.wei$objective-pnll.sigma(1, t = t, d = d, X = X)))

1-pchisq(chi2,df=1)

    
## Directly in R
(mod2 <- survreg(Surv(t, dead) ~ factor(group), data = dat,
                 dist = "weibull"))
summary(mod2)
opt.wei$par

(mod1 <- survreg(Surv(t, dead) ~ factor(group), data = dat,
                  dist="exponential"))
summary(mod1)

anova(mod1,mod2) ## Likelihood ratio test

##################################################
## Dianostic plot
par(mfrow=c(1,1))
dat$CS.wei <- exp((log(dat$t)-mod2$linear.predictors)/mod2$scale)
surv.wei <- survfit(Surv(CS.wei, dead==1)~1 , data = dat)
plot(surv.wei$time, -log(surv.wei$surv))
abline(a=0, b=1)
### So better than the exponential model

##################################################
## Compare expected survival (under model) and KM
par(mfrow=c(1,2))
t <- seq(0,320,by=1)
n <- length(t)
h0.wei <- exp(log.h.w(opt.wei$par[1:2], opt.wei$par[3],
                      t, X= cbind(1,rep(0,n))))
h1.wei <- exp(log.h.w(opt.wei$par[1:2], opt.wei$par[3],
                      t, X= cbind(1,rep(1,n))))
H0.wei <- H.w(opt.wei$par[1:2], opt.wei$par[3],
              t, X= cbind(1,rep(0,n)))
H1.wei <- H.w(opt.wei$par[1:2], opt.wei$par[3],
              t, X= cbind(1,rep(1,n)))


## Hazard
matplot(t,cbind(h0.wei,h1.wei),col=1:2, type="l",lwd=2)
matplot(t,cbind(h0.wei/h1.wei),col=1:2, type="l",lwd=2)

## Survival
plot(Surv.Ex,mark.time=TRUE,col=1:2)
matlines(t,exp(-cbind(H0.wei,H1.wei)),lty=2,lwd=2)
## So clearly much better than exponential


##################################################
## log-logistic model
mod3 <- survreg(Surv(t, dead) ~ factor(group), data = dat,
                dist = "loglogistic")
summary(mod3)

## Likelihood estimation (based on hazad functions)
log.h.log <- function(theta, sigma, t, X){
    mu <- X%*%theta
    z <- (log(t) - mu) / sigma
    - log(sigma) + z - log(1+exp(z)) - log(t)
}

H.log <- function(theta,sigma,t,X){
    mu <- X%*%theta
    z <- (log(t) - mu) / sigma
    log(1+exp(z))
}

nll.log <- function(theta,t,d,X){
    p <- length(theta) - 1
    sigma <- theta[p+1]
    theta <- theta[1:p]
    - sum(d * log.h.log(theta, sigma, t, X) - H.log(theta, sigma, t, X))
}

nll.log(c(1,1,1),t=dat$t,d=dat$dead,X)
## MLE
opt.log <- nlminb(c(1,1,1), nll.log, t=dat$t, d = dat$dead, X=X,
              lower=c(-Inf,-Inf,0.01))

opt.log
opt.log$par
mod3

-opt.log$objective
-opt.wei$objective
## log-logistic slightly better

## log-logis diagnostic
par(mfrow=c(1,1))
dat$z <- (log(dat$t) - mod3$linear.predictors)/mod3$scale
dat$CS.log <- log(1+exp(dat$z))
surv.log <- survfit(Surv(CS.log, dead==1)~1 , data = dat)
plot(surv.log$time, -log(surv.log$surv))
abline(a=0, b=1)


## Kaplan-Meier plot
par(mfrow=c(1,2))
t <- seq(0,350,by=1)

## Hazard
z0 <- (log(t) - coef(mod3)[1])/mod3$scale
z1 <- (log(t) - coef(mod3)[1]-coef(mod3)[2])/mod3$scale
h0 <- 1/(1 + exp(-z0)) * 1 / mod3$scale * 1/t
h1 <- 1/(1 + exp(-z1)) * 1 / mod3$scale * 1/t

h0 <- exp(log.h.log(coef(mod3), mod3$scale, t, cbind(1,rep(0,length(t)))))
h1 <- exp(log.h.log(coef(mod3), mod3$scale, t, cbind(1,rep(1,length(t)))))

H0 <- exp(-H.log(coef(mod3), mod3$scale, t, cbind(1,rep(0,length(t)))))
H1 <- exp(-H.log(coef(mod3), mod3$scale, t, cbind(1,rep(1,length(t)))))


par(mfrow=c(1,2))
## plot hazard
n <- length(t)
matplot(t,cbind(h0, h1),col=1:2,lty=2,lwd=2,type="l")


## plot survival
plot(Surv.Ex,col=1:2)
matlines(t,cbind(H0,H1),col=1:2, lty=2, lwd=2,type="l")
## Not clear if weibull or log-logistic is better

c(AIC(mod1),AIC(mod2),AIC(mod3))

anova(mod1,mod2,mod3)



## Can we now see the effect of groups?
mod3R <- survreg(Surv(t, dead) ~ 1, data = dat,
                dist = "loglogistic")
anova(mod3R,mod3)

##################################################
## WHAS
##################################################

par(mfrow=c(1,1))
WHAS <- read.delim("WHAS.txt")
WHAS$YearFol <- WHAS$lenfol/365.25
Surv.Bysex <- survfit(Surv(YearFol, fstat == 1) ~ gender, 
                     conf.type = "log-log", data = WHAS)
plot(Surv.Bysex, conf.int = FALSE, las = 1, xlab = "Years since admission", 
     ylab = "Estimated Survival Prob.", col=2:3, lwd = 2, mark.time=TRUE)
legend("bottomleft", col = 2:3, c("Male","Female"), lwd = 2)


##################################################
## Exponential model 
mod1 <- survreg(Surv(YearFol, fstat) ~ gender, data = WHAS,
                dist = "exponential")
summary(mod1)
confint(mod1)


##################################################
## Time ratio
exp(c(coef(mod1)[2],confint(mod1)[2, ]))

##################################################
## Hazard ratio
exp(-c(coef(mod1)[2],confint(mod1)[2, ]))

##################################################
## A more complicated model

mod2 <- survreg(Surv(YearFol, fstat) ~ gender + age + bmi, data = WHAS,
                dist = "exponential")
summary(mod2)
confint(mod2)

##################################################
##
##################################################

##################################################
## Time ratio
(TR.age <- exp(c(coef(mod2)["age"],confint(mod2)["age", ])))
(TR.bmi <- exp(c(coef(mod2)["bmi"],confint(mod2)["bmi", ])))
(TR.gender <- exp(c(coef(mod2)["gender"],confint(mod2)["gender", ])))

confint(mod1)
confint(mod2)
## notice that gender is no longer significant
summary(WHAS[WHAS$gender==1,c(6:9)])
summary(WHAS[WHAS$gender==0,c(6:9)])
##################################################

##################################################
## Weibull 
##################################################


weib <- function(x,sigma,beta0){1/sigma*(x^(1/sigma-1))*exp(-beta0/sigma)}
curve(weib(x,2,1),0,8, xlab="Time", ylab="h(t)", col=1, lwd=2)
curve(weib(x,0.8,1),0,8, add=TRUE, col=2, lwd=2)
curve(weib(x,1,1),0,8, add=TRUE, col=3, lwd=2)
curve(weib(x,0.25,1),0,8, add=TRUE, col=4, lwd=2)
legend("topright", col=1:4, c("sigma=2","sigma=0.8","sigma=1","sigma=0.25"), lwd=2)


mod3 <- survreg(Surv(YearFol, fstat) ~ gender + age + bmi,
                data = WHAS, dist = "weibull")
summary(mod3)

confint(mod3)


mod4 <- survreg(Surv(YearFol, fstat) ~  age + bmi,
                data = WHAS, dist = "weibull")
summary(mod4)



##################################################
## Log logistic
##################################################


mod4 <- survreg(Surv(YearFol, fstat) ~ gender + age + bmi,
                data = WHAS, dist = "loglogistic")
summary(mod4)

## Time ratio
exp(cbind(coef(mod4)[2],confint(mod4)[2, 1], confint(mod4)[2, 2]))

## Hazard ratio
exp(coef(mod4)[2]/mod4$scale)

##################################################
## Dianostic plot
WHAS$z <- (log(WHAS$YearFol)-mod4$linear.predictors)/mod4$scale
WHAS$CS4 <- log(1+exp(WHAS$z))
surv4 <- survfit(Surv(CS4, fstat==1)~1 , data = WHAS)
plot(surv4$time, -log(surv4$surv))
abline(a=0, b=1, lty=2)


## Estimated survival functions
xrange <- range(WHAS$YearFol)
t <- seq(xrange[1],xrange[2],length=500)
#Exponential
coef2 <- mod2$coefficients
S21 <- exp(-t/exp(coef2[1]+coef2[2]+coef2[3]*71+coef2[4]*27))
S20 <- exp(-t/exp(coef2[1]+coef2[3]*71+coef2[4]*27))
#WEIBULL
coef3 <- mod3$coefficients
z31 <- (log(t)-(coef3[1]+coef3[2]+coef3[3]*71+coef3[4]*27))/mod3$scale
z30 <- (log(t)-(coef3[1]+coef3[3]*71+coef3[4]*27))/mod3$scale
S31 <- exp(-exp(z31))
S30 <- exp(-exp(z30))
#LOG LOGISTIC
coef4 <- mod4$coefficients
z41 <- (log(t)-(coef4[1]+coef4[2]+coef4[3]*71+coef4[4]*27))/mod4$scale
z40 <- (log(t)-(coef4[1]+coef4[3]*71+coef4[4]*27))/mod4$scale
S41 <- (1+exp(z41))^-1
S40 <- (1+exp(z40))^-1

# get the range for the x and y axis
yrange <- range(S41)

# set up the plot
plot(xrange, yrange, type="n", xlab="Years since admission",
     ylab="Probability of survival (age=71, BMI=27)",las=1) 
#PLOT THE SURVIVAL FUNCTIONS
lines(t, S41, type="l", col=1, lty=2, lwd=2)
lines(t, S40, type="l", col=2, lwd=2)
lines(t, S31, type="l", col=3 , lty=2, lwd=2)
lines(t, S30, type="l", col=4, lwd=2)
lines(t, S21, type="l", col=5, lty=2, lwd=2)
lines(t, S20, type="l", col=6, lwd=2)
legend(x="topright", lwd=2, col=1:6, 
       legend=c("LL Fem","LL Male","Wei Fem","Wei Male",
                "Exp Fem","Exp Male"))
##################################################
##
##################################################
