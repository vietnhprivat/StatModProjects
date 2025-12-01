rm(list=ls())
setwd("~/Kurser/02418/2024Fall/lect06/")
## Example 11.5:
xdat<- scan('rat.dat',skip=1,
            what=list(group=0,surv=0,status=0))
xdat

xdat$surv[xdat$group==1]
xdat$surv[xdat$group==2]


## A quick look at data (ignoring censoring)
par(mfrow=c(1,1))
plot(ecdf(xdat$surv[xdat$group==1]),xlim=range(xdat$surv))
lines(ecdf(xdat$surv[xdat$group==2]),col=2)


##################################################
## likelihood estimation: Exponential model

## Exponential model 
nll.exp <- function(theta, delta, y){
    - sum(dexp(y[delta==1],rate=1/theta,log=TRUE)) - 
        sum(log(1-pexp(y[delta==0],rate=1/theta))) 
    ## Negative log-likelihood
}


## Eg. group 1:
optimize(nll.exp,c(0,600),
         delta=xdat$status[xdat$group==1],
         y= (xdat$surv[xdat$group==1]))

## commom mean for both models
(OPT0 <- optimize(nll.exp,c(0,600),
                 delta=xdat$status,
                 y= (xdat$surv)))

## Different means in the two models
nll.exp2 <- function(theta, delta, y, group){
     nll.exp(theta[1], delta[group==1], y[group==1]) +
        nll.exp(theta[1]+theta[2], delta[group==2], y[group==2])
}


(OPT1 <- nlminb(c(OPT0$minimum,0),nll.exp2,
               delta=xdat$status,
               y= (xdat$surv),group=xdat$group))

## Parameters
OPT0$minimum
OPT1$par

## negative log-likelihood
OPT0$objective
OPT1$objective


## LRT
chi2 <- - 2 * ((OPT1$objective)-OPT0$objective) ## test statistics
1 - pchisq(chi2, df = 1) ## p-value
## Hence no difference

##################################################
## Parameter unceartainty

## Quadratic approximation
library(numDeriv)
H <- hessian(nll.exp2,OPT1$par,delta=xdat$status,
             y= (xdat$surv),group=xdat$group)
sqrt(diag(solve(H)))
## CI for difference
(CI.diff <- OPT1$par[2] + c(-1,1) * qnorm(0.975) * sqrt(diag(solve(H)))[2])

## Profile the parameter
pll <- function(theta2, delta, y, group){
    f.tmp <- function(theta1, theta2, delta, y, group){
        nll.exp2(c(theta1,theta2), delta, y, group)
    }
    optimize(f.tmp,c(0,350),theta2=theta2, delta = delta,
             y = y, group = group)$objective
}


theta2 <- seq(-180,280)
matplot(theta2, -(sapply(theta2,pll,delta=xdat$status,
                       y= (xdat$surv),group=xdat$group) -
                  OPT1$objective),type="l",ylab="ll")
lines(range(theta2), - c(1,1)*qchisq(0.95,df=1)/2,col=2,lty=2)
rug(CI.diff)
CI.diff ## Wald interval (quadratic approximation)
## The profile likelihood is read off the plot


## Model check (is the exponential assumption reasonable?)
## Compare with fig 11.3 in IAL
n1 <- sum(xdat$group==1 & xdat$status==1)
a1 <- -log(1-seq(1/(2*n1),by=1/n1,len=n1))
diff(pexp(a1))
x1 <- xdat$surv[xdat$group==1 & xdat$status==1]
plot(sort(a1),sort(x1), xlab='Exponential quantile',
     ylab='Ordered uncensored data',pch="1",ylim = c(0,max(x1)))
n2 <- sum(xdat$group==2 & xdat$status==1)
a2 <- -log(1-seq(1/(2*n2),by=1/n2,len=n2))
x2 <- xdat$surv[xdat$group==2 & xdat$status==1]
points(sort(a2),sort(x2),pch='2')
## Does not seem exponential! (should follow straight
## lines from (0,0))
##################################################
##
##################################################


##################################################
## Kaplan Meier estimator
## Rat example

## Start by group 1
surv <- xdat$surv[xdat$group==1]
status <- xdat$status[xdat$group==1]

I <- sort(surv,index.return=TRUE)$ix
cbind(surv,status)[I, ]
## Kaplan Meier (by "hand" calculations)
(n <- length(surv)) ## no of rats
status <- status[sort(surv,index.return=TRUE)$ix]
t <- sort(surv) ## time of events

## setting up the table
t.new <- t[c(TRUE,diff(t)!=0)]
n2 <- length(t.new)
j <- 1:n2
dt <- numeric(n2) ## death at time t(j)
ct <- numeric(n2) ## Censor at time t(j)
for(i in 1:n2){
    dt[i] <- sum(status[t==t.new[i]])
    ct[i] <- sum(1-status[t==t.new[i]])
}

Rt <- n - c(0,cumsum((dt + ct)[-n2])) ## number at risk
S <- cumprod(1-dt/Rt) ## Survival function
round(cbind(j,t.new,Rt,d=dt,c = ct, CProb= 1-dt/Rt, S = S),digits=2)


## Kaplan-Meier plot
plot(t.new,S,type="s",ylim=c(-0.05,1.025),xlim=c(100,350))

## Variance of logS (and naive CI)
VlogS <- cumsum(dt/(Rt*(Rt-dt)))
VS <- VlogS * S^2
lines(t.new,S + qnorm(0.975) * sqrt(VS),lty=2,type="s")
lines(t.new,S - qnorm(0.975) * sqrt(VS),lty=2,type="s")
## Note how lines cross 1 and 0

## Directly in R
library(survival)
dat <- data.frame(t=xdat$surv,dead=xdat$status,group=xdat$group)
dat
Surv.Ex <- survfit(Surv(t, dead==1)~1, conf.type = "plain",
                   data = dat[dat$group==1, ])
summary(Surv.Ex) ## R simply truncate
plot(Surv.Ex)


## Kaplan-Meier plot
plot(t.new,S,type="s",ylim=c(0,1))

## Unceartainty using log-log
VloglogS <- 1/(log(S))^2 * cumsum(dt/(Rt*(Rt-dt)))
cu <- log(-log(S)) +  qnorm(0.975) * sqrt(VloglogS)
cl <- log(-log(S)) -  qnorm(0.975) * sqrt(VloglogS)
lines(t.new,exp(-exp(cu)), lty=2,col=2,type="s")
lines(t.new,exp(-exp(cl)), lty=2,col=2,type="s")
## intervals in (0,1)

## Using R (log-log)
Surv.Ex <- survfit(Surv(t, dead==1)~1, conf.type = "log-log",
                   data = dat[dat$group==1, ])
summary(Surv.Ex)
plot(Surv.Ex,xlim=c(100,320))



## Estimating quantiles
lines(c(100,320), 0.5 * c(1,1),col=2,lty=2)
lines(c(100,320), 0.25 * c(1,1),col=3,lty=2)
lines(c(100,320), 0.75 * c(1,1),col=3,lty=2)

quantile(Surv.Ex,c(0.25,0.5,0.75))

## 1 - survival function
plot(Surv.Ex,fun=function(x){1-x})



##################################################
## Comparing survival functions
##################################################
by(dat$t,dat$group,summary)
## overlap in IQR

Surv.Bygroup <- survfit(Surv(t,dead == 1) ~ group, conf.type = "log-log",
                        data = dat)
Surv.Bygroup ## CI from KM plot

## plot without confidence interval
## Small overlap, but are they different?

plot(Surv.Bygroup, conf.int = FALSE, las = 1, xlab = "Time", 
     ylab = "Estimated Survival Prob.", col=2:3, lwd = 2,xlim=c(100,350),ylim=c(0,1))
legend("bottomleft", col = 2:3, c("group 1","group 2"), lwd = 2)
lines(c(100,350), c(1,1)*0.5)

## plot with confidence interval
plot(Surv.Bygroup, conf.int = TRUE, las = 1, xlab = "Time", 
     ylab = "Estimated Survival Prob.", col=2:3, lwd = 2,xlim=c(100,350))
legend("bottomleft", col = 2:3, c("group 1","group 2"), lwd = 2)
lines(c(100,350), c(1,1)*0.5)
##################################################



##################################################
## Testing the difference
survdiff(Surv(t,dead == 1) ~ group, data = dat)


##################################################
## The test by hand calculaions (skip...)
dat
Surv.Ex1 <- survfit(Surv(t, dead==1)~1, data = dat[dat$group==1,])
Surv.Ex2 <- survfit(Surv(t, dead==1)~1, data = dat[dat$group==2,])
Surv.Ex <- survfit(Surv(t, dead==1)~1, data = dat)

summary(Surv.Ex1)
summary(Surv.Ex2)
summary(Surv.Ex)
## Are the two first different form the common surv func.
names(Surv.Ex)

## Setting up a table for comparing the functions
## (i.e. calculate Q)
tab <- cbind(t = Surv.Ex$time, Ri=Surv.Ex$n.risk,
             di=Surv.Ex$n.event,R1i=Surv.Ex1$n.risk[1],
             d1i=0, E1i = 0, R2i=0,d2i = 0, E2i = 0)
## fill in for group 1
for(i in 1:length(Surv.Ex1$time)){
    I1 <- Surv.Ex$time >= Surv.Ex1$time[i] 
    tab[I1,"R1i"] <- Surv.Ex1$n.risk[i] ## number at risk at time ti
    I1 <- Surv.Ex$time==Surv.Ex1$time[i]
    tab[I1,"d1i"] <- Surv.Ex1$n.event[i] ## number dead at time ti
    I1 <- Surv.Ex$time>Surv.Ex1$time[i]
    tab[I1,"R1i"] <- Surv.Ex1$n.risk[i] - Surv.Ex1$n.event[i]
    ## number at risk after time ti
}
## calculating the number in group 2
tab[ ,"R2i"] <- tab[ ,"Ri"]-tab[ ,"R1i"]
tab[ ,"d2i"] <- tab[ ,"di"]-tab[ ,"d1i"]

## calculate the expected number at each time (and group)
tab[ ,"E1i"] <- tab[ ,"R1i"] * tab[ ,"di"] / tab[ ,"Ri"]
tab[ ,"E2i"] <- tab[ ,"R2i"] * tab[ ,"di"] / tab[ ,"Ri"]

## Total expected death in each group
E1 <- sum(tab[ ,"E1i"])
E2 <- sum(tab[ ,"E2i"])

## Total number of dead in each group
O1 <- sum(tab[ ,"d1i"])
O2 <- sum(tab[ ,"d2i"])


cbind(c(O1,O2), c(E1,E2), c((E1-O1)^2/E1, (E2-O2)^2/E2))

survdiff(Surv(t,dead == 1) ~ group, data = dat)

## Calculation the variance using wi=1
v <- tab[ ,"R1i"] * tab[ ,"R2i"] * tab[ ,"di"] *
    (tab[ ,"Ri"]-tab[ ,"di"]) /
    (tab[ ,"Ri"]^2 * (tab[ ,"Ri"]-1))

v[is.na(v)]<-0
sum(v)

(Q <- round((sum(tab[ ,"E1i"]-tab[ ,"d1i"]))^2/sum(v),digits=4))
round((sum(tab[ ,"E2i"]-tab[ ,"d2i"]))^2/sum(v),digits=4)

s<- survdiff(Surv(t,dead == 1) ~ group, data = dat,rho=0)
s
## P-value
1-pchisq(Q,df=1)

##################################################
## Harrington and Flemming:
## Survival function
S <- c(1,Surv.Ex$surv[-length(Surv.Ex$surv)])
S
round((sum(S*(tab[ ,"E1i"]-tab[ ,"d1i"])))^2/sum(S^2*v),digits=4)

survdiff(Surv(t,dead == 1) ~ group, data = dat,rho=1)
## Same conclusion (different p-value)
##################################################

##################################################
## WHAS
##################################################
setwd("~/Kurser/02418/2022Fall/lect06/")

WHAS <- read.delim("WHAS.txt")
head(WHAS, n = 10)
dim(WHAS)


##  I CANNOT SUBTRACT THEM AS THEY ARE!
WHAS$foldate- WHAS$admitdate


## NEED TO CONVERT TO DAY FORMAT
WHAS$InDate <- as.Date(WHAS$admitdate, "%d-%m-%Y")
WHAS$OutDate <- as.Date(WHAS$foldate, "%d-%m-%Y")
WHAS$DaysDif <- WHAS$OutDate-WHAS$InDate
head(WHAS[ ,c(2:3,5,10:12)], n = 5)
WHAS$YearFol <- WHAS$lenfol/365.35

## Summarizing the survival time
by(WHAS$YearFol, WHAS$gender, summary)

library(survival)
## survival analysis
Surv.WHAS <- survfit(Surv(YearFol, fstat == 1)~1, 
                     conf.type = "log-log",data = WHAS)
summary(Surv.WHAS)

## Graphical presentation:
## Total surv function
par(mfrow=c(1,1))
plot(Surv.WHAS, conf.int = FALSE, las = 1, xlab = "Years since admission", 
     main="Kaplan-Meier plot",ylab = "Estimated Survival Prob.")

## Including CI
plot(Surv.WHAS, conf.int = T, las = 1, xlab = "Years since admission", 
     main="Kaplan-Meier plot",ylab = "Estimated Survival Prob.")


abline(h = 0.5, col = "green", lwd = 2) ## Median survival time
abline(v = 6, col = "green", lwd = 2) ## how many survived 6 years?


plot(Surv.WHAS, conf.int = T, las = 1, xlab = "Years since admission", 
     main="Kaplan-Meier plot",ylab = "Estimated Survival Prob.")
abline(h = 0.75, col = "blue", lwd = 2) ## lower quartile (of dead)
abline(v = c(1.4726, 2.1185), col = 5:6, lwd = 2) ## any number in this
                                                  ## interval eg midpoint

## Median survival including CI
plot(Surv.WHAS, conf.int = T, las = 1, xlab = "Years since admission", 
     main="Kaplan-Meier plot",ylab = "Estimated Survival Prob.")
abline(h = 0.5, col = "green", lwd = 2)
abline(v = 6.02, col = "green", lwd = 2)
abline(v = c(4.45, 7.42), col = 6, lwd = 2, lty=2)

## calculated by R
quantile(Surv.WHAS, c(0.25, 0.50, 0.75))


#PLOT 1-KM Cumulative failure function
plot(Surv.WHAS, conf.int = F, fun=function(x) { 1- x }, las = 1, 
     xlab = "Years since admission", 
     ylab = "Estimated Failure Prob.", lwd = 2)

## Compare number dead and censored
by(WHAS$fstat, WHAS$gender, sum) ## number dead (not censored)
by(WHAS$fstat==0, WHAS$gender, sum) ## number censored

## Comparing survival function
Surv.Bysex <- survfit(Surv(YearFol, fstat == 1) ~ gender, 
                     conf.type = "log-log", data = WHAS)
Surv.Bysex

## Graphically
plot(Surv.Bysex, conf.int = FALSE, las = 1, xlab = "Years since admission", 
     ylab = "Estimated Survival Prob.", col=2:3, lwd = 2)
legend("bottomleft", col = 2:3, c("Male","Female"), lwd = 2)

    ## including CI
plot(Surv.Bysex, conf.int = TRUE, las = 1, xlab = "Years since admission", 
     ylab = "Estimated Survival Prob.", col=2:3, lwd = 2)
legend("bottomleft", col = 2:3, c("Male","Female"), lwd = 2)


## Test (log-rank)
survdiff(Surv(YearFol, fstat == 1) ~ gender, data = WHAS)

## More wheight on early survival (HF)
survdiff(Surv(YearFol, fstat == 1) ~ gender, data = WHAS, 
         rho = 1)
##################################################
## End !!
##################################################
