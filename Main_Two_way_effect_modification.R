################################################################################
#### Interactive effects of air pollution and air temperature on mortality #####
#### Using PM10 and total mortality in Chicago during 1987-2000 as an example ##
#### Kai Chen, December 2017                                               #####
#### Institute of Epidemiology, Helmholtz Zentrum Muenchen                 #####
################################################################################
##load the pacakges
library(mgcv); library(dlnm); library(splines)

##set the working direction
## setwd("Your own working folder")
###Load the function crossreduce() in crossreduce_int.R using 2 AP_cats
source("crossreduce_int_2APcats.R")

##load the predefined function to calculate the modification effect by temperature using 3 categories
source("Modbytemp_3cat.R")
## Modbytemp_3cat <-function(formula,data,poll,temp_cat,incr,pp,...){...}
  
#####################
### Preparation #####
#####################
##(a)Load the city data
data <- chicagoNMMAPS

##(b)Define the formula using penalized DLNM
form.bytemp <- death~ns(time,df=8*14)+dow+cbgam
form.byAP <- death~ns(time,df=8*14)+dow+AP+cbgam*AP_cat

##(c)Define the crossbasis using Penalized DLNM, based on Gasparrini et al. Biometrics 2017
# PS: EXCLUDE THE DEFAULT PENALTY FROM THE LAG-RESPONSE FUNCTION WITH fx=T
dfPS <- 7
maxlag <- 21
argvar <- list(fun="ps",df=dfPS-1)
cbgam <- crossbasis(data$temp,lag=maxlag,argvar=argvar,arglag=list(fun="ps", df=dfPS, fx=T))
##onebasis without lag
b1temp <- onebasis(data$temp, fun="ps", df=dfPS-1)

# Define the doubly varying penalty matrices
# Lag using varying difference penalty
C <- do.call('onebasis',c(list(x=0:maxlag,fun="ps",df=dfPS,intercept=T)))
D <- diff(diag(maxlag+1),diff=2)
P <- diag((seq(0,maxlag-2))^2)
Slag1 <- t(C) %*% t(D) %*% P %*% D %*% C
# Coefficients using varying ridge penalty
Slag2 <- diag(rep(0:1,c(3,4)))
cbgamPen <- cbPen(cbgam,addSlag=list(Slag1,Slag2))
PP <- list(cbgam=cbgamPen)

##(d)Define the temp_cat based on 25th and 75th percentile
Tempmean25th <- quantile(data$temp, 0.25, na.rm = T)
Tempmean75th <- quantile(data$temp, 0.75, na.rm = T)
data$temp_cat <- ifelse(data$temp <= Tempmean25th, 1, 
                        ifelse(data$temp < Tempmean75th, 2, 3))

##(e)Define AP using ozone, and set AP_cat for 2 categories
data$AP <- data$pm10
data$AP_cat <- as.factor(ifelse(data$AP <= median(data$AP, na.rm = T), 1, 2))

#############################
###I. Modification by temp ##
#############################
result.pm10 <- Modbytemp_3cat(form.bytemp, data=data, poll=pm10, temp_cat =temp_cat, incr=10,pp=PP)
result.pm10 #Highest effect at high temperatures

#############################
###II. Modification by AP  ##
#############################
model2 <- gam(form.byAP, data=data, family=quasipoisson,paraPen=PP, method='REML') ##Using 'REML' as AG's paper
summary(model2)

# Reduction to the overall estimates
# Using the modified function crossreduce_int for interaction analysis with AP_cat (3 categories)
red <- crossreduce_int_2APcats(cbgam,model2, cen=21) ##using 21°C as the cen value
#cat1
coef.APcat1 <- red$coef_APcat1
vcov.APcat1 <- red$vcov_APcat1
#cat2
coef.APcat2 <- red$coef_APcat2
vcov.APcat2 <- red$vcov_APcat2

##prediction
cp.APcat1 <- crosspred(b1temp, coef=coef.APcat1, vcov=vcov.APcat1, 
                       model.link = "log", cen=21)
cp.APcat2 <- crosspred(b1temp, coef=coef.APcat2, vcov=vcov.APcat2, 
                       model.link = "log", cen=21)

#plot the overall CRF
plot(cp.APcat2, "overall", ylab="Overall cumulative RR", xlab="Temperature",lwd=2, 
     ci.arg=list(density=20,col="salmon"), col="red", ylim=c(0.5,2.5), main="Total mortality")
lines(cp.APcat1,"overall",ci="area",ylim=c(0.5,2.5), lwd=2, ci.arg=list(density=20,angle=-45,col="dodgerblue"), col="blue")
title("Chicago", cex=2,adj=0)
legend("top",c("Lowly polluted", "Highly polluted"),col=c("blue","red"),lwd=2, bty = "n", 
       y.intersp = 1.5, title="PM10")

###Calculate the heat effect (29°C vs 21°C)
quantile(data$temp, 0.99, na.rm = T) ##28.6°C
##Low AP
cp.APcat1$allRRfit[names(cp.APcat1$allRRfit)=="29"] 
cp.APcat1$allRRlow[names(cp.APcat1$allRRlow)=="29"] 
cp.APcat1$allRRhigh[names(cp.APcat1$allRRhigh)=="29"] 
##high AP
cp.APcat2$allRRfit[names(cp.APcat2$allRRfit)=="29"] 
cp.APcat2$allRRlow[names(cp.APcat2$allRRlow)=="29"] 
cp.APcat2$allRRhigh[names(cp.APcat2$allRRhigh)=="29"] 

###Calculate the cold effect (-15°C vs 21°C)
quantile(data$temp, 0.01, na.rm = T) ##-15°C
##Low AP
cp.APcat1$allRRfit[names(cp.APcat1$allRRfit)=="-15"] 
cp.APcat1$allRRlow[names(cp.APcat1$allRRlow)=="-15"] 
cp.APcat1$allRRhigh[names(cp.APcat1$allRRhigh)=="-15"] 

##High AP
cp.APcat2$allRRfit[names(cp.APcat2$allRRfit)=="-15"] 
cp.APcat2$allRRlow[names(cp.APcat2$allRRlow)=="-15"] 
cp.APcat2$allRRhigh[names(cp.APcat2$allRRhigh)=="-15"] 


