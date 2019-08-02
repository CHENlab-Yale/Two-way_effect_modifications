#############################################################################################################
### Function for calculating the air pollution effects modified by temperature categories 
### Used in single-pollutant QuasiPoisson Generalized Additive Models
### Dr. Kai Chen and Dr. Regina Hampel in 2017. 
### Function description:
#Calculation of air pollution effects modified by temperature categories (<1. Quartile vs. 1.-3.Quartile vs. >3.Quartile)
#Calculation of percent changes

#formula: formula with confounders, e.g. death~ ns(trend, df) + dow
#data: analysis data set
#poll: vector with pollutants, an interaction between air pollutant and temperature will be calculated,  e.g. c(m24pm25_lag0, m24pnc_lag1)
#temp_cat: temperature in categories 1 (<1.quartile), 2 (1.-3.quartile), 3 (>3 quartile)

#effects are presented as percent change
#incr: air pollution increment
#pp: paraPen
#############################################################################################################

##Predefined function to calculate the modification effect by temperature using 3 categories
Modbytemp_3cat <-function(formula,data,poll,temp_cat,incr,pp,...){
  
  polleff<-matrix(rep(NA,3*7),ncol=7) 
  colnames(polleff)<-c("estimate","se","p-value","increment","Pchange","Plower", "Pupper")
  rnames<-vector()
  increment<-rep(incr,3)
  
  #model without interaction
  mod<-gam(formula=formula, data=data,na.action=na.omit,family=quasipoisson,paraPen = pp, method='REML',...)
  ##set the pollname and tempname(with lag structures in the names)
  pollname <- substitute(poll)
  tempname<-  substitute(temp_cat)
  
  #include interaction
  tempfac<-paste("as.factor(",tempname,")",sep="")
  interaction<-paste(tempfac,pollname,sep="*")
  fff<-as.formula(paste(paste(".~.", interaction, sep="+"),"+1",sep=""))
  mod2<-update(mod,fff)
  
  #calculate exposure effects for temperature groups
  coeff<-as.data.frame(mod2$coeff)
  p<-which(row.names(coeff)==pollname)
  
  cov<-summary(mod2)$cov.scaled
  index<-which(row.names(cov)==pollname)
  
  #cat=1
  polleff[(1),1]<-coeff[p,1]
  var1<-cov[index,index]
  polleff[(1),2]<-sqrt(var1)
  #cat=2
  polleff[(2),1]<-coeff[p+1,1]+coeff[p,1]
  polleff[(2),2]<-sqrt(var1+cov[index+1,index+1]+cov[index,index+1]+cov[index+1,index])
  #cat=3
  polleff[(3),1]<-coeff[p+2,1]+coeff[p,1]
  polleff[(3),2]<-sqrt(var1+cov[index+2,index+2]+cov[index,index+2]+cov[index+2,index])
  
  #increments
  polleff[(1:3),4]<-increment
  
  #calculation of single p-values
  polleff[,3]<-round((1-pnorm(abs(polleff[,1]/polleff[,2])))*2,digits=9)
  
  #calculation of percent changes
  
  polleff[,5]<-(exp(polleff[,1]*polleff[,4])-1)*100
  polleff[,6]<-(exp((polleff[,1]-1.96*polleff[,2])*polleff[,4])-1)*100
  polleff[,7]<-(exp((polleff[,1]+1.96*polleff[,2])*polleff[,4])-1)*100
  
  #output file
  polleff2<-as.data.frame(polleff)
  polleff2$pollutant <- rep(as.character(pollname),3)
  polleff2$city<-rep(data$city[1],3)
  polleff2$outcome<-rep(as.character(formula[[2]]),3)
  polleff2$temp_cat<-c('Low temperature','Medium temperature', 'High temperature')
  
  ##### results
  return(polleff2)
}

