#ST7016 Final Project

#min-max scaling function
normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

#multivariate normal distribution posterior draws
rmvnorm<-
  function(n,mu,Sigma) {
    p<-length(mu)
    res<-matrix(0,nrow=n,ncol=p)
    if( n>0 & p>0 ) {
      E<-matrix(rnorm(n*p),n,p)
      res<-t(  t(E%*%chol(Sigma)) +c(mu))
    }
    res
  }

#generating stationarity plots for assessing the convergence
stationarity.plot<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  
  boxplot(x~group,...)               }

#Import the data
abalone<-read.csv("~/Desktop/ST7016_Final_Project/abalone.data", header = FALSE, col.names = c("sex", "length","diameter","height", "whole weight", "shucked weight","viscera weight","shell weight", "rings"))

#A brief view of dataset
str(abalone)
summary(abalone)
dim(abalone)

#To check if there is null values
sum(is.na(abalone))

#To transfer sex into a quantitative variable
abalone$sex<-factor(abalone$sex)
abalone$sex<-unclass(abalone$sex)
abalone$sex<-as.numeric(abalone$sex)
abalone$rings<-as.numeric(abalone$rings)

#To dectect the spread of the response variable
table(abalone$rings)

#To scale the transformed dataset
aba<-normalize(abalone)

#To briefly detect the correlationships between predictors and found the collinearity exists
cor(aba)

#Referring to the existed library from lecture materials
source("~/Desktop/regression_gprior.R")

index <- sample(1:nrow(aba), size = 0.4 * nrow(aba))
train<-aba[index,]
test<-aba[-index,]

Y<-train[,9]
n<-length(Y)
X<-train[,-9]
X<-cbind(rep(1,n),X)
X<-as.matrix(X)
colnames(X)<-c("Intercept","X1","X2","X3","X4","X5","X6","X7","X8")

Y.te<-test[,9]
X.te<-test[,-9]
X.te<-cbind(rep(1,length(Y.te)),X.te)
X.te<-as.matrix(X.te)
colnames(X.te)<-c("Intercept","X1","X2","X3","X4","X5","X6","X7","X8")

#main effects
X2<-X[,3]
X3<-X[,4]
X4<-X[,5]
X5<-X[,6]
X6<-X[,7]
X7<-X[,8]
X8<-X[,9]

#two-way interactions

X2X3<-X2*X3
X2X4<-X2*X4
X2X5<-X2*X5
X2X6<-X2*X6
X2X7<-X2*X7
X2X8<-X2*X8

X3X4<-X3*X4
X3X5<-X3*X5
X3X6<-X3*X6
X3X7<-X3*X7
X3X8<-X3*X8

X4X5<-X4*X5
X4X6<-X4*X6
X4X7<-X4*X7
X4X8<-X4*X8

X5X6<-X5*X6
X5X7<-X5*X7
X5X8<-X5*X8

X6X7<-X6*X7
X6X8<-X6*X8

X7X8<-X7*X8

X<-cbind(X,X2X3,X2X4,X2X5,X2X6,X2X7,X2X8,X3X4,X3X5,X3X6,X3X7,X3X8,X4X5,X4X6,X4X7,X4X8,X5X6,X5X7,X5X8,X6X7,X6X8,X7X8)
X<-as.matrix(X)

#main effects
Xte2<-X.te[,3]
Xte3<-X.te[,4]
Xte4<-X.te[,5]
Xte5<-X.te[,6]
Xte6<-X.te[,7]
Xte7<-X.te[,8]
Xte8<-X.te[,9]

#two-way interactions

Xte2Xte3<-Xte2*Xte3
Xte2Xte4<-Xte2*Xte4
Xte2Xte5<-Xte2*Xte5
Xte2Xte6<-Xte2*Xte6
Xte2Xte7<-Xte2*Xte7
Xte2Xte8<-Xte2*Xte8

Xte3Xte4<-Xte3*Xte4
Xte3Xte5<-Xte3*Xte5
Xte3Xte6<-Xte3*Xte6
Xte3Xte7<-Xte3*Xte7
Xte3Xte8<-Xte3*Xte8

Xte4Xte5<-Xte4*Xte5
Xte4Xte6<-Xte4*Xte6
Xte4Xte7<-Xte4*Xte7
Xte4Xte8<-Xte4*Xte8

Xte5Xte6<-Xte5*Xte6
Xte5Xte7<-Xte5*Xte7
Xte5Xte8<-Xte5*Xte8

Xte6Xte7<-Xte6*Xte7
Xte6Xte8<-Xte6*Xte8

Xte7Xte8<-Xte7*Xte8

X.te<-cbind(X.te,Xte2Xte3,Xte2Xte4,Xte2Xte5,Xte2Xte6,Xte2Xte7,Xte2Xte8,Xte3Xte4,Xte3Xte5,Xte3Xte6,Xte3Xte7,Xte3Xte8,Xte4Xte5,Xte4Xte6,Xte4Xte7,Xte4Xte8,Xte5Xte6,Xte5Xte7,Xte5Xte8,Xte6Xte7,Xte6Xte8,Xte7Xte8)
X.te<-as.matrix(X.te)
#sigma OLS
model<-lm(Y~-1+X)
summary(model)$sigma^2

p<-dim(X)[2] 
S<-1000

BETA<-Z<-matrix(NA,S,p) 
S2<-NULL
z<-rep(1,dim(X)[2] ) 
lpy.c<-lpy.X(Y,X[,z==1,drop=FALSE]) 

for(s in 1:S)
{
  for(j in sample(2:p)) 
  {
    zp<-z ; zp[j]<-1-zp[j]
    lpy.p<-lpy.X(Y,X[,zp==1,drop=FALSE])
    r<- (lpy.p - lpy.c)*(-1)^(zp[j]==0)
    z[j]<-rbinom(1,1,1/(1+exp(-r)))
    if(z[j]==zp[j]) {lpy.c<-lpy.p}
  }
  
  beta<-z
  if(sum(z)>0){
    temp<-lm.gprior(Y,X[,z==1,drop=FALSE],S=1)
    beta[z==1]<-temp$beta
    s2<-temp$s2}
  
  Z[s,]<-z
  BETA[s,]<-beta
  S2<-c(S2,s2)
  print(s)
}
#8:55 to 10:46
#9:17 to 10:03

#Generate posterior predictive values and store them
y.pred<-NULL 
for (s in 1:S){
  beta<-BETA[s,]
  s2<-S2[s]
  y.pred<-rbind(y.pred,c(X%*%beta)+rnorm(n,0,sqrt(s2)))
}

#Use prediction error as the assessing criteria to evaluate the model
beta.bma<-apply(BETA,2,mean)
y.te.bma<-X.te%*%beta.bma
mean( (Y.te-y.te.bma)^2)

#Generate the plot of posterior probabilities that each coefficient is non-zero
plot(apply(Z,2,mean,na.rm=TRUE),ylim=c(0,1),xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2)


#Residual plot diagnostics
e.pred<-NULL
for (i in 1:S){
  e.pred<-rbind(e.pred,Y-y.pred[s,])
}
e.pred.avg<-apply(e.pred,2,mean)
plot(e.pred.avg,pch=19)
abline(h=0)

##QQ plot diagnostics 
qqnorm(e.pred.avg)
qqline(e.pred.avg)


##Create acf plot 
par(mfrow=c(4,4))
t<-seq(1,30,1)
for (i in 1:30){
  xname<-bquote(beta[.(t[i])])
  acf(BETA[,i],xlab=xname,cex.axis=0.8)
}
acf(S2,xlab=expression(sigma^2),cex.axis=0.8)

#effective sizes
library(coda)
apply(BETA,2,function(x) effectiveSize(x))
effectiveSize(S2)

##Create stationarity plot
par(mfrow=c(4,4))
for (i in 1:30){
  stationarity.plot(BETA[,i],xlab="iteration",ylab=expression(beta[i]),xaxt="n")
}
stationarity.plot(S2,xlab="iteration",ylab=expression(sigma^2),xaxt="n")

##Thinning process
thin<-thin<-c(1,(1:50)*(S/50))

par(mfrow=c(4,4))
##Create acf plot after thinning
t<-seq(1,30,1)
for (i in 1:30){
  xname<-bquote(beta[.(t[i])])
  acf(BETA[thin,i],xlab=xname,cex.axis=0.8, ylim=c(0,1))
}
acf(S2[thin],xlab=expression(sigma^2),cex.axis=0.8)


##effective sizes after thinning
library(coda)
apply(BETA[thin,],2,function(x) effectiveSize(x))
effectiveSize(S2[thin])

##Create stationarity plot after thinning
par(mfrow=c(4,4))
for (i in 1:30){
  stationarity.plot(BETA[thin,i],xlab="iteration",ylab=expression(beta[i]),xaxt="n")
}
stationarity.plot(S2[thin],xlab="iteration",ylab=expression(sigma^2),xaxt="n")

Y<-train[,9]
X<-train[,-9]
n<-length(Y)
X<-cbind(rep(1,n),X)
X<-as.matrix(X)
colnames(X)<-c("Intercept","X1","X2","X3","X4","X5","X6","X7","X8")

Y.te<-test[,9]
X.te<-test[,-9]
X.te<-cbind(rep(1,length(Y.te)),X.te)
X.te<-as.matrix(X.te)
colnames(X.te)<-c("Intercept","X1","X2","X3","X4","X5","X6","X7","X8")

#sigma OLS
model<-lm(Y~-1+X)
summary(model)$sigma^2

p<-dim(X)[2] 
S<-5000

BETA<-Z<-matrix(NA,S,p) 
S2<-NULL
z<-rep(1,dim(X)[2] ) 
lpy.c<-lpy.X(Y,X[,z==1,drop=FALSE]) 

for(s in 1:S)
{
  for(j in sample(2:p)) 
  {
    zp<-z ; zp[j]<-1-zp[j]
    lpy.p<-lpy.X(Y,X[,zp==1,drop=FALSE])
    r<- (lpy.p - lpy.c)*(-1)^(zp[j]==0)
    z[j]<-rbinom(1,1,1/(1+exp(-r)))
    if(z[j]==zp[j]) {lpy.c<-lpy.p}
  }
  
  beta<-z
  if(sum(z)>0){
    temp<-lm.gprior(Y,X[,z==1,drop=FALSE],S=1)
    beta[z==1]<-temp$beta
    s2<-temp$s2}
  
  Z[s,]<-z
  BETA[s,]<-beta
  S2<-c(S2,s2)
  print(s)
}

#Generate posterior predictive values and store them
y.pred<-NULL 
for (s in 1:S){
  beta<-BETA[s,]
  s2<-S2[s]
  y.pred<-rbind(y.pred,c(X%*%beta)+rnorm(n,0,sqrt(s2)))
}

#Use prediction error as the assessing criteria to evaluate the model
beta.bma<-apply(BETA,2,mean)
y.te.bma<-X.te%*%beta.bma
mean( (Y.te-y.te.bma)^2)

par(mfrow=c(1,1))
#Generate the plot of posterior probabilities that each coefficient is non-zero
plot(apply(Z,2,mean,na.rm=TRUE),ylim=c(0,1),xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2)


#Residual plot diagnostics
e.pred<-NULL
for (i in 1:S){
  e.pred<-rbind(e.pred,Y-y.pred[s,])
}
e.pred.avg<-apply(e.pred,2,mean)
plot(e.pred.avg,pch=19)
abline(h=0)

##QQ plot diagnostics 
qqnorm(e.pred.avg)
qqline(e.pred.avg)


##Create acf plot 
par(mfrow=c(3,4),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
t<-seq(1,9,1)
for (i in 1:9){
  xname<-bquote(beta[.(t[i])])
  acf(BETA[,i],xlab=xname,cex.axis=0.8)
}
acf(S2,xlab=expression(sigma^2),cex.axis=0.8)

#effective sizes
library(coda)
apply(BETA,2,function(x) effectiveSize(x))
effectiveSize(S2)

##Create stationarity plot
par(mfrow=c(3,4),mar=c(2.75,2.75,.5,.5),mgp=c(1.7,.7,0))
for (i in 1:9){
  stationarity.plot(BETA[,i],xlab="iteration",ylab=expression(beta[i]),xaxt="n")
}
stationarity.plot(S2,xlab="iteration",ylab=expression(sigma^2),xaxt="n")

##Thinning process
thin<-thin<-c(1,(1:100)*(S/100))

par(mfrow=c(4,3))
##Create acf plot after thinning
t<-seq(1,9,1)
for (i in 1:9){
  xname<-bquote(beta[.(t[i])])
  acf(BETA[thin,i],xlab=xname,cex.axis=0.8, ylim=c(0,1))
}
acf(S2[thin],xlab=expression(sigma^2),cex.axis=0.8)


##effective sizes after thinning
library(coda)
apply(BETA[thin,],2,function(x) effectiveSize(x))
effectiveSize(S2[thin])

##Create stationarity plot after thinning
par(mfrow=c(4,3))
for (i in 1:9){
  stationarity.plot(BETA[thin,i],xlab="iteration",ylab=expression(beta[i]),xaxt="n")
}
stationarity.plot(S2[thin],xlab="iteration",ylab=expression(sigma^2),xaxt="n")


y.pred.avg<-apply(y.pred,1,mean)
mean(y.pred.avg > mean(Y))

y.pred.median<-apply(y.pred,1,median)
mean(y.pred.median > median(Y))


library(e1071)
test1<-mean(Y) #mean
test2<-quantile(Y,0.5) #median


multi.fun<-function(x){
  rr<-c(pred.test1<-mean(x),pred.test2<-median(x))
  names(rr)<-c("mean","median")
  return(rr)
}

pred.test<-apply(y.pred, 1, multi.fun)

mean(pred.test[1,]>test1)
1-mean(pred.test[2,]>test2)



par(mfrow=c(1,2))
hist(pred.test[1,],xlab=expression(bar(y)),main="",breaks=20)
abline(v=test1,col="red")
hist(pred.test[2,],xlab="median(y)",main="",breaks=20)
abline(v=test2,col="red")





