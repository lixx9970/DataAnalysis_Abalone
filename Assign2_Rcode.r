#Source code from Ch 7 - Hoff
###
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
###



###
ldmvnorm<-function(y,mu,Sig){  # log mvn density
  c(  -(length(mu)/2)*log(2*pi) -.5*log(det(Sig)) -.5*
        t(y-mu)%*%solve(Sig)%*%(y-mu)   )  
}
####

### sample from the Wishart distribution
rwish<-function(n,nu0,S0)
{
  sS0 <- chol(S0)
  S<-array( dim=c( dim(S0),n ) )
  for(i in 1:n)
  {
    Z <- matrix(rnorm(nu0 * dim(S0)[1]), nu0, dim(S0)[1]) %*% sS0
    S[,,i]<- t(Z)%*%Z
  }
  S[,,1:n]
}

#### From the `sbgcop` package
plotci.sA<-function(sA, ylabs = colnames(sA[, , 1]), mgp = c(1.75, 0.75, 
                                                             0)) 
{
  qA <- qM.sM(sA)
  p <- dim(qA)[1]
  tmp <- c(qA)
  tmp <- tmp[tmp != 1]
  par(mgp = mgp)
  for (j in 1:p) {
    plot(0, 0, type = "n", ylim = range(c(tmp), na.rm = TRUE), 
         xlim = c(1, p), ylab = ylabs[j], xaxt = "n", xlab = "")
    points((1:p)[-j], qA[j, -j, 2], pch = 16, cex = 0.6)
    segments(x0 = (1:p)[-j], y0 = qA[j, -j, 1], x1 = (1:p)[-j], 
             y1 = qA[j, -j, 3])
    abline(h = 0, col = "gray")
    abline(v = j, col = "gray")
  }
  axis(side = 1, at = 1:p, labels = colnames(qA[, , 1]), las = 2,cex.axis=0.6)
}

sR.sC<-function(sC) 
{
  p <- dim(sC)[1]
  s <- dim(sC)[3]
  sR <- array(dim = c(p, p, s))
  dimnames(sR) <- dimnames(sC)
  for (l in 1:s) {
    C <- sC[, , l]
    R <- C * NA
    for (j in 1:p) {
      R[j, -j] <- C[j, -j] %*% solve(C[-j, -j])
    }
    sR[, , l] <- R
  }
  sR
}

#### A plot for evaluating lack of convergence
stationarity.plot<-function(x,...){
  
  S<-length(x)
  scan<-1:S
  ng<-min( round(S/100),10)
  group<-S*ceiling( ng*scan/S) /ng
  
  boxplot(x~group,...)               }




#Problem 1
#a)
tennis<-read.csv("~/Desktop/tennis.csv",header=TRUE)

Grp1<-tennis[1:50,]
Grp2<-tennis[51:100,]

n1<-nrow(Grp1)
n2<-nrow(Grp2)

#Group 1 (top50) multivariate analysis
Y<-Grp1
n<-n1
ybar<-apply(Y,2,mean) #sample means for each variable
Sigma<-cov(Y) #sample covariance matrix
THETA_Grp1<-SIGMA_Grp1<-NULL #objects to store posterior draws in

##Prior parameters
mu0<-ybar
L0<-S0<-Sigma
nu0<-dim(Y)[2]+2

set.seed(1)

for (s in 1:10000){
  ###update theta
  Ln<-solve(solve(L0)+n*solve(Sigma))
  mun<-Ln%*%(solve(L0)%*%mu0+n*solve(Sigma)%*%ybar)
  theta<-rmvnorm(1,mun,Ln)
  
  ##update Sigma
  Sn<-S0+(t(Y)-c(theta))%*%t(t(Y)-c(theta))
  Sigma<-solve(rwish(1,nu0+n,solve(Sn)))
  
  ##save results
  THETA_Grp1<-rbind(THETA_Grp1,theta)
  SIGMA_Grp1<-rbind(SIGMA_Grp1,c(Sigma))
}

#Group 2 (rank 51-100) multivariate analysis
Y<-Grp2
n<-n2
ybar<-apply(Y,2,mean) #sample means for each variable
Sigma<-cov(Y) #sample covariance matrix
THETA_Grp2<-SIGMA_Grp2<-NULL #objects to store posterior draws in

##Prior parameters
mu0<-ybar
L0<-S0<-Sigma
nu0<-dim(Y)[2]+2

set.seed(1)

for (s in 1:10000){
  ###update theta
  Ln<-solve(solve(L0)+n*solve(Sigma))
  mun<-Ln%*%(solve(L0)%*%mu0+n*solve(Sigma)%*%ybar)
  theta<-rmvnorm(1,mun,Ln)
  
  ##update Sigma
  Sn<-S0+(t(Y)-c(theta))%*%t(t(Y)-c(theta))
  Sigma<-solve(rwish(1,nu0+n,solve(Sn)))
  
  ##save results
  THETA_Grp2<-rbind(THETA_Grp2,theta)
  SIGMA_Grp2<-rbind(SIGMA_Grp2,c(Sigma))
}

#b)
p<-dim(Y)[2]
#Posterior mean estimate of correlation matrix
COR_Grp1 <- array( dim=c(p,p,10000) )
for(s in 1:10000)
{
  Sig<-matrix( SIGMA_Grp1[s,] ,nrow=p,ncol=p)
  COR_Grp1[,,s] <- Sig/sqrt( outer( diag(Sig),diag(Sig) ) )
}

apply(COR_Grp1,c(1,2),mean)
colnames(COR_Grp1)<-rownames(COR_Grp1)<-c("1stServe%","1stServePts%","2ndServePts%","BPtsSaved%","ServGamesWon%")

COR_Grp2 <- array( dim=c(p,p,10000) )
for(s in 1:10000)
{
  Sig<-matrix( SIGMA_Grp2[s,] ,nrow=p,ncol=p)
  COR_Grp2[,,s] <- Sig/sqrt( outer( diag(Sig),diag(Sig) ) )
}

apply(COR_Grp2,c(1,2),mean)

colnames(COR_Grp2)<-rownames(COR_Grp2)<-c("1stServe%","1stServePts%","2ndServePts%","BPtsSaved%","ServGamesWon%")

#plots of 95% posterior confidence intervals of correlations
par(mfcol=c(5,2),mar=c(1,2.75,1,1),mgp=c(1.75,.75,0),oma=c(1.5,0,0,0))
plotci.sA(COR_Grp1)

REG_Grp1<-sR.sC(COR_Grp1)
plotci.sA(REG_Grp1)

par(mfcol=c(5,2),mar=c(1,2.75,1,1),mgp=c(1.75,.75,0),oma=c(1.5,0,0,0))
plotci.sA(COR_Grp2)

REG_Grp2<-sR.sC(COR_Grp2)
plotci.sA(REG_Grp2)

#c) 
par(mar=c(3,3,.25,1),mgp=c(1.75,.75,0))
par(mfrow=c(3,2))
for (i in 1:5){
  xname<-bquote(theta[.(names(Y)[i])])
  xl<-c(min(range(THETA_Grp1[,i])[1],range(THETA_Grp2[,i])[1]),max(range(THETA_Grp1[,i])[2],range(THETA_Grp2[,i])[2]))
  d_Grp1<-density(THETA_Grp1[,i])
  d_Grp2<-density(THETA_Grp2[,i])
  yl<-c(min(min(d_Grp1[2]$y),min(d_Grp2[2]$y)),max(max(d_Grp1[2]$y),max(d_Grp2[2]$y)))
  plot(d_Grp1,xlab=xname,ylab="",main="",xlim=xl,ylim=yl,col="red",cex.axis=0.8)
  lines(d_Grp2)
  legend(x="topleft",c("Grp1","Grp2"),col=c("red","black"),lty=c(1,1),cex=0.5)
}

#Compute marginal posterior probabilities 
pp<-NULL
for (i in 1:5){
  pp<-c(pp,mean(THETA_Grp1[,i]>THETA_Grp2[,i]))
  
}
names(pp)<-names(Y)
pp

#d) convergence diagnostics

#effective sizes
library(coda)
apply(THETA_Grp1, 2, function(x) effectiveSize(x))
apply(THETA_Grp2, 2, function(x) effectiveSize(x))
apply(SIGMA_Grp1, 2, function(x) effectiveSize(x))
apply(SIGMA_Grp2, 2, function(x) effectiveSize(x))

#autocorrelation plots
par(mar=c(3,3,.25,1),mgp=c(1.75,.75,0))
par(mfrow=c(3,2))
for (i in 1:5){
  xname<-bquote(theta[.(names(Y)[i])])
  acf(THETA_Grp1[,i],xlab=xname,cex.axis=0.8)
}

par(mfrow=c(3,2))
for (i in 1:5){
  xname<-bquote(theta[.(names(Y)[i])])
  acf(THETA_Grp2[,i],xlab=xname,cex.axis=0.8)
}

par(mfrow=c(5,5))
t<-c(11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45,51,52,53,54,55)
for (i in 1:25){
  xname<-bquote(Sigma[.(t[i])])
  acf(SIGMA_Grp1[,i],xlab=xname)
}

par(mfrow=c(5,5))
t<-c(11,12,13,14,15,21,22,23,24,25,31,32,33,34,35,41,42,43,44,45,51,52,53,54,55)
for (i in 1:25){
  xname<-bquote(Sigma[.(t[i])])
  acf(SIGMA_Grp2[,i],xlab=xname)
}

#stationarity plots
par(mar=c(3,3,.25,1),mgp=c(1.75,.75,0))
par(mfrow=c(3,2))
for (i in 1:5){
  yname<-bquote(theta[.(names(Y)[i])])
 stationarity.plot(THETA_Grp1[,i],xlab="iteration",ylab=expression(theta))
}

par(mfrow=c(3,2))
for (i in 1:5){
  yname<-bquote(theta[.(names(Y)[i])])
  stationarity.plot(THETA_Grp2[,i],xlab="iteration",ylab=expression(theta))
}

par(mfrow=c(5,5))
for (i in 1:25){
  stationarity.plot(SIGMA_Grp1[,i],xlab="iteration",ylab=expression(Sigma))
}

par(mfrow=c(5,5))
for (i in 1:25){
  stationarity.plot(SIGMA_Grp2[,i],xlab="iteration",ylab=expression(Sigma))
}


#Problem 2
source("~/Desktop/regression_gprior.R")

Y<-tennis[1:100,5]
X<-tennis[1:100,1:4]
n<-length(Y)
X<-cbind(c(rep(1,n)),X)
colnames(X)<-c("Intercept","X1","X2","X3","X4")
#main effects
X1<-X[,2]
X2<-X[,3]
X3<-X[,4]
X4<-X[,5]
#two-way interactions
X1X2<-X1*X2
X1X3<-X1*X3
X1X4<-X1*X4
X2X3<-X2*X3
X2X4<-X2*X4
X3X4<-X3*X4
X1sq<-X1^2
X2sq<-X2^2
X3sq<-X3^2
X4sq<-X4^2


X<-cbind(X,X1X2,X1X3,X1X4,X2X3,X2X4,X3X4)
X<-as.matrix(X)

#sigma OLS
model<-lm(Y~-1+X)
summary(model)$sigma^2

p<-dim(X)[2] # number of regression coefficients including intercept
S<-20000

BETA<-Z<-matrix(NA,S,p) #object to store posterior draws in 
S2<-NULL
z<-rep(1,dim(X)[2] ) #starting values of z
lpy.c<-lpy.X(Y,X[,z==1,drop=FALSE]) #starting marginal log likelihood of data

for(s in 1:S)
{
  for(j in sample(2:p)) #always include intercept term
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
}

#c)

par(mfrow=c(1,1))
plot(apply(Z,2,mean,na.rm=TRUE),ylim=c(0,1),xlab="regressor index",ylab=expression(
  paste( "Pr(",italic(z[j] == 1),"|",italic(y),",X)",sep="")),type="h",lwd=2)

apply(BETA,2,mean)

#d)
y.pred<-NULL #object to store posterior predictive values
for (s in 1:S){
  beta<-BETA[s,]
  s2<-S2[s]
  y.pred<-rbind(y.pred,c(X%*%beta)+rnorm(n,0,sqrt(s2)))
  
}
#Residual plot diagnostics
e.pred<-NULL
for (i in 1:S){
  e.pred<-rbind(e.pred,Y-y.pred[s,])
}

e.pred.avg<-apply(e.pred,2,mean)
par(mfrow=c(1,1))
plot(e.pred.avg,pch=19)
abline(h=0)

#e)

#ACF plots
par(mfrow=c(3,4))
t<-seq(1,11,1)
for (i in 1:11){
  xname<-bquote(beta[.(t[i])])
  acf(BETA[,i],xlab=xname,cex.axis=0.8)
}

thin<-thin<-c(1,(1:1000)*(S/1000))

par(mfrow=c(3,4))
t<-seq(1,11,1)
for (i in 1:11){
  xname<-bquote(beta[.(t[i])])
  acf(BETA[thin,i],xlab=xname,cex.axis=0.8)
}
acf(S2[thin],xlab=expression(sigma^2),cex.axis=0.8)

#effective sizes
apply(BETA[thin,],2,function(x) effectiveSize(x))
effectiveSize(S2[thin])

#Stationarity plots (not required)

par(mfrow=c(3,4))

for (i in 1:11){
  stationarity.plot(BETA[thin,i],xlab="iteration",ylab=expression(beta[i]),xaxt="n")
}
stationarity.plot(S2[thin],xlab="iteration",ylab=expression(sigma^2),xaxt="n")


#Problem 3
athlete<-read.csv("athlete.csv",header=TRUE)
athlete$y<-athlete$Post-athlete$Pre
#subset data by Group
IHE<-subset(athlete,Group=="IHE")
Heat<-subset(athlete,Group=="Heat")
Placebo<-subset(athlete,Group=="Placebo")

m<-3 #number of groups

#list object to store ferritin loss by group
Y<-list()
Y[[1]]<-IHE$y
Y[[2]]<-Heat$y
Y[[3]]<-Placebo$y

#summary statistics by group
y1<-mean(IHE$y)
y2<-mean(Heat$y)
y3<-mean(Placebo$y)

sv1<-var(IHE$y)
sv2<-var(Heat$y)
sv3<-var(Placebo$y)

n1<-length(IHE$y)
n2<-length(Heat$y)
n3<-length(Placebo$y)

ybar<-c(y1,y2,y3)
sv<-c(sv1,sv2,sv3)
n<-c(n1,n2,n3)


## weakly informative priors
nu0<-1  ; s20<-1
eta0<-1 ; t20<-1
mu0<-0 ; g20<-1000

#starting values
theta<-ybar
sigma2<-mean(sv)
mu<-mean(theta)
tau2<-var(theta)

## setup MCMC
set.seed(1)
S<-100000
THETA<-matrix( nrow=S,ncol=m) #object to store posterior draws of theta_j in
MST<-matrix( nrow=S,ncol=3) #object to store posterior draws of sigma2, mu, tau

## MCMC algorithm
for(s in 1:S) 
{
  
  # sample new values of the thetas
  for(j in 1:m) 
  {
    vtheta<-1/(n[j]/sigma2+1/tau2)
    etheta<-vtheta*(ybar[j]*n[j]/sigma2+mu/tau2)
    theta[j]<-rnorm(1,etheta,sqrt(vtheta))
  }
  
  #sample new value of sigma2
  nun<-nu0+sum(n)
  ss<-nu0*s20;for(j in 1:m){ss<-ss+sum((Y[[j]]-theta[j])^2)}
  sigma2<-1/rgamma(1,nun/2,ss/2)
  
  #sample a new value of mu
  vmu<- 1/(m/tau2+1/g20)
  emu<- vmu*(m*mean(theta)/tau2 + mu0/g20)
  mu<-rnorm(1,emu,sqrt(vmu)) 
  
  # sample a new value of tau2
  etam<-eta0+m
  ss<- eta0*t20 + sum( (theta-mu)^2 )
  tau2<-1/rgamma(1,etam/2,ss/2)
  
  #store results
  THETA[s,]<-theta
  MST[s,]<-c(mu,sigma2,tau2)
} 

#c) diagnostic plots


#ACF
par(mfrow=c(2,3))
acf(THETA[,1],xlab=expression(theta[IHE]))
acf(THETA[,2],xlab=expression(theta[Heat]))
acf(THETA[,3],xlab=expression(theta[Placebo]))
acf(MST[,1],xlab=expression(mu))
acf(MST[,2],xlab=expression(sigma^2))
acf(MST[,3],xlab=expression(tau^2))

#thin - take every 50th iteration
thin<-c(1,(1:2000)*(S/2000))

par(mfrow=c(2,3))
acf(THETA[thin,1],xlab=expression(theta[IHE]))
acf(THETA[thin,2],xlab=expression(theta[Heat]))
acf(THETA[thin,3],xlab=expression(theta[Placebo]))
acf(MST[thin,1],xlab=expression(mu))
acf(MST[thin,2],xlab=expression(sigma^2))
acf(MST[thin,3],xlab=expression(tau^2))


#convergence
par(mfrow=c(2,3))
stationarity.plot(THETA[thin,1],xlab="iteration",xaxt="n",ylab=expression(theta[IHE]))
stationarity.plot(THETA[thin,2],xlab="iteration",xaxt="n",ylab=expression(theta[Heat]))
stationarity.plot(THETA[thin,3],xlab="iteration",xaxt="n",ylab=expression(theta[Placebo]))
stationarity.plot(MST[thin,1],xlab="iteration",xaxt="n",ylab=expression(mu))
stationarity.plot(MST[thin,2],xlab="iteration",xaxt="n",ylab=expression(sigma^2))
stationarity.plot(MST[thin,3],xlab="iteration",xaxt="n",ylab=expression(tau^2))

#effective size
apply(MST[thin,],2,function(x) effectiveSize(x))


#d)
min.prop<-max.prop<-rep(0,m)
for (i in 1:m){
  min.ind<-apply(THETA[thin,-i],1,min)
  max.ind<-apply(THETA[thin,-i],1,max)
  min.prop[i]<-mean(THETA[thin,i]<min.ind)
  max.prop[i]<-mean(THETA[thin,i]>max.ind)
}
min.prop
max.prop

apply(THETA[thin,],2,mean)

#e)
cred<-apply(THETA[thin,],2,function(x) quantile(x,c(0.05,0.95)))
cred
cred[2,]-cred[1,]

#f)
R<-MST[thin,3]/(MST[thin,2]+MST[thin,3])
par(mar=c(3,3,1,1),mgp=c(1.75,.75,0))
par(mfrow=c(1,1))
plot(density(R,adj=2),xlab="R",cex.main=0.8,main="Monte Carlo estimate of posterior density of R")
quantile(R,c(0.05,0.90))

