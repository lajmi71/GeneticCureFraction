library(PredictABEL)
library(mvtnorm)

###########################################################################
## Generate Genotypes for independent individuals                        ##
## n: number of individuals                                              ##
## freq.min and freq.max : minimum and maximum minor allele frequencies  ##
## n.var : number of variants                                            ##
###########################################################################

Generate.Genotype=function(n,freq.min,freq.max,n.var)
{
ORfreq=cbind(runif(n.var,min=1,max=1.5),rep(1,n.var),runif(n.var,freq.min,freq.max),rep(1,n.var))
popRisk=0.5
popSize=n
Data=simulatedDataset(ORfreq=ORfreq, poprisk=popRisk, popsize=popSize)
return(as.matrix(Data[,-(1:4)]))
}



est.freq=function(Y) {sum(Y==1)+2*sum(Y=2)}

###########################################################################
## Generate survival traits with a cure fraction                         ##
## n: number of individuals                                              ##
## mu : regression coefficients for the logistic part                    ##
## alpha : regression coefficients for the conditional cox model         ##
## freq.min and freq.max : minimum and maximum minor allele frequencies  ##
## tau, phi, rho : variance components under the alternative hypothesis  ##
## phi has to be between 0 and 1                                         ##
## rho has to be between -1 and 1                                        ##
## tau = 0 corresponds to the null hypothesis                            ##
###########################################################################

Sim.c=function(n,mu,alpha,freq.min,freq.max,tau,phi,rho,n.var)
{
G=Generate.Genotype(n,freq.min,freq.max,n.var)

freq.MAF=apply(G,2,est.freq)/n
w = (dbeta(freq.MAF,1,25))^2
W = diag(w)

V=rbind(cbind(tau*phi*W,rho*tau*sqrt(phi*(1-phi))*W),cbind(rho*tau*sqrt(phi*(1-phi))*W,tau*(1-phi)*W))
betagamma=as.vector(rmvnorm(n=1,mean=rep(0,2*n.var),sigma=V))
beta=betagamma[1:(n.var)]
gamma=betagamma[(n.var+1):(2*n.var)]

X1=runif(n,min=-0.5,max=0.5)
X2=as.numeric(runif(n)<0.2)
X=cbind(rep(1,n),X1,X2)
Z1=X1
Z2=X2
Z=cbind(Z1,Z2)

covx=exp(as.vector(X%*%mu+G%*%beta))
covz=exp(as.vector(Z%*%alpha+G%*%gamma))

P=covx/(covx+1)

U=runif(n)

Y=ifelse((U<P),qweibull((1-(U/P))^(1/covz),shape=2,scale=2,lower.tail=FALSE),Inf)
C=rweibull(n,shape=2,scale=2)

T=pmin(Y,C)
delta=as.numeric(Y<C)

return(list(T=T,delta=delta,C=C,X1=X1,X2=X2,Z1=Z1,Z2=Z2,G=G))
}

###########################################################################
## Generate survival traits with no cure fraction                        ##
## n: number of individuals                                              ##
## alpha : regression coefficients for the cox model                     ##
## freq.min and freq.max : minimum and maximum minor allele frequencies  ##
## tau, phi: variance components under the alternative hypothesis        ##
## phi has to be between 0 and 1                                         ##
## tau = 0 corresponds to the null hypothesis                            ##
###########################################################################

Sim.nc=function(n,alpha,freq.min,freq.max,tau,phi,n.var)
{
G=Generate.Genotype(n,freq.min,freq.max,n.var)

freq.MAF=apply(G,2,est.freq)/n
w = (dbeta(freq.MAF,1,25))^2
W = diag(w)

V=tau*(1-phi)*W
gamma=as.vector(rmvnorm(n=1,mean=rep(0,n.var),sigma=V))

Z1=runif(n,min=-0.5,max=0.5)
Z2=as.numeric(runif(n)<0.2)
Z=cbind(Z1,Z2)

covz=exp(as.vector(Z%*%alpha+G%*%gamma))

U=runif(n)

Y=qweibull((U)^(1/covz),shape=2,scale=2,lower.tail=FALSE)
C=rweibull(n,shape=2,scale=2)

T=pmin(Y,C)
delta=as.numeric(Y<C)

return(list(T=T,delta=delta,C=C,X1=Z1,X2=Z2,Z1=Z1,Z2=Z2,G=G))
}

