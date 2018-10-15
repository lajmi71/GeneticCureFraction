

####################################################################
## Compute D, P, R, Lambda0 and covz under the null hypothesis    ##
## Data is a data.frame and not a list                            ##
####################################################################

Prep.Score=function(Data)
{
fit=smcure(Surv(T,delta)~Z1+Z2,cureform=~X1+X2,model="ph",Var=FALSE,data=Data)

Z1=Data$Z1
Z2=Data$Z2
X1=Data$X1
X2=Data$X2

T=Data$T
delta=Data$delta
C=Data$C
S=fit$s

S=pmax(S,min(S[S>0]))

Lambda=-log(S)

Z=cbind(Z1,Z2)
X=cbind(rep(1),X1,X2)

mu.hat=fit$b
alpha.hat=fit$beta

covx=as.vector(exp(X%*%mu.hat))
P=covx/(1+covx)

covz=as.vector(exp(Z%*%alpha.hat))

A=P*(S^covz)
B=1-P

D=delta+(1-delta)*A/(A+B)

R=D*(delta-Lambda*covz)

rr = order(T)
sorted.T = sort(T)
sorted.Lambda = c(0,Lambda[rr])
Lambda0 = stepfun(sorted.T,sorted.Lambda)
Lambda0.C = Lambda0(C)

return(list(D=D,P=P,R=R,Lambda0.C=Lambda0.C,covz=covz))
}

################################################################################

Compute.Var.D = function(Prep)
{
P = Prep$P
Lambda0.C = Prep$Lambda0.C
covz = Prep$covz

term1 = 0 

for(i in 1:n)
{
covz.i = covz[i]
A = exp(-Lambda0.C*covz.i)
term1[i]=mean(A)
}

term2 = 0

for(i in 1:n)
{
covz.i = covz[i]
P.i = P[i]
A = exp(-2*Lambda0.C*covz.i)
B = exp(-Lambda0.C*covz.i)
C = (P.i*A)/(1-P.i+P.i*B)
term2[i] = mean(C)
}

return(P*(1-term1+term2)-P*P)

}

Compute.Var.R = function(Prep)
{
P = Prep$P
Lambda0.C = Prep$Lambda0.C
covz = Prep$covz

term1 = P

term2 = 0
for(i in 1:n)
{
covz.i = covz[i]
A = exp(-Lambda0.C*covz.i)
term2[i] = 1-mean(A)
}

term3 = P*(1-P)

term4 = 0
for(i in 1:n)
{
P.i = P[i]
covz.i = covz[i]
A = (Lambda0.C*covz.i)^2*exp(-Lambda0.C*covz.i)
B = 1 - P.i + P.i*exp(-Lambda0.C*covz.i)
term4[i] = mean(A/B)
}

return(term1*term2-term3*term4)
}

Compute.Cov.D.R =  function(Prep)
{
P = Prep$P
Lambda0.C = Prep$Lambda0.C
covz = Prep$covz

term1 = P*(1-P)

term2 = 0
for(i in 1:n)
{
P.i = P[i]
covz.i = covz[i]
A = (Lambda0.C*covz.i)*exp(-Lambda0.C*covz.i)
B = 1 - P.i + P.i*exp(-Lambda0.C*covz.i)
term2[i] = mean(A/B)
}

return(term1*term2)
}



Compute.test.M1 = function(data,phi,rho)
{

Prep = Prep.Score(data.frame(T=data$T,delta=data$delta,C=data$C,X1=data$X1,X2=data$X2,Z1=data$Z1,Z2=data$Z2))
vect = c(Prep$D-Prep$P,Prep$R)

G=data$G
freq.MAF=apply(G,2,mean)/2
w = (dbeta(freq.MAF,1,25))^2
W = diag(w)
K = G%*%W%*%t(G)
tuning.params.matrix = matrix(c(phi,rho*sqrt(phi*(1-phi)),rho*sqrt(phi*(1-phi)),1-phi),ncol=2)
K.tilde = kronecker(tuning.params.matrix,K)
test.stat = (t(vect)%*%K.tilde%*%vect)[1,1]

Var.D = Compute.Var.D(Prep)
Var.R = Compute.Var.R(Prep)
Cov.D.R = Compute.Cov.D.R(Prep)

V = rbind(cbind(diag(Var.D),diag(Cov.D.R)),cbind(diag(Cov.D.R),diag(Var.R)))


V.sqrt = Compute.V.sqrt(V)

mat = V.sqrt%*%K.tilde%*%V.sqrt
ev = eigen(mat)
theta = ev$values[ev$values>0]

p.value = davies(test.stat,theta)$Qq

return(data.frame(test.stat=test.stat,p.value=p.value))
}

###########################################################################################

###############################################
## Compute the quadratic forms Q1, Q2 and Q3 ##
###############################################

Compute.Q=function(Prep,K)
{
A=t(Prep$D-Prep$P)%*%K
Q1=A%*%(Prep$D-Prep$P)
Q2=t(Prep$R)%*%K%*%(Prep$R)
Q3=A%*%(Prep$R)
return(c(Q1,Q2,Q3))
}

##############################################################
## Compute the permuted quadratic forms Q1*b, Q2*b and Q3*b ##
############################################################## 

Prep.Perm=function(Prep,ind,K)
{
R=Prep$R
D=Prep$D
P=Prep$P

A1=matrix((D-P)[ind],ncol=n)
A2=matrix(R[ind],ncol=n)
B1=A1%*%K
B2=A2%*%K
Q1=sapply(1:B,mult,A=A1,B=B1)
Q2=sapply(1:B,mult,A=A2,B=B2)
Q3=sapply(1:B,mult,A=A2,B=B1)

return(cbind(Q1,Q2,Q3))
}

Compute.test.M2.M3=function(data,B)
{

Prep = Prep.Score(data.frame(T=data$T,delta=data$delta,C=data$C,X1=data$X1,X2=data$X2,Z1=data$Z1,Z2=data$Z2))
G=data$G
freq.MAF=apply(G,2,mean)/2
w = (dbeta(freq.MAF,1,25))^2
W = diag(w)
K = G%*%W%*%t(G)

n = dim(G)[1]

ind = create.indices(n,B)

Q.obs=Compute.Q(Prep,K)
Q.perm=Prep.Perm(Prep,ind,K)

Mat = Compute.Grid()

score.obs=as.vector(Q.obs%*%t(Mat))
score.perm=Q.perm%*%t(Mat)

len=dim(Mat)[1]
p.value.obs=sapply(1:len,rank.vect.mat,vect=score.obs,matr=score.perm)

p.value.perm=1-sapply(1:len,rank.mat,M=score.perm)/(B+1)

p.min.obs=min(p.value.obs)
p.min.perm=apply(p.value.perm,1,min)
p.value.M2=mean(p.min.obs>p.min.perm)

p.log.obs=-log(p.min.obs)
p.log.perm=-log(p.min.perm)

m=mean(p.log.perm)
v=mean((p.log.perm-m)^2)
c=mean((p.log.perm-m)^4)
gam=(c/v^2)-3
df=12/gam
A=((p.log.obs-m)*sqrt(2*df)/sqrt(v))+df
p.value.M3=1-pchisq(A,df=df)

return(data.frame(test.stat=p.min.obs,p.value.M2=p.value.M2,p.value.M3=p.value.M3))
}


