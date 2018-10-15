library(CompQuadForm)
library(smcure)
library(survival)

Compute.V.sqrt = function(V)
{
d = dim(V)[1]
n = d/2

V.sqrt = matrix(0,nrow=d,ncol=d)

for (i in 1:n)
{
A = V[i,i]
B = V[i,(i+n)]
C = B
D = V[(i+n),(i+n)]
del = A*D-B*C
s = sqrt(del)
tau = A+D
t = sqrt(tau+(2*s))
V.sqrt[i,i]=(A+s)/t
V.sqrt[i,(i+n)]=B/t
V.sqrt[(i+n),i]=C/t
V.sqrt[(i+n),(i+n)]=(D+s)/t
} 

return(V.sqrt)
}





#######################################################################
## Return a B by n matrix of permuted indices                        ##
## The bth line of the returned matrix is a permutation of {1,...,n} ##
####################################################################### 


create.indices.i = function(i,n) {sample(n,n)}

create.indices = function(n,B) {t(sapply(1:B,create.indices.i,n=n))}


#########################################################
## Compute the ith element of the diagonal of t(A)%*%B ##
#########################################################

mult=function(i,A,B) {return(t(A[i,])%*%B[i,])}

Compute.Grid=function()
{
rho.vect=seq(-0.9,0.9,0.1)
phi.vect=seq(0.1,0.9,0.1)
l.rho=length(rho.vect)
l.phi=length(phi.vect)
Mat=NULL
for(i in 1:l.rho)
{
for(j in 1:l.phi)
{
v=c(phi.vect[j],1-phi.vect[j],2*rho.vect[i]*sqrt(phi.vect[j]*(1-phi.vect[j])))
Mat=rbind(Mat,v)
}
}
return(Mat)
}

rank.mat=function(i,M) {rank(M[,i])}

rank.vect.mat = function(i,vect,matr) {mean(vect[i]<matr[,i])}