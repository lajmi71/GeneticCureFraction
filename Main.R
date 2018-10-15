source("Functions.R")
source("Sim.R")
source("Score.R")

n=1000   ### sample size
mu=c(-0.2,1,1) ### Regression coefficients for the logistic part
alpha=c(1,1) ### Regression coefficients for the survival part
freq.min=0.001 ### Minimum minor allele frequency
freq.max=0.05 ### Maximum minor allele frequency
tau=0
phi=0.5
rho=0.5
n.var=5 ### Number of genetic variants

B=1000 ### Number of permutations

data=Sim.c(n,mu,alpha,freq.min,freq.max,tau,phi,rho,n.var) #### Generate data 	
u=Compute.test.M1(data,phi,rho)  ### Compute test statistic and p.value for M1.
v=Compute.test.M2.M3(data,B) ### Compute test statistic and p.values for M2 and M3.
