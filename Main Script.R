#Fitting procedure for logloss data and Monte Carlo Simulation
#The Monte Carlo simulation assumes a poisson distribution for the # of claims and a mixed lognormal distribution 
#for the size of loss based on the fitting distribution.
#MaryAnn Case 2022-02-17

library(MASS)
library(lattice)
library(ExcelFunctionsR)

#Main procedure
#Step 1: Build -> Set Working Directory to location of logloss file.

logloss<-as.matrix(read.csv("logloss.csv",header=FALSE))
logloss
theta0<-c(5.0,7.8,1.1,1.9,0.4)
theta0
mixnormal = function(x,theta0) {
  part1= (1-theta0[5])*dnorm(x,theta0[1],theta0[3])   
  #part2= theta0[5]*dnorm(x,theta0[2],theta0[4])
  gam = part2/(part1+part2)
  denom1<- sum(1-gam)
  denom2<- sum(gam)
  mu1<-sum((1-gam)*x)/denom1
  mu2<-sum(gam*x)/denom2
  sig1<-sqrt(sum((1-gam)*((x-mu1)^2))/denom1)
  sig2<-sqrt(sum(gam*((x-mu2) *(x-mu2))/denom2))  
  p<-mean(gam)
  mixnormal=c(mu1,mu2,sig1,sig2,p)
  mixnormal
}
mixnormal(logloss,theta0)
for( I in 1:1000)
{theta0= mixnormal(logloss,theta0)
mixnormal(logloss,theta0)}
k=theta0
k
a<-rnorm((1-k[5])*50000,mean=k[1],sd=k[3])
b<-rnorm((k[5])*50000,mean=k[2],sd=k[4]) 
c<-append(a,b)
plot( density(logloss), xlim=range( c(logloss),c),main="",xlab="")
lines(density(c), col=3)
lines(density(c), col=3)

#Start of the Monte Carlo Simulation


#parameters used

#Assuming a poisson distribution

lamb <- 240 #expected # of claims

#Assuming a mixed lognormal distribution 


n <- 50 #Number of iterations = accident years

Unlimited <- rep(NA,n)
Lim_100K <-rep(NA,n)
Lim_1M <- rep(NA,n)
TotalCounts <- rep(NA,n)
Total <- 0
Total_Lim <-0
Total_1M <-0

for(i in 1:n)
{
  #intialize Totals back to 0.
  Total <- 0
  Total_Lim <-0
  Total_1M <-0
  counts <- rpois(1,lamb) #rpois(n,lambda)
  counts 
  
  for (j in 1:counts){
    loss <- rlnorm(1,k[1],k[3])*(1-k[5]) + rlnorm(1,k[2],k[4])*k[5]
    loss
    Total <- Total + loss
    Total_Lim <- Total_Lim + min(loss,100000)
    Total_1M <- Total_1M + min(loss, 1000000)
  }
  
  TotalCounts[i]<-counts
  Unlimited[i]<- Total
  Lim_100K[i] <- Total_Lim
  Lim_1M[i] <- Total_1M}

data <- data.frame(TotalCounts,Unlimited,Lim_100K, Lim_1M)
write.csv(data,"data.csv") #data will be saved to the working directory.


