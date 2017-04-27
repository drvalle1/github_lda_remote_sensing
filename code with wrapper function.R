rm(list=ls(all=T))
library('MCMCpack')
library(Rcpp)
set.seed(1)

setwd('U:\\pedro\\lda remote sensing\\github_lda_remote_sensing')
sourceCpp('aux1.cpp')
source('gibbs functions.R')
source('gibbs.R')

#each row corresponds to a different pixel and 
#each column represents the reflectance in a particular band
remote=data.matrix(read.csv('fake data remote.csv',as.is=T))

#input
#nloc: number of locations (i.e., pixels)
#nbands: number of remote sensing bands
#ndig.values: maximum pixel value (8 bit = 256)
#ngibbs: number of gibbs iterations
#ncommun: maximum number of communities/endmembers
#gamma: prior for truncated stick-breaking prior (smaller = more parsimonius = less endmembers)
#a.omega, b.omega: beta distribution hyper-prior parameters for omega matrix (set to 1 implies an uniform distribution)

nloc=nrow(remote)
nbands=ncol(remote)
ncommun=10
ngibbs=10000
  
results=gibbs.RS(nloc=nloc,nbands=nbands,
                 ndig.values=256,
                 ngibbs=ngibbs,
                 ncommun=ncommun,
                 gamma=0.01,a.omega=1,b.omega=1)

#output: 2 matrices (theta and omega) and each row represents a sample from the posterior distribution

#remove burn-in period and calculate average
nburnin=5000
seq1=nburnin:ngibbs
theta=matrix(colMeans(results$theta[seq1,]),nloc,ncommun)
omega=matrix(colMeans(results$omega[seq1,]),ncommun,nbands)

#theta: has nloc rows and ncommun columns. 
#This matrix contains the proportion of each endmember in each
#location/pixel. Rows should sum to 1

#omega: has ncommun rows and nbands columns. 
#Each row of this matrix contains the spectral signature of a endmember 
#represented as probabilities. Rows DO NOT have to sum to 1
#Expected reflectance is given by: ndig.values*omega

#------------------------------
#It seems that only the 5 first groups are important
boxplot(theta)
plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:5) lines(theta[,i],col=i)

