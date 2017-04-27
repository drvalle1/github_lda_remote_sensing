fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  log(jold)-log(jnew) #add this to pnew
}
#----------------------------------------------------------------------------------------------
tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#-------------------------------
print.adapt = function(accept1z,jump1z,accept.output){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<100
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------------------------------
update.omega=function(param,jump,ncommun,nbands,ndig.values,a.omega,b.omega){
  omega.orig=omega.old=param$omega
  tmp=tnorm(nbands*ncommun,lo=0,hi=1,mu=omega.old,sig=jump)
  novos=matrix(tmp,ncommun,nbands)

  tmp=fix.MH(lo=0,hi=1,omega.old,novos,jump)
  ajuste=matrix(tmp,ncommun,nbands)

  prior.old=matrix(dbeta(omega.old,a.omega,b.omega,log=T),ncommun,nbands)
  prior.new=matrix(dbeta(novos,a.omega,b.omega,log=T),ncommun,nbands)
  
  for (i in 1:ncommun){
    omega.new=omega.old
    omega.new[i,]=novos[i,]
      
    prob.old=param$theta%*%omega.old
    llk.old=colSums(dbinom(remote,size=ndig.values,prob=prob.old,log=T))
      
    prob.new=param$theta%*%omega.new
    llk.new=colSums(dbinom(remote,size=ndig.values,prob=prob.new,log=T))
      
    k=acceptMH(llk.old+prior.old[i,],
               llk.new+prior.new[i,]+ajuste[i,],
              omega.old[i,],omega.new[i,],F)
    omega.old[i,]=k$x
  }
  list(omega=omega.old,accept=omega.old!=omega.orig)
}
#------------------------------------------------
# make.from.SB.mat=function(vmat){
#   n=ncol(vmat)
#   mat=matrix(NA,nrow(vmat),n)
#   mat[,1]=vmat[,1]
#   mat[,2]=vmat[,2]*(1-vmat[,1])
#   for (i in 3:n){
#     mat[,i]=vmat[,i]*apply(1-vmat[,1:(i-1)],1,prod)
#   }
#   mat
# }

# vmat=matrix(rnorm(nloc*ncommun),nloc,ncommun)
# res=make.from.SB.mat(param$v)
# res1=convertSBtoNormal(vmat=param$v,
#                        ncol=ncol(param$v),nrow=nrow(param$v),
#                        prod=rep(1,nloc))
# range(res1-res)
#------------------------------------------------
update.theta=function(param,jump,ncommun,nloc,ndig.values){
  v.orig=v.old=param$v
  tmp=tnorm(nloc*(ncommun-1),lo=0,hi=1,mu=v.old[,-ncommun],sig=jump[,-ncommun])
  novos=cbind(matrix(tmp,nloc,ncommun-1),1)
  ajuste=matrix(fix.MH(lo=0,hi=1,v.old,novos,jump),nloc,ncommun)
  
  prior.old=matrix(dbeta(v.old,1,param$gamma,log=T),nloc,ncommun)
  prior.new=matrix(dbeta(novos,1,param$gamma,log=T),nloc,ncommun)
  
  for (j in 1:(ncommun-1)){ #last column has to be 1
    v.new=v.old
    v.new[,j]=novos[,j]
      
    theta.old=convertSBtoNormal(vmat=v.old,ncol=ncommun,nrow=nloc,prod=rep(1,nloc))
    theta.new=convertSBtoNormal(vmat=v.new,ncol=ncommun,nrow=nloc,prod=rep(1,nloc))

    #contribution from reflectance data
    pold=theta.old%*%param$omega
    pnew=theta.new%*%param$omega
    p1.old=rowSums(dbinom(remote,size=ndig.values,pold,log=T))
    p1.new=rowSums(dbinom(remote,size=ndig.values,pnew,log=T))
      
    k=acceptMH(p1.old+prior.old[,j],
               p1.new+prior.new[,j]+ajuste[,j],
               v.old[,j],v.new[,j],F)
    v.old[,j]=k$x
  }
  theta=convertSBtoNormal(vmat=v.old,ncol=ncommun,nrow=nloc,prod=rep(1,nloc))
  list(theta=theta,v=v.old,accept=v.old!=v.orig)
}
