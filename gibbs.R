gibbs.RS=function(nloc,nbands,ndig.values,ngibbs,ncommun,
                  gamma,a.omega,b.omega){
  
  #initial values
  omega=matrix(runif(ncommun*nbands),ncommun,nbands)
  theta=matrix(1/ncommun,nloc,ncommun)
  v=theta
  v[,ncommun]=1
  
  #stuff for gibbs sampling
  param=list(theta=theta,omega=omega,v=v,gamma=gamma)
  vec.theta=matrix(NA,ngibbs,nloc*ncommun)
  vec.omega=matrix(NA,ngibbs,ncommun*nbands)
  
  #stuff for MH algorithm
  jump1=list(omega=matrix(1,ncommun,nbands),
             v=matrix(0.3,nloc,ncommun))
  accept1=list(omega=matrix(0,ncommun,nbands),
               v=matrix(0,nloc,ncommun))
  accept.output=50
  
  for (i in 1:ngibbs){
    print(i)
    tmp=update.theta(param,jump1$v,ncommun,nloc,ndig.values)
    param$theta=tmp$theta #theta.true#
    param$v=tmp$v
    accept1$v=accept1$v+tmp$accept
    
    tmp=update.omega(param,jump1$omega,ncommun,nbands,ndig.values,a.omega,b.omega)
    param$omega=tmp$omega #omega.true#
    accept1$omega=accept1$omega+tmp$accept
    
    if (i%%accept.output==0 & i<1000){
      k=print.adapt(accept1,jump1,accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
    
    vec.theta[i,]=param$theta
    vec.omega[i,]=param$omega
  }
  list(theta=vec.theta,omega=vec.omega)
}
  
