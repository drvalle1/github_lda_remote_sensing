boxplot(param$theta)

plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(param$theta[,i],col=i)

plot(vec.gamma,type='l')

plot(param$omega,omega.true)
# unique(param$omega-omega.true)

