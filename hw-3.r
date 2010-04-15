## Question 2
## data input
y<-c(28,8,-3,7,-1,1,18,12)
sigma<-c(15,10,16,11,9,11,10,18)
ngrid<-1000

## Marginal log-Posterior for tau^2

logpost.tau<-function(tausq,y,sigma){
  vsq<-sum(1/(sigma^2+tausq))
  beta<-sum(y/(sigma^2+tausq))*vsq
  a<--1/2*log(tausq)
  b<-sum(log(dnorm(y,1,sqrt(sigma^2+tausq))))
  c<-log(dnorm(1,beta,vsq))
  value<-a+b-c
  return(value)
}
tausq<-seq(0.001,10,len=ngrid)
logpost.tausq<-rep(NA,ngrid)
for (i in 1:ngrid){
  logpost.tausq[i]<-logpost.tau(tausq[i],y,sigma)
}

post.tausq<-exp(logpost.tausq-max(logpost.tausq))
post.tausq<-post.tausq/sum(post.tausq)


## Conditional posterior of mu
post.condmu<-function(tau,sigma,y){
  V<-1/sum(1/(sigma^2+tau^2))
  beta<-sum(y/(sigma^2+tau^2))*V
  mu<-rnorm(1,beta,sqrt(V))
  value<-dnorm(mu,beta,sqrt(V))
  return(value)
}
post.cond.mu<-rep(NA,ngrid)
for (i in 1:ngrid){
  post.cond.mu[i]<-post.condmu(tausq[i],sigma,y)
}

## Joint posterior of mu and tau^2
logjointpost<-function(mu,sigma,tausq,y){
  a<--1/2log(tausq)
  b<-sum(dnorm(y,mu,sigma^2+tausq))
  value<-a+b
  return(value)
}
#joint.post<-rep(NA,ngrid)
#for (i in 1:ngrid){
#  joint.post[i]<-post.tausq[i]*post.cond.mu[i]
#}
#joint.post<-joint.post/sum(joint.post)

## Question 3

## grid sampling: sample 1000 values of tau^2

tausq.sample<-sample(tausq,size=1000,replace=T,prob=post.tausq)

## Sampling mu given sampled tau^2

mu.sample<-rep(NA,ngrid)
for (i in 1:ngrid){
  vsq<-sum(1/(sigma^2+tausq.sample[i]))
  beta<-sum(y/(sigma^2+tausq.sample[i]))*vsq
  mu.sample[i]<-rnorm(1,beta,sqrt(vsq))
}

mean(mu.sample)
quantile(mu.sample,c(0.025,0.5,0.975))
mean(tausq.sample)
quantile(tausq.sample,c(0.025,0.5,0.975))

## Question 4
theta.sample<-function(mu,tausq,sigma,y){
  alpha<-rep(NA,8)
  beta<-rep(NA,8)
  value<-rep(NA,8)
  for (i in 1:8){
    beta[i]<-1/(1/tausq+1/sigma[i]^2)
    alpha[i]<-(mu/tausq+y[i]/sigma[i]^2)*beta[i]
    value[i]<-rnorm(1,alpha[i],sqrt(beta[i]))
  }
  return(value)
}
theta.samp<-matrix(NA,nrow=ngrid,ncol=8)
for (i in 1:ngrid){
  theta.samp[i,]<-theta.sample(mu.sample[i],tausq.sample[i],sigma,y)
}
mean<-apply(theta.samp,2,mean)

## Question 6
tausq0<-median(tausq.sample)
theta.new.samp<-matrix(NA,nrow=ngrid,ncol=8)
mu.new.samp<-rep(NA,ngrid)
new.samp<-matrix(NA,nrow=ngrid,ncol=9)
start<-c(8,5,5,5,5,5,5,5,5)
#start<-c(8,10,10,10,10,10,10,10,10)  #different starting value
new.sample<-function(tausq0,y,sigma){
  param<-matrix(NA,ngrid,9)
  colnames(param)<-c("mu","theta1","theta2","theta3","theta4","theta5","theta6","theta7","theta8")
  xi<-rep(NA,8)
  gamma<-rep(NA,8)
  param[1,]<-start
  for (i in 2:ngrid){
    param[i,1]<-dnorm(1,sum(param[i-1,2:9])/8,tausq0/8)
    for (j in 1:8){
      gamma[j]<-1/(1/tausq0+1/sigma[j]^2)
      xi[j]<-(param[i-1,1]/tausq0+y[j]/sigma[j]^2)*gamma[j]
      param[i,j+1]<-rnorm(1,xi[j],sqrt(gamma[j]))
    }
  }
  return(param)
}
new.samp<-new.sample(tausq0,y,sigma)
t(apply(new.samp,2,quantile,c(0.025,0.5,0.975)))
apply(new.samp,2,mean)

par(mfrow=c(2,4))
for (i in 2:9){
  plot(new.samp[,i],type="l")
}

par(mfrow=c(2,4))
for (i in 2:9){
  acf(new.samp[,i])
}
t(apply(new.samp,2,quantile,c(0.025,0.5,0.975)))
apply(new.samp,2,mean)

## Trying different starting value
start<-c(8,10,10,10,10,10,10,10,10)


## Question 7
# draw 1000 samples for each school
new.ngrid<-1000
y.pred<-matrix(NA,nrow=new.ngrid,ncol=8)
theta.samp<-new.samp[,2:9]
for (i in 1:1000){
  for (j in 1:8){
  y.pred[i,j]<-rnorm(1,theta.samp[i,j],sqrt(sigma[j]))
  }
}
n<-0
for (i in 1:new.ngrid){
  if (y.pred[i,1]>max(y.pred[i,2:8]))
  n<-n+1
}
n/new.ngrid

## Simulation 1000 for the situation above
new.ngrid<-1000
y1.pred<-matrix(NA,ngrid,ngrid)
y2.pred<-matrix(NA,ngrid,ngrid)
y3.pred<-matrix(NA,ngrid,ngrid)
y4.pred<-matrix(NA,ngrid,ngrid)
y5.pred<-matrix(NA,ngrid,ngrid)
y6.pred<-matrix(NA,ngrid,ngrid)
y7.pred<-matrix(NA,ngrid,ngrid)
y8.pred<-matrix(NA,ngrid,ngrid)
prediction<-function(theta,sigma){
  value<-rnorm(1,theta,sqrt(sigma))
  return(value)
}
for (i in 1:ngrid){
  for (j in 1:ngrid){
    y1.pred[i,j]<-prediction(theta.samp[i,1],sqrt(sigma[1]))
  }
}
for (i in 1:ngrid){
  for (j in 1:ngrid){
    y2.pred[i,j]<-prediction(theta.samp[i,2],sqrt(sigma[2]))
  }
}
for (i in 1:ngrid){
  for (j in 1:ngrid){
    y3.pred[i,j]<-prediction(theta.samp[i,3],sqrt(sigma[3]))
  }
}
for (i in 1:ngrid){
  for (j in 1:ngrid){
    y4.pred[i,j]<-prediction(theta.samp[i,4],sqrt(sigma[4]))
  }
}
for (i in 1:ngrid){
  for (j in 1:ngrid){
    y5.pred[i,j]<-prediction(theta.samp[i,5],sqrt(sigma[5]))
  }
}
for (i in 1:ngrid){
  for (j in 1:ngrid){
    y6.pred[i,j]<-prediction(theta.samp[i,6],sqrt(sigma[6]))
  }
}
for (i in 1:ngrid){
  for (j in 1:ngrid){
    y7.pred[i,j]<-prediction(theta.samp[i,7],sqrt(sigma[7]))
  }
}
for (i in 1:ngrid){
  for (j in 1:ngrid){
    y8.pred[i,j]<-prediction(theta.samp[i,8],sqrt(sigma[8]))
  }
}
best.times<-rep(0,ngrid)
for (j in 1:ngrid){
  for (i in 1:ngrid){
    if (y1.pred[i,j]>max(y2.pred[i,j],y3.pred[i,j],y4.pred[i,j],y5.pred[i,j],y6.pred[i,j],y7.pred[i,j],y8.pred[i,j]))
    {best.times[j]<-best.times[j]+1}
  }
}
## probability that school A offers the best program
prob<-best.times/ngrid
mean(prob)
quantile(prob,c(0.025,0.5,0.975))

##Question 11
y<-c(74,99,58,70,122,77,104,129,308,119)
m<-length(y)

##posterior log-likelohood of alpha and beta
loglike<-function(phi,theta){
  if (phi[1]<=0 | phi[2]<=0) return(1e-10)
  else {
  value<-sum(dgamma(theta,phi[1],rate=phi[2],log=T))
  return(value)
  }
}
## Gibbs Sampler 1: Independent log-Normal Distribution
nsim<-10000
##phi<-c(3,3)
phi<-c(10,10)
theta<-y
##sd.phi<-c(.5,.5)
##sd.phi<-c(.35,.35)
sd.phi<-c(.75,.75)
nacc<-0
theta.norm<-matrix(NA,nsim,m)
phi.norm<-matrix(NA,nsim,2)
for (i in 1:nsim){
  ##update theta
  theta<-rgamma(m,phi[1]+y,phi[2]+1)
  ##update alpha and beta
  new<-exp(rnorm(2,log(phi),sd.phi))
  log.new<-loglike(new,theta)+sum(log(new))
  log.old<-loglike(phi,theta)+sum(log(phi))
  log.diff<-log.new-log.old
  logu<-log(runif(1))
  if (logu<log.diff) {phi<-new; nacc<-nacc+1}
  theta.norm[i,]<-theta
  phi.norm[i,]<-phi
}
#phi1<-phi.norm
#phi2<-phi.norm
phi3<-phi.norm
labels<-expression(alpha,beta)

## Acf plot for phi1
par(mfrow=c(2,1))
for (i in 1:2) acf(phi1[,i],lag.max<-100)

## acf for phi2
par(mfrow=c(2,1))
for (i in 1:2) acf(phi2[,i],lag.max<-100)

## acf for phi3
par(mfrow=c(2,1))
for (i in 1:2) acf(phi3[,i],lag.max<-100)


## Compare chains
par(mfrow=c(2,1))
plot(phi1[,1],type="l",ylim=range(phi1[,1],phi2[,1],phi3[,1]))
lines(phi2[,1],col=2)
lines(phi3[,1],col=3)

plot(phi1[,2],type="l",ylim=range(phi1[,2],phi2[,2],phi3[,2]))
lines(phi2[,2],col=2)
lines(phi3[,2],col=3)

summary(phi1[,1])
summary(phi2[,1])
summary(phi3[,1])
quantile(phi1[,1],c(0.025,0.0975))
quantile(phi2[,1],c(0.025,0.0975))
quantile(phi3[,1],c(0.025,0.0975))
summary(phi1[,2])
summary(phi2[,2])
summary(phi3[,2])
quantile(phi1[,2],c(0.025,0.0975))
quantile(phi2[,2],c(0.025,0.0975))
quantile(phi3[,2],c(0.025,0.0975))

##Running a longer chain

nsim<-101000
##phi<-c(3,3)
phi<-c(10,10)
theta<-y
##sd.phi<-c(.5,.5)
##sd.phi<-c(.35,.35)
sd.phi<-c(.75,.75)
nacc<-0
theta.norm<-matrix(NA,nsim,m)
phi.norm<-matrix(NA,nsim,2)
for (i in 1:nsim){
  ##update theta
  theta<-rgamma(m,phi[1]+y,phi[2]+1)
  ##update alpha and beta
  new<-exp(rnorm(2,log(phi),sd.phi))
  log.new<-loglike(new,theta)+sum(log(new))
  log.old<-loglike(phi,theta)+sum(log(phi))
  log.diff<-log.new-log.old
  logu<-log(runif(1))
  if (logu<log.diff) {phi<-new; nacc<-nacc+1}
  theta.norm[i,]<-theta
  phi.norm[i,]<-phi
}

labels<-expression(alpha,beta)
par(mfrow=c(2,1))
plot(phi.norm[,1],type="l",ylab=labels[1])
plot(phi.norm[,2],type="l",ylab=labels[2])

par(mfrow=c(2,1))
acf(phi.norm[,1],lag.max=100)
acf(phi.norm[,2],lag.max=100)


alpha<-phi.norm[,1]
beta<-phi.norm[,2]
mu<-alpha/beta
sd<-sqrt(alpha/beta^2)

par(mfrow=c(2,2))
hist(alpha)
hist(beta)
hist(mu)
hist(sd)

##Examining Shrinkage
par(mfrow=c(1,1))
theta.postmean<-apply(theta.norm,2,mean)
mu.postmean<-mean(mu)
plot(1:m,y,main="shrinkage of sample proportions",pch=19)
abline(h=mu.postmean,col=4,lwd=2)
points(1:m,theta.postmean,pch=19,col=2)
legend(2,250,c("Sample","Post Mean","Overal Mean"),pch=19,col=c(1,2,4))


temp<-seq(1001,101000,by=100)
plot(phi.norm[temp,1],type="l",ylab=labels[1])
plot(phi.norm[temp,2],type="l",ylab=labels[2])

acf(phi.norm[temp,1],lag.max=100)
acf(phi.norm[temp,2],lag.max=100)




## 1000 samples of alpha and beta
alpha.samp<-alpha[temp]
beta.samp<-beta[temp]
mean(alpha.samp)
quantile(alpha.samp,c(.025,.5,.975))
mean(beta.samp)
quantile(beta.samp,c(.025,.5,.0975))

## Question 12
##posterior predictive for a new residential street with a bike lane

theta.new<-rgamma(1000,alpha.samp,rate=beta.samp)
y.new<-rpois(1000,theta.new)
par(mfrow=c(2,1))
hist(theta.new)
abline(v=mean(theta.new),col=2)
hist(y.new)
abline(v=mean(y.new),col=2)
mean(y.new)
quantile(y.new,c(.025,.50,.975))
