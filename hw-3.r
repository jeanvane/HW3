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

Question 11
y<-c(74,99,58,70,122,77,104,129,308,119)
m<-length(y)
## posterior log-likelihood of alpha and beta
loglike<-function(alpha,beta,theta){
  if (beta<=0) return(-1e10)
  if (beta>0){
  n<-length(theta)
  value<-n*alpha*log(beta)-n*digamma(alpha)
  value<-value+(alpha-1)*sum(log(theta))-beta*sum(theta)
  return(value)
  }
}
nsim<-10000
alpha.samp<-rep(NA,nsim)
beta.samp<-rep(NA,nsim)
theta.samp<-matrix(NA,10000,m)
alpha.samp[1]<-3
beta.samp[1]<-0.03
for(j in 1:m){
  theta.samp[1,j]<-rgamma(1,alpha.samp[1]+y[j],rate=beta.samp[1]+1)
}
for (i in 2:nsim){
  alpha.star<-rnorm(1,alpha.samp[i-1],1)
  beta.star<-rnorm(1,beta.samp[i-1],1)
  log.l<-loglike(alpha.samp[i-1],beta.samp[i-1],theta.samp[i-1,])
  log.star<-loglike(alpha.star,beta.star,theta.samp[i-1,])
  logr<-log.star-log.l
  u<-runif(1)
  logu<-log(u)
  if (logu<=logr){
    alpha.samp[i]<-alpha.star
    beta.samp[i]<-beta.star
  }
  if (logu>logr){
    alpha.samp[i]<-alpha.samp[i-1]
    beta.samp[i]<-beta.samp[i-1]
  }
  for (j in 1:m){
    theta.samp[i,j]<-rgamma(1,alpha.samp[i]+y[j],rate=beta.samp[i]+1)
  }
}

##Examing the samples
par(mfrow=c(2,1))
plot(1:nsim,alpha.samp,type="l")
plot(1:nsim,beta.samp,type="l")

