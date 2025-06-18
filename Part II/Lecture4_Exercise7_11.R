set.seed(40)
x1=rnorm(100)
x2=rnorm(100)*1.5+2-x1/2
y=1+x1-0.2*x2+0.01*x1**3+rnorm(100)
par(mfrow=c(1,1))
plot(x1,y)
plot(x2,y)
plot(x2,x1)

b0=rep(0,1000)
b1=rep(0,1000)
b2=rep(0,1000)
beta1=10
for (i in 1:1000){
  a =y - beta1 * x1
  beta2 = lm ( a~x2 ) $coef [2]
  a =y - beta2 * x2
  lmT = lm ( a~x1 )
  beta1=lmT$coef[2]
  b0[i]=lmT$coef[1]
  b1[i]=beta1
  b2[i]=beta2
}
par(mfrow=c(1,1))
plot(1:10,b0[1:10],col=1,type="l",ylim=c(-2.5,1.5),ylab="Coefficient value",xlab="iteration")#it is basically stable after 10 iterations
abline(h=0,col=1,lty=2)
lines(1:10,b1[1:10],col=2)
lines(1:10,b2[1:10],col=3)
legend(8,-1.5,legend=c("B0","B1","B2"),col=c(1,2,3),lty=1,cex=0.7,bg="grey95",box.lty=0)
summary(lm(y~x1+x2))
b0[1000];b1[1000];b2[1000]
          