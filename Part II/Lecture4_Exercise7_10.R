library(ISLR) 
library(leaps) 
library(gam)

#A. Split the data, and forward stepwise selection

set.seed(1)
dim(College)
train = sample(1:dim(College)[1],dim(College)[1]*1/2) 
test = (-train)
summary(College)
pairs(College[,1:10])

regfit.fwd=regsubsets(Outstate~. ,College[train,],method="forward")
summary(regfit.fwd)
plot(summary(regfit.fwd)$bic,type="b",ylab="BIC")#we select 6 variables
plot(regfit.fwd, scale="bic")
lm1 = lm(Outstate~Private +Room.Board +Terminal+ perc.alumni+Expend+Grad.Rate ,data=College,subset=train)
predlm1 = predict(lm1, newdata=College[test,]) 
mselm1 = mean((predlm1-College[test,"Outstate"])^2)#worse as expected

#B-D. Fit a GAM, plot the results, evaluate the model. Are there non-linear effects? 
gam1 = gam(Outstate~Private +s(Room.Board,4) +s(Terminal,4)+ s(perc.alumni,4)+s(Expend,4)+s(Grad.Rate,4) ,data=College,subset=train) 
#ignore the gam "non-list contrasts" warning; it's a (harmless) bug
par(mfrow=c(2,3))
plot(gam1,se=TRUE,col="purple")#Room.board, alumni, Grad.rate mostly linear, Terminal non-linear effect has rel. high error, expend looks non-linear
summary(gam1)#only expend seems to have non-linear effect
predgam = predict(gam1, newdata=College[test,]) 
msegam1 = mean((predgam-College[test,"Outstate"])^2)
gam2 = gam(Outstate~Private +Room.Board +Terminal+ perc.alumni+s(Expend,4)+Grad.Rate ,data=College,subset=train) 
plot(gam2,se=TRUE,col="purple")#Room.board mostly linear, PhD non-linear effect has v. high error, alumni linear, expend looks non-linear, grad.rate mostly linear
summary(gam2)#only expend seems to have non-linear effect
predgam2 = predict(gam2, newdata=College[test,]) 
msegam2 = mean((predgam2-College[test,"Outstate"])^2)#essentially identical to msegam1 for my train/test
anova(gam1,gam2)#simplification justified as expected.


lm2 = lm(Outstate~Private +Room.Board +Terminal+ perc.alumni+ns(Expend,4)+Grad.Rate ,data=College,subset=train)
predlm2 = predict(lm2, newdata=College[test,]) 
mselm2 = mean((predlm2-College[test,"Outstate"])^2)#marginally worse than msegam2


