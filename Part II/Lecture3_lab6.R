

#Labs 6.6 and 6.7
library(ISLR) #hitter dataset
library(glmnet)#lasso and ridge
library(pls)#pcr and plsr
library(ggplot2)
library(reshape2)
#install.packages("ISLR")
#install.packages("ggplot2")
#install.packages("leaps")
#install.packages("glmnet")
#install.packages("pls")
names(Hitters)
dim(Hitters)#dimensionality of the data frame
sum(is.na(Hitters$Salary))#this shows how many NAs there are
Hitters=na.omit(Hitters)#this removes all NA lines. 
dim(Hitters)#dimensionality of the data frame without NA values
#check the variables ranges and types
summary(Hitters)
#are there correlations between variables?
cormat =cor(Hitters[,c(-15,-14,-20)])

ggplot(data = melt(cormat), aes(x=Var1, y=Var2, fill=value)) + 
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation")+ geom_tile()
pairs(Hitters[,1:3])

#RIDGE
x=model.matrix(Salary~.,Hitters)[,-1]
y=Hitters$Salary
#compare x to the original file. 
x[1:3,]
Hitters[1:3,]

#define series of values for lambda over large range, then perform ridge fit; alpha=0, lasso would be alpha=1
grid=10^seq(10,-2,length=100) 
ridge.mod=glmnet(x,y,alpha=0,lambda=grid)

dim(coef(ridge.mod)) #includes intercepts
ridge.mod$lambda[50] #50th lambda value
coef(ridge.mod)[,50] #respective coefficients

coef=coef(ridge.mod)[-1,] #no intercept
l2=sqrt(apply(coef*coef,2,sum))
l2[50] #length of the coefficient vector
coef[,50]
l2[60] # smaller lambda so larger coefficients
coef[,60]
plot(log(ridge.mod$lambda),l2,type='l')
#see something odd? 
#Let's have a look at the coefficients
plot(ridge.mod) #coefficients don't go to 0

#get coefficients for a new lambda value
predict(ridge.mod,s=50,type="coefficients")[1:20,]

#make test and training set
set.seed(1)
train=sample(1:nrow(x),nrow(x)/2)
test=(-train)
y.test=y[test]

#learn model with training data
ridge.mod=glmnet(x[train,],y[train],alpha=0,lambda=grid,
		thresh=1e-12)
ridge.pred=predict(ridge.mod,s=4,newx=x[test,],type="response")#predictions for test and lambda=4
mean((ridge.pred-y.test)^2) #performance in MSE on test
mean((mean(y[train])-y.test)^2)#MSE for a null model with only intercept
ridge.pred=predict(ridge.mod,s=1e10,newx=x[test,])#a null model with only intercept,by forcing coefficients to 0 with large lambda
mean((ridge.pred-y.test)^2) #MSE
ridge.pred=predict(ridge.mod,s=0,newx=x[test,],x=x[train,],y=y[train],exact=TRUE) #least squares fit, lambda =0
mean((ridge.pred-y.test)^2) #MSE
lm(y~x,subset=train) #least squares fit coefficients
predict(ridge.mod,s=0,exact=TRUE,x=x[train,],y=y[train],type="coefficients")[1:20,] #least squares fit coefficients by setting lambda 0

#find optimal value for lambda with cross validattion, ten-fold by default
set.seed(1)
cv.out=cv.glmnet(x[train,],y[train],alpha=0)
plot(cv.out)#the variability in the CV is shown as well as two thresholds
minlam=cv.out$lambda.min
selam=cv.out$lambda.1se
ridge.pred=predict(ridge.mod,s=minlam,newx=x[test,])
mean((ridge.pred-y.test)^2)
ridge.pred=predict(ridge.mod,s=selam,newx=x[test,])
mean((ridge.pred-y.test)^2)

out=glmnet(x,y,alpha=0)
predict(out,type="coefficients",s=minlam)[1:20,]#coefficients fit w. best lambda on full dataset

#LASSO
lasso.mod=glmnet(x[train,],y[train],alpha=1,lambda=grid) #alpha is 1 for lasso!
#observe the % of deviance explained over different value of lambda and # of components
print(lasso.mod)
plot(lasso.mod)#the coefficents go to zero
plot(lasso.mod,xvar="lambda",label=TRUE)#behaviour coef plotted against log lambda
plot(lasso.mod,xvar="dev",label=TRUE)#behaviour coef plotted against percentage deviance explained
#perform cross-validation for lambda, then measure performance and look at selected coefficients
set.seed(1)
cv.out=cv.glmnet(x[train,],y[train],alpha=1)
plot(cv.out)
bestlam=cv.out$lambda.min
lasso.pred=predict(lasso.mod,s=bestlam,newx=x[test,])
mean((lasso.pred-y.test)^2)
out=glmnet(x,y,alpha=1,lambda=grid)
lasso.coef=predict(out,type="coefficients",s=bestlam)[1:20,]
lasso.coef
#PCR
set.seed(2)
pcr.fit=pcr(Salary~.,data=Hitters,scale=TRUE,
		validation="CV")
summary(pcr.fit)#see the variance explained by the PC's
explvar(pcr.fit) #extract the variance explained per PC
plot(pcr.fit, plottype = "scores", comps = 1:3)#scores for samples per PC;do you see outliers, or groups, other patterns 

validationplot(pcr.fit,val.type="MSEP") #shows the MSE to the number of coefficients

#you can observe the loadings across the 19 variables for first 2 PC's
plot(pcr.fit, "loadings", comps = 1:2, legendpos = "topright",labels = 1:2)
abline(h = 0)

set.seed(1)
pcr.fit=pcr(Salary~.,data=Hitters,subset=train,scale=TRUE,
		validation="CV")
validationplot(pcr.fit,val.type="MSEP") #select 7 components
coefplot(pcr.fit)
pcr.pred=predict(pcr.fit,x[test,],ncomp=7)
mean((pcr.pred-y.test)^2)
pcr.fit=pcr(y~x,scale=TRUE,ncomp=7)
summary(pcr.fit)



