library(leaps)
library(glmnet)
#Chapter 6 exercise 8
#In this exercise, we will generate simulated data, and will then use this data to perform best subset selection.
#(a) Use the rnorm() function to generate a predictor X of length n = 100, as well as a noise vector of length n = 100.
set.seed(1)
X <- rnorm(100)
error <- rnorm(100,sd=4)

#(b) Generate a response vector Y of length n = 100 for a cubic polynomial, with coefficients for all parts.
#Choose values for the coefficients
B0= 1
B1= 2
B2= 3
B3=-4
Y <- B0 + B1*X + B2*X^2 + B3*X^3 + error
plot(X,Y) 
#Pretty! You can vary the influence of the error (increase the sd above) to explore when the approach fails!

#(c) Use regsubsets() to perform best subset selection in order to choose the best model.
#What is the best model obtained according to Cp, BIC, and adjusted R2? 
#Show some plots to provide evidence for your answer, and report the coefficients of the best model ob- tained.

data <- data.frame(x=X,Y=Y)
subset <- regsubsets(Y~poly(X,10,raw=TRUE), data = data)
#Selecting raw in the poly function prevents poly from generating orthogal polynomials!
#Reverting to raw basically generates the polynomials you expect.

subsetSummary <- summary(subset)
plot(subsetSummary$bic)
subsetSummary$bic # 3 predictors, which is correct! 
plot(subsetSummary$cp)
subsetSummary$cp # 3 predictors, which is correct!
plot(subsetSummary$adjr2)
subsetSummary$adjr2 # 4 predictors, which is not correct!
coef(subset,3)#Quite close to the real values! 
coef(subset,4)#Oops

#8d.Repeat (c), using forward stepwise selection and also using backwards stepwise selection. 
#How does your answer compare to the results in (c)?
fwd <- regsubsets(Y~poly(X,10,raw=TRUE),data=data,method="forward")
fwdSUM <- summary(fwd)
plot(fwd)
plot(fwdSUM$bic)
subsetSummary$bic # 3 predictors, which is correct! 
plot(fwdSUM$cp)# 4 predictors, which is not correct!
plot(fwdSUM$adjr2)# 4 predictors, which is correct!
coef(fwd,3)
bwd=regsubsets(Y~poly(x,10,raw=TRUE),data=data,method="backward")
bwdSUM <- summary(bwd)
plot(bwd)
#in this example you see that the second coefficient is accidentally excluded early in the procedure
#this means it remains excluded as you go to simpler models
plot(bwdSUM$bic)#3 correct number
plot(bwdSUM$cp) #3 correct number
plot(bwdSUM$adjr2) #5! too high
coef(bwd,3)
#Again some statistics find the correct complexity, others don't! However the 2nd variable is incorrectly excluded
#The results will change if you change the error amplitude.

#(e) Fit a lasso model. Use cross-validation to select the optimal value of λ. Create plots of the cross-validation error as a function of λ. Report the resulting coefficient estimates.
#Discuss the results obtained.
#mat <- model.matrix(Y~poly(X,10,raw=TRUE),data=data)[,-1]
set.seed(10)
mat=model.matrix(Y~poly(X,10,raw=TRUE))
#The model.matrix() function is particularly useful for creating x; not only does it produce a matrix corresponding to the 19 predictors but it also automatically transforms any qualitative variables into dummy variables. The latter property is important because glmnet() can only take numerical, quantitative inputs.
lasso <- cv.glmnet(mat, Y, alpha = 1)
lambda <- lasso$lambda.1se #select the simplest model, using the most restrictive lambda with good performance
plot(lasso)
lassoModel <- glmnet(mat, Y, alpha = 1)
predict(lassoModel, s =lambda, type = "coefficients")
#depending on the error amplitude and seed here and above you will find the correct set of variables or not!
#here a small coefficient is given to X^5 but not to X^1
#Let's plot this fit now, how good is this?

Xr=seq(from=-2.5,to=2.5,by=0.1)
Yr=B0 + B1*Xr + B2*Xr^2 + B3*Xr^3
Yl=1.1129360+2.7315556*Xr^2+-2.501508*Xr^3+-0.2018969*Xr^5
plot(X,Y,pch=19,col="grey") #points in grey
lines(Xr,Yr,col=1,lty=1) #real formula in solid black 
lines(Xr,Yl,col=2,lty=2) #lasso fit in dotted red, looking pretty good
#Is this a consequence of shrinking the coefficients? The lasso fit has a coefficient vector with an l1 norm of 5.43 vs 9 for the real coefficients

#(f)Now generate a response vector Y according to the model Y = β0 + β7X7 + ε,
#and perform best subset selection and the lasso. Discuss the results obtained.

B7 <- 0.1
Y <- B0 + B7*X^7 + error
plot(X,Y) #almost a straight horizontal, with a few observations to note the deviations from linear at the edges
data <- data.frame(x=X,Y=Y)
subset <- regsubsets(Y~poly(X,10,raw=TRUE), data = data, nvmax=10)
subsetSummary <- summary(subset)
names(subsetSummary)
plot(subsetSummary$cp)
plot(subsetSummary$bic)
plot(subsetSummary$adjr2)
#all select a model with a single value
coefficients(subset, id = 1) #this is an amazingly good fit! Nearly perfect.

#Lasso 
mat <- model.matrix(Y~poly(X,10,raw=TRUE),data=data)[,-1]
lasso <- cv.glmnet(mat, Y, alpha = 1)
lambda <- lasso$lambda.1se
plot(lasso)
lassoModel <- glmnet(mat, Y, alpha = 1)
plot(lassoModel)
predict(lassoModel, s =lambda, type = "coefficients")
#This is quite a good fit! Again in this case not the perfect model is fit for lasso, with a small coefficient for X^9.
#the l1 norm for the B is 0.074 for lasso vs 0.1 in reality
#It is true that the X^7 and X^9 are highly correlated, but X^9 grows faster and requires a smaller coefficient
plot(Xr^7,Xr^9)
cor(Xr^7,Xr^9)

