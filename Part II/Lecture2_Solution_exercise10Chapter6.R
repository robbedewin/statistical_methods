#Exercise 10 Chapter 6
#We have seen that as the number of features used in a model increases, the training error will necessarily decrease, but the test error may not. 
#We will now explore this in a simulated data set.
#(a) Generate a data set with p = 20 features, n = 1,000 observations, and an associated quantitative response vector generated according to the model
#Y =Xβ+ε,
#where β has some elements that are exactly equal to zero.
set.seed(1)
library (leaps)


p = 20
n = 1000
x = matrix(rnorm(n*p),n,p)
B = rnorm(p) 
B[c(2,3,8,9,10,15,20)] = 0 
e = rnorm(n)
y = x %*% B + e
#(b) Split your dataset into a trainingset containing 100 observations 
#and a test set containing 900 observations.
train = sample(seq(1000),100,replace=F)
#(c) Perform best subset selection on the training set.
#Plot the training set MSE associated with the best model of each size.
subset = regsubsets(y~.,data=data.frame(x=x[train,],y=y[train]),nvmax=p)
plot(summary(subset)$rss/100)#The training MSE (=RSS/n) for different sizes of the subset
#(d) Plot the test set MSE associated with the best model of each size.
newdata=data.frame(x=x[-train,],y=y[-train])
#the following function is from the lab
predict.regsubsets =function (object ,newdata ,id ,...){
	form=as.formula(object$call[[2]])
	mat=model.matrix(form,data=newdata)
	coefi=coef(object ,id=id)
	xvars=names(coefi)
	mat[,xvars]%*%coefi
}
val.errors = rep(0,p)
for (i in 1:p) {
	pred = predict(subset,newdata=newdata, id = i)# as in the lab 
	val.errors[i] = mean((y[-train]-pred)^2)
}
plot(val.errors ,ylab="test MSE",pch=19,type="b") # it's beautiful

#(e) For which model size does the test set MSE take on its minimum value? 
#Comment on your results.
which.min(val.errors) #14 variables, it should be 13!
#From the evolution of the test MSE it is also clear the increase of the test error is not large with more variables. So what is going on?

#(f)How does the model at which the test set MSE is minimized compare to the true model used to generate the data? Comment on the coefficient values.
#Let's compare the coefficients of the best model with the true B 
coef(subset,id=which.min(val.errors)) 
B
#We find an intercept and variable 10 has a small coefficient instead of 0.
#Some smaller coefficients maybe missed in your run, their signal probably small compared to the noise.
#The overlapping coefficients show a  good fit however.

#(g) Create a plot displaying j=1(βj − βj^r ) for a range of values
#of r, where βˆjr is the jth coefficient estimate for the best model containing r coefficients. Comment on what you observe.
#How does this compare to the test MSE plot from (d)?

names=paste("x.",1:20,sep="")
result=rep(-1,p)
for (i in 1:p) {
	coefi = coef(subset,id=i)[-1]
	vals=rep(0,p)
	vals[match(names(coefi),names)]=coefi
	result[i]=sqrt(sum((B-vals)^2))
}
plot(result,type="b",xlab="Subset size", ylab="Coefficient Error")
coef(subset,id=13) 
#The plot shows a very similar pattern to the one shown in the test MSE, with a tiny dip in coefficient error when
#the 14 parameters are fit. For the correct number of parameters the overall error in the coefficients is slightly higher.
#The parameters are harder to fit correctly as the number of parameters increases. In lower dimensions the parameters can be fit with less variance,
#which in this cases wins over a potential decrease in bias.   
#The similarity to the test MSE plot is no surprise as the plots reflects how closely the estimated formula approximates 
#the true underlying function, and this should correlate with the predictive performance.
