# This program contains all the code used to generate the outcomes employed in the lecture of 
# the  Mising Data Course. Including the different analysis of the case study and all the 
# simulations.

## Run only one time to install the packages

install.packages("mice")
install.packages("lattice")
install.packages("VIM")
install.packages("aod") # neccesary for the Wald test
install.packages("xtable")
install.packages("BaM")  # neccesary for generating multivariate normals in the simulations
install.packages("MASS")
install.packages("nnet")

## Run everyy time to upload the packages

library(mice)
library(lattice)
library(VIM)
library(aod)
library(xtable)
library(BaM)

## Defining the working directory

setwd("C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Bioinformatics\\Missing-Data\\R-code")


##################################### This does not  need to be done anymore
#########################################################################################

## Reading the data

titanic.missing<-read.table("http://lib.stat.cmu.edu/S/Harrell/data/ascii/titanic.txt",sep=',', header=TRUE) 
titanic.missing<-subset(titanic.missing,select=c('survived','pclass','sex','age'))
titanic.missing$sex<-as.numeric(titanic.missing$sex)-1

## Making some sex missing

set.seed(4321)
titanic.missing$sex[sample(nrow(titanic.missing), 10)] <- NA

## Saving the data

write.table(titanic.missing, file="titanicmissing.txt", sep=",")

#########################################################################################
#########################################################################################

## Reading the data
 
titanic.missing <- read.table("titanicmissing.txt", header=T, sep=",")
head(titanic.missing,15)


## Exploring the missingness
## Document  VIM-graphical-tools.pdf in the literature

## numbers = TRUE, prop = FALSE plots the frequencies
## numbers = TRUE, prop = TRUE plots the proportions
## numbers = FALSE neither frequencies nor proportions
## combined=FALSE adds a histogram besides with the frequency/proportion of missigness per variable
## combined=TRUE does  not adds a histogram besides 

## Plot with the histogram besides and the frequencies

titanic.missing.aggr=aggr(titanic.missing,numbers = TRUE, prop = FALSE, ylab=c("Histogram of missing data","Pattern"))
titanic.missing.aggr

##  Plot  without  the histogram and the proportions

aggr(titanic.missing, combined=TRUE, numbers = TRUE, prop = TRUE,cex.numbers=0.87, varheight = FALSE)

## barMiss functio

## Barplot with highlighting of missing/imputed values in other variables by splitting each bar
## into two parts. Additionally, information about missing/imputed values in the variable of 
## interest is shown on the right hand side. 

## Click in the left margin to switch to the previous variable or in the
## right margin to switch to the next variable. To regain use of the VIM GUI
## and the R console, click anywhere else in the graphics window.

## If more than one variable is supplied, the bars for the variable of interest are split 
## according to missingness/number of imputed missings in the additional variables. 

## http://rpackages.ianhowson.com/cran/VIM/man/barMiss.html

## Amount of missigness in age for each sex group
## If we click on the right we get a histogram of the amount of missigness in sex for each 
## age group

barMiss(titanic.missing[,c("sex","age")])

## Bart chart or histogram of all variables in the data set. Each colunm is divided in two
## to  show the amount of missing values (in any other variable) for  that colunm

barMiss(titanic.missing)

## histMiss as bartMiss but with histograms

histMiss(titanic.missing)


matrixplot(titanic.missing)


## Fitting a logistic regression model for the complete cases

titanic.logistic.omit<-glm(survived ~ pclass + sex + age, family=binomial, data = titanic.missing)
summary(titanic.logistic.omit)
confint(titanic.logistic.omit)

## Global effect of class

wald.test(b=coef(titanic.logistic.omit), Sigma=vcov(titanic.logistic.omit), Terms=2:3)

## Creating tex files

titanic.logistic.omit.table<-xtable(titanic.logistic.omit,digits=2)
print(titanic.logistic.omit.table)

## Odds ratios

odds.point.omit<-exp(titanic.logistic.omit$coefficients)
odds.CI.omit<-exp(confint(titanic.logistic.omit))

exp(cbind(OR =titanic.logistic.omit$coefficients, confint(titanic.logistic.omit)))

odds.point.omit.table<-xtable(as.data.frame(odds.point.omit),digits=2)
print(odds.point.omit.table)

odds.CI.omit.table<-xtable(odds.CI.omit,digits=2)
print(odds.CI.omit.table)

## Predicted probabilities of survival

prob.survival.male<-data.frame(sex=1,age=mean(titanic.omit$age),pclass=c('1st','2nd','3rd'))
prob.survival.male$ProbSur<-predict(titanic.logistic.omit,newdata=survival.male,type="response")

prob.survival.female<-data.frame(sex=0,age=mean(titanic.omit$age),pclass=c('1st','2nd','3rd'))
prob.survival.female$ProbSur<-predict(titanic.logistic.omit,newdata=prob.survival.female,type="response")

prob.survival.sex<-data.frame(sex=c(0,1,0,1,0,1),age=mean(titanic.omit$age),pclass=c('1st','1st','2nd','2nd','3rd','3rd'))
prob.survival.sex$ProbSur<-predict(titanic.logistic.omit,newdata=prob.survival.sex,type="response")

###### Titanic multiple imputation

## Studying the patterns of missiness

pattern=md.pattern(titanic.missing)
pairs=md.pairs(titanic.missing)

## Imputing the missing values

#titanic.missing$sex=as.factor(titanic.missing$sex)

imp <- mice(titanic.missing, m=100)

imp <- mice(titanic.missing, meth = c("", "", "logreg", "pmm"), m=100)

plot(imp)

## Imputed values for age. Each row corresponds to a missing entry in age. 
## The columns contain the multiple imputations.

imp$imp$age[1:10,1:5]

## Imputed values for sex. Each row corresponds to a missing entry in sex. 
## The columns contain the multiple imputations.

imp$imp$sex[1:10,1:5]

## The complete data combine observed and imputed data.
## The first completed data set can be obtained as (only first 10 passenger shown)

complete(imp,1)[1:10,]

## It is often useful to inspect the distributions of original and the imputed data.
## The complete() function extracts the original and the imputed data sets from the imp object
## as a long (row-stacked) matrix. The col vector separates the observed (blue) and imputed (red) 
## data for age

com <- complete(imp, "long", inc=T)
col <- rep(c("blue","red")[1+as.numeric(is.na(imp$data$age))],101)
stripplot(age~.imp, data=com, jit=TRUE, fac=0.8, col=col, pch=20, cex=1.4, 
          xlab="Imputation number")

## Analyzing the imputed data sets

fit <- with(data=imp, exp=glm(survived ~ pclass + sex + age, family=binomial))
summary(fit)

## Creating a data set with the results of all the analysis

MI.matrix<-matrix(0,100,5)
for(k in 1:100){ MI.matrix[k,]<-coefficients(fit$analyses[[k]])}
MI.results=data.frame(Intercept=MI.matrix[,1], pclass2=MI.matrix[,2], pclass3=MI.matrix[,3], 
                      sex=MI.matrix[,4], age=MI.matrix[,5])
MI.results[1:10,]

## Combining the results using Rubin's rule

est <- pool(fit)
summary(est)
est$qbar

###### Titanic IPW

## Creating the missing data indicator variable r 

titanic.missing$r<-as.numeric(!is.na(titanic.missing$age))*as.numeric(!is.na(titanic.missing$sex))
head(titanic.missing,15)

## Fitting the logistic regression model to  calculate the probabilities of being complete

titanic.ipw.glm<-glm(r ~ pclass + survived, data=titanic.missing,family=binomial)
summary(titanic.ipw.glm)

## Calculating the weights: Inverse Probabilities

titanic.missing$w<-1/fitted(titanic.ipw.glm)
head(titanic.missing,15)

## Fitting the weighted logistic regression

titanic.results.ipw<- glm(survived ~ pclass + sex + age, data=titanic.missing, weights=titanic.missing$w, family=binomial)
summary(titanic.results.ipw)

##############################################################################################
############################# Presentation Simulations I #####################################
##############################################################################################

z<-rmultnorm(50000, matrix(c(0,0),2,1), vmat=matrix(c(1,0.5,0.5,1),2,2))

########### Missing Completely at Random

id<-runif(50000, min=0, max=1)
id1<-(id<0.5)
z.mcar<-z
z.mcar[id1,1]<-NA

# Creating tex files

z.table<-xtable(z[1:10,],digits=2)
z.mcar.table<-xtable(z.mcar[1:10,],digits=2)

#Listwise deletion
mean(z.mcar[,1], na.rm=T)
mean(z.mcar[!id1,2])

cov(z.mcar[,1],z.mcar[,2], use = "complete.obs")
var(z.mcar[,1], na.rm=T)
var(z.mcar[!id1,2])


# Pairwise deletion
mean(z.mcar[,1], na.rm=T)
mean(z.mcar[,2], na.rm=T)

cov(z.mcar[,1],z.mcar[,2], use = "pairwise.complete.obs")
var(z.mcar[,1],use = "complete.obs")
var(z.mcar[,2])


# Mean Imputation

mean.mcar<-z.mcar[,1]
mean.mcar[is.na(mean.mcar)]<-mean(z.mcar[,1], na.rm=T)

mean(mean.mcar, na.rm=T)
mean(z.mcar[,2], na.rm=T)

cov(mean.mcar,z.mcar[,2])
var(mean.mcar)
var(z.mcar[,2])

########### Missing At Random

z.mar<-z
z.mar[z.mar[,2]<0,1]<-NA


# Creating tex files

z.table<-xtable(z[1:10,],digits=2)
z.mar.table<-xtable(z.mar[1:10,],digits=2)

#Listwise deletion

mean(z.mar[,1], na.rm=T)
mean(z.mar[z.mar[,2]>0,2])

cov(z.mar[,1],z.mar[,2], use = "complete.obs")
var(z.mar[,1], na.rm=T)
var(z.mar[z.mar[,2]>0,2])


# Pairwise deletion

mean(z.mar[,1], na.rm=T)
mean(z.mar[,2])

cov(z.mar[,1],z.mcar[,2], use = "pairwise.complete.obs")
var(z.mar[,1],use = "complete.obs")
var(z.mar[,2])


# Mean Imputation

mean.mar<-z.mar[,1]
mean.mar[is.na(mean.mar)]<-mean(z.mar[,1], na.rm=T)

mean(mean.mar)
mean(z.mar[,2])

cov(mean.mar,z.mar[,2])
var(mean.mar)
var(z.mar[,2])


########### Missing Not At Random

z.mnar<-z
z.mnar[z.mnar[,1]<0,1]<-NA


#Listwise deletion

cov(z.mnar[,1],z.mnar[,2], use = "complete.obs")
cor(z.mnar[,1],z.mnar[,2], use = "complete.obs")
var(z.mnar[,1], na.rm=T)
var(z.mnar[!is.na(z.mnar[,1]),2])


# Pairwise deletion

cov(z.mnar[,1],z.mnar[,2], use = "pairwise.complete.obs")
cor(z.mnar[,1],z.mnar[,2], use = "complete.obs")
var(z.mnar[,1],use = "complete.obs")
var(z.mnar[,2])


# Mean Imputation

mean.mnar<-z.mnar[,1]
mean.mnar[is.na(mean.mnar)]<-mean(z.mnar[,1], na.rm=T)
cov(mean.mnar,z.mnar[,2])
var(mean.mnar, na.rm=T)
var(z.mnar[,2])

##############################################################################################
##################### Simulation binary imitating Titanic data ###############################
##############################################################################################

#### Simulation binary imitating Titanic data

# Obtaining the parameters to generate the data based on the case study

titanic.simu<-titanic.omit 
titanic.simu$pclassbinary<-(titanic.omit$pclass=='1st')*1
titanic.surv.simu<-glm(survived ~ pclassbinary + sex + age, family=binomial, data = titanic.simu)
titanic.class.simu<-glm(pclassbinary ~ sex + age, family=binomial, data = titanic.simu)
mean(log(titanic.simu$age))
var(log(titanic.simu$age))
hist(log(titanic.simu$age))
plot(density(log(titanic.simu$age)))

# Generating the data

B <- 2500 # number of simulation replications
coefs.complete <- matrix(0, B, 4)
coefs.incomplete <- matrix(0, B, 4)
for (i in 1:B) {

    logage <- rnorm(1000,3.26,0.568)
    age<-exp(logage)
    sex <- rbinom(1000, 1, prob=0.5)

    p.class<- 1/(1+exp(-(-2.68446-0.52206*sex+0.07288*age)))
    class <- rbinom(1000, 1, prob=p.class)

    p.sur<- 1/(1+exp(-(2.18092+1.93246*class-3.03775*sex-0.04085*age)))
    survived <- rbinom(1000, 1, prob=p.sur)

    data.complete<-data.frame(survived,class,sex,age)

    p.missp.age<- 1/(1+exp(-(2.11-1.5*class-2.85*survived)))
    # p.missp.age<- 1/(1+exp(-(-5.11+0.2*age)))
    miss.age <- rbinom(1000, 1, prob=p.missp.age)

    data.incomplete<-data.frame(survived,class,sex,age)
    data.incomplete$age[miss.age==1]<-NA
    #sum(is.na(data.incomplete$age))
    data.incomplete<-na.omit(data.incomplete) 


    coefs.complete[i, ] <- coef(glm(survived ~ class + sex + age, data=data.complete,family=binomial))
    coefs.incomplete[i, ] <- coef(glm(survived ~ class + sex + age, data=data.incomplete,family=binomial))
}


par(mfrow=c(1,2))

hist.complete <- hist(coefs.complete[,2], plot = FALSE, freq=FALSE)
plot(hist.complete, border = "dark blue", col = "light blue",
     main = expression(paste("Histogram of ", hat(beta)[1])), sub="Complete Data", ylim=c(0,500), xlim=c(1,3), xlab = expression(hat(beta)[1]))
abline(v = 1.9, col = "red")

hist.incomplete <- hist(coefs.incomplete[,2], plot = FALSE)
plot(hist.incomplete, border = "dark blue", col = "light blue",
     main = expression(paste("Histogram of ", hat(beta)[1])), sub="Missing Data LD", xlim=c(0,2.5), xlab = expression(hat(beta)[1]))
abline(v = 1.9, col = "red")

##############################################################################################
################## Simulation binary imitating Titanic data IPW ##############################
##############################################################################################

#### Simulation binary imitating Titanic data IPW

B.ipw <- 2500 # number of simulation replications
coefs.complete.ipw <- matrix(0, B.ipw, 4)
coefs.incomplete.ipw <- matrix(0, B.ipw, 4)
for (i in 1:B.ipw) {

    logage.ipw <- rnorm(1000,3.26,0.568)
    age<-exp(logage.ipw)
    sex <- rbinom(1000, 1, prob=0.5)

    p.class.ipw<- 1/(1+exp(-(-2.68446-0.52206*sex+0.07288*age)))
    class <- rbinom(1000, 1, prob=p.class.ipw)

    p.sur.ipw<- 1/(1+exp(-(2.18092+1.93246*class-3.03775*sex-0.04085*age)))
    survived <- rbinom(1000, 1, prob=p.sur.ipw)

    data.complete.ipw<-data.frame(survived,class,sex,age)

    p.missp.age.ipw<- 1/(1+exp(-(2.11-1.5*class-2.85*survived)))
    # p.missp.age<- 1/(1+exp(-(-5.11+0.2*age)))
    miss.age.ipw <- rbinom(1000, 1, prob=p.missp.age.ipw)

    data.incomplete.ipw<-data.frame(survived,class,sex,age)
    data.incomplete.ipw$age[miss.age.ipw==1]<-NA
    data.incomplete.ipw$r<-as.numeric(!is.na(data.incomplete.ipw$age))
    #sum(is.na(data.incomplete$age))
    
    ipw.glm<-glm(r ~ class + survived, data=data.incomplete.ipw,family=binomial)
    data.incomplete.ipw$w<-1/fitted(ipw.glm)

    coefs.complete.ipw[i, ] <- coef(glm(survived ~ class + sex + age, data=data.complete.ipw,family=binomial))
    coefs.incomplete.ipw[i, ] <- coef(glm(survived ~ class + sex + age, data=data.incomplete.ipw, weights=data.incomplete.ipw$w, family=binomial))
}


par(mfrow=c(1,2))

hist.complete.ipw <- hist(coefs.complete.ipw[,2], plot = FALSE, freq=FALSE)
plot(hist.complete.ipw, border = "dark blue", col = "light blue",
     main = expression(paste("Histogram of ", hat(beta)[1])), sub="Complete Data", ylim=c(0,500), xlim=c(1,3), xlab = expression(hat(beta)[1]))
abline(v = 1.9, col = "red")

hist.incomplete.ipw <- hist(coefs.incomplete.ipw[,2], plot = FALSE)
plot(hist.incomplete.ipw, border = "dark blue", col = "light blue",
     main = expression(paste("Histogram of ", hat(beta)[1])), sub="IPW", xlim=c(0.5,3.5), xlab = expression(hat(beta)[1]))
abline(v = 1.9, col = "red")


##############################################################################################
############################# Simulation multiple imputation #################################
##############################################################################################

###### Simulation multiple imputation


B.mi <- 500 # number of simulation replications
coefs.complete.mi <- matrix(0, B.mi, 4)
coefs.incomplete.mi <- matrix(0, B.mi, 4)
for (i in 1:B.mi) {

    logage <- rnorm(1000,3.26,0.568)
    age<-exp(logage)
    sex <- rbinom(1000, 1, prob=0.5)

    p.class<- 1/(1+exp(-(-2.68446-0.52206*sex+0.07288*age)))
    class <- rbinom(1000, 1, prob=p.class)

    p.sur<- 1/(1+exp(-(2.18092+1.93246*class-3.03775*sex-0.04085*age)))
    survived <- rbinom(1000, 1, prob=p.sur)

    data.complete.mi<-data.frame(survived,class,sex,age)

    p.missp.age<- 1/(1+exp(-(2.11-1.5*class-2.85*survived)))
    # p.missp.age<- 1/(1+exp(-(-1.11-1.09*sex-1.85*class)))
    miss.age <- rbinom(1000, 1, prob=p.missp.age)

    data.incomplete.mi<-data.frame(survived,class,sex,age)
    data.incomplete.mi$age[miss.age==1]<-NA

    coefs.complete.mi[i, ] <- coef(glm(survived ~ class + sex + age, data=data.complete.mi,family=binomial))
    
    imputation <- mice(data.incomplete.mi, m=5)
    allimplogreg  <- with(data=imputation, exp=glm(survived ~ class + sex + age, family=binomial))
    allimplogreg.pooled <- pool(allimplogreg)

    coefs.incomplete.mi[i, ]<-allimplogreg.pooled$qbar
    
}


par(mfrow=c(1,2))

hist.complete.mi <- hist(coefs.complete.mi[,2], plot = FALSE, freq=FALSE)
plot(hist.complete.mi, border = "dark blue", col = "light blue",
     main = expression(paste("Histogram of ", hat(beta)[1])), sub="Complete Cases", ylim=c(0,130), xlim=c(1,3), xlab = expression(hat(beta)[1]))
abline(v = 1.9, col = "red")

hist.incomplete.mi <- hist(coefs.incomplete.mi[,2], plot = FALSE)
plot(hist.incomplete.mi, border = "dark blue", col = "light blue",
     main = expression(paste("Histogram of ", hat(beta)[1])), sub="MI", ylim=c(0,150), xlim=c(1,3), xlab = expression(hat(beta)[1]))
abline(v = 1.9, col = "red")



