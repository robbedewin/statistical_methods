## Multilevel Models: Longitudinal

# Home

setwd("C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Bioinformatics\\Multilevel-Models\\Data-R-code")

update.packages()

# Packages

install.packages("foreign")
install.packages("UsingR")
install.packages("gplots")
install.packages("xtable")
install.packages("Rcmdr")
install.packages("stats")
install.packages("asbio")
install.packages("mvtnorm")
install.packages("multcomp")
install.packages("lawstat")
install.packages("gmodels")
install.packages("doBy")
install.packages("faraway")
install.packages("Hmisc")
install.packages("RLRsim")



install.packages("corrplot") 

install.packages("lme4") 
install.packages("nlme")
install.packages("arm")
install.packages("pbkrtest")
install.packages("LMERConvenienceFunctions")
install.packages("languageR")
install.packages("CorrMixed")
install.packages("lsmeans")
install.packages("lmerTest")
install.packages("varTestnlme")
 
## Practical Regression and Anova using R by Julian J. Faraway

# install.packages("faraway")

library('corrplot')
library(CorrMixed)

library(nlme)
library(doBy)
library(foreign)
library(UsingR)
library(gplots)
library(xtable)
library(stats)
library("asbio")
library(graphics)
library(mvtnorm)
library(matrixcalc)
library(multcomp)
library(lawstat)
library(MASS)
library(faraway) # library of the book of linearodels in r. You have this book.
library(Rcmdr) # This one opens a windows

library(lme4)
library(lattice)
library("lsmeans")
library(arm)
library(car)
library(pbkrtest)
library(LMERConvenienceFunctions)
library(coda)
library(languageR)
library(Hmisc)

##############################################################
#           General function to plot error bars              #
##############################################################

errbar=function(x,y,height,width,lty=1,col="black")
{arrows(x,y,x,y+height,angle=90,length=width,lty=lty,
col=col)
arrows(x,y,x,y-height,angle=90,length=width,lty=lty,
col=col)}

#############################################################################
# General to create an Spaghetti Plot with the mean o median function on it #
#############################################################################

## tto is the treatment variable or any other variable defining the  groups
## if bytto is true then the means or medians are calculated for  every group
## colme is a  vector given the colors for every mean/median curve

Spaghetti.Plot.new=
function (Dataset, Outcome, Time, Id, tto, Add.Profiles = TRUE, Add.Mean = TRUE, 
    Add.Median = FALSE, Col = 8, colme, ltyme, Lwd.Me = 3, xlim, ylim, bytto=FALSE,...) 
{
    Object <- x <- Dataset
    Outcome <- x[, paste(substitute(Outcome))]
    Time <- x[, paste(substitute(Time))]
    Id <- x[, paste(substitute(Id))]
    tto<-x[, paste(substitute(tto))]
    Data <- data.frame(cbind(Outcome, Time, Id,tto))
    ttovalue=unique(tto)
    max_val_y <- max(Data$Outcome, na.rm = TRUE)
    max_val_x <- max(Data$Time, na.rm = TRUE)
    if (missing(xlim) == TRUE) {
        xlim = c(0, max_val_x)
    }
    if (missing(ylim) == TRUE) {
        ylim = c(0, max_val_y)
    }
    plot(y = Outcome, x = Time, type = "n", ylim = ylim, xlim = xlim, 
        ...)
    if (Add.Profiles == TRUE) {
        for (i in 1:length(unique(Id))) {
            lines(y = Data$Outcome[Data$Id == unique(Data$Id)[i]], 
                x = Data$Time[Data$Id == unique(Data$Id)[i]], 
                col = Col)
        }
    }
    if (Add.Mean == TRUE) 
	{
    	 if (bytto==TRUE)
		{
	       #mean.tto=rep(0,length(ttovalue))
    		 for (i in 1:length(ttovalue)) {
		 mean.tto <- tapply(Data$Outcome[Data$tto==ttovalue[i]], INDEX = Data$Time[Data$tto==ttovalue[i]], FUN = mean, 
                               na.rm = TRUE)
            lines(mean.tto, x = unique(Data$Time), lwd = Lwd.Me, col=colme[i], lty=ltyme[i])
		}}else
		  {
    		   mean <- tapply(Data$Outcome, INDEX = Data$Time, FUN = mean, 
                       na.rm = TRUE)
              lines(mean, x = unique(Data$Time), lwd = Lwd.Me, lty=ltyme[i])
		  }
    }
    if (Add.Median == TRUE) 
	{
    	 if (bytto==TRUE)
		{
	       #median.tto=rep(0,length(ttovalue))
    		 for (i in 1:length(ttovalue)) {
		 median.tto <- tapply(Data$Outcome[Data$tto==ttovalue[i]], INDEX = Data$Time[Data$tto==ttovalue[i]], FUN = median, 
                                 na.rm = TRUE)
            lines(median.tto, x = unique(Data$Time), lwd = Lwd.Me,col=colme[i], lty=ltyme[i])
		}}else
		  {
    		   median<- tapply(Data$Outcome, INDEX = Data$Time, FUN = median, 
                       na.rm = TRUE)
              lines(median, x = unique(Data$Time), lwd = Lwd.Me, lty=ltyme[i])
		  }
    }
}


########################################################################
######################### Examples #####################################
########################################################################

#############################################################################################
############################## Early dietary intervention ###################################
#############################################################################################


#########################
############################## Analysis of the data.
#########################

## Reading the data

early.int1 <- read.table("earlyint.txt", header=T, sep=",")
rownames(early.int1)=NULL
early.int1$age0<-early.int1$age-1
head(early.int1)

early.int1.table <- xtable(early.int1[1:24,])
print(early.int1.table)

## Attach data to the search path

attach(early.int1)

## Spaghettiplot

n=length(unique(id))
interaction.plot(age,id,cog, xlab="Agein years", ylab="IQ", legend=F) 

# Plot individual profiles + mean with function

Spaghetti.Plot.new(Dataset=early.int1, Outcome=cog, Time=age, Id=id, tto=program,
			  Add.Profiles = TRUE, Add.Mean = TRUE, Add.Median = FALSE, Col = 8, 
			  colme=c(2,1), ltyme=c(1,1), Lwd.Me = 3, xlim=c(1,2), ylim=c(60,150), bytto=TRUE) 

legend("topright", inset=.05, ,legend = c("Control","Intervertion"),horiz=TRUE
         , text.col = c("black","red")
         , pt.bg = c("black","red")
         , pch = c(16,16),col = c("black","red"))

## Adding a loess curve per group. There is a problem probably due to small n
with (early.int1, {
lines (loess.smooth (age[program==0], cog[program==0], family = "gaussian"),lty = 3,col=4,lwd = 4)
lines (loess.smooth (age[program==1], cog[program==1], family = "gaussian"),lty = 4,col=5,lwd = 4)
})

# Plot individual profiles + mean  by hand

plot(age, cog, type = "n")

for (i in 1:length(unique(id)))
{ lines(y =early.int1$cog[early.int1$id == unique(early.int1$id)[i]], 
                x =early.int1$age[early.int1$id == unique(early.int1$id)[i]], 
                col = 8)
}

mean.prog0=tapply(cog[program==0], INDEX =age[program==0], FUN = mean, na.rm = TRUE)
lines(mean.prog0, x = unique(age), lwd = 3)

mean.prog1=tapply(cog[program==1], INDEX =age[program==1], FUN = mean, na.rm = TRUE)
lines(mean.prog1, x = unique(age), lwd = 4,col=2)


## Descriptives

## Mean:
early.mean=tapply(cog,list(age,program),mean)

## Standard deviation:
early.sd=tapply(cog,list(age,program),sd)

## Variance:
early.var=tapply(cog,list(age,program),var)

## Frequency:
early.n=table(age,program)

## Boxplots:

boxplot(cog~age,xlab="Age (in years)",ylab="IQ", col = "gray")

## Boxplots per program

par(mfrow=c(2,1))
boxplot(cog[program==0]~age[program==0],main="No intervention",main="No intervention",xlab="Age (in years)",ylab="IQ", col = "gray")
boxplot(cog[program==1]~age[program==1],main="Intervention",main="No intervention",xlab="Age (in years)",ylab="IQ", col = "gray")

## Plotting mean evolutions

plot(age[id==1],early.mean[,1],type="b",xlim=c(1,2),
ylim=c(40,160),xlab="Age (in years)",ylab="IQ",axes=F,
main="Mean evolution (with 1 SE intervals)")
axis(side=1,at=c(1,1.5,2),labels=c(1,1.5,2))
axis(side=2,at=seq(40,160,20))

box()
points(age[id==1],early.mean[,2],type="b",col="red")
errbar(age[id==1]-.005,early.mean[,1],early.sd[,1],.1)
errbar(age[id==1]+.005,early.mean[,2],early.sd[,2],.1,col="red")

legend("topright", inset=.05, ,legend = c("Control","Intervertion"),horiz=TRUE
         , text.col = c("black","red")
         , pt.bg = c("black","red")
         , pch = c(16,16),col = c("black","red"))

## Correlations

## Reshaping the data into a wide form

early.int2 <- reshape(early.int1[,c(1,2,3,4)], 
  timevar = "age", idvar = c("id", "program"), direction = "wide")

## Correlation between the IQ scores  at  different ages

cor(early.int2[,3:5])
cor(early.int2[early.int2$program==1,3:5])
cor(early.int2[early.int2$program==0,3:5])

## Linear regression per person + histograms

## Displaying the linear regression per person:

cf<-sapply(early.int1$id, function(x) 
    coef(lm(cog~age0, data=subset(early.int1, id==x))))

plot(cf[1,],cf[2,],xlab="Individual regression intercept",
     ylab="Individual regression slope", main="Individual regression intercept versus slope")
identify(cf[1,],cf[2,],n=1)

Sx<-reorder(early.int1$id, cf[1,])

xyplot(cog ~ age0|Sx,groups=program,data=early.int1,
 type=c('p','r'),auto.key=T,aspect="xy",
 par.settings=list(axis.text=list(cex=0.6),
 fontsize=list(text=8, points=10)),
 scales=list(
    x=list(
      at=c(0,0.5,1),
      labels=c("0","0.5","1")))
)


## Linear regression per participant of cog on age

## Coefficients

lin.reg.coef <- by(early.int1, early.int1$id, 
          function(data) coef(lm(cog ~ age0, data=data)))
lin.reg.coef1 <- unlist(lin.reg.coef)
names(lin.reg.coef1) <- NULL 
lin.reg.coef2=matrix(lin.reg.coef1,length(lin.reg.coef1)/2,2,byrow = TRUE)

## R squared

lin.reg.r.squared <- by(early.int1, early.int1$id, 
          function(data) summary(lm(cog ~ age, data=data))$r.squared )
lin.reg.r.squared1<- as.vector(unlist(lin.reg.r.squared))

## Histograms

par(mfrow=c(3,1))
hist(lin.reg.coef2[,1],xlab="Intercept",col="lightblue",main="Histogram of individual intercepts")
hist(lin.reg.coef2[,2],xlab="Slope",col="lightblue",main="Histogram of individual slopes")
hist(lin.reg.r.squared1,xlab="R squared",col="lightblue",main="Histogram of individual R squared")

##   Correlations

int.slope.corr=cor(lin.reg.coef2[,1],lin.reg.coef2[,2])
plot(lin.reg.coef2[,1],lin.reg.coef2[,2],xlab="Intercept", ylab="Slope", main="Intercept versus Slope")

## Plotting individual regression lines per group

reg.coef=cbind(lin.reg.coef2, early.int1[early.int1$age==1,]$program)

mean.int<-tapply(reg.coef[,1],reg.coef[,3],mean)
mean.slope<-tapply(reg.coef[,2],reg.coef[,3],mean)

par(mfrow=c(1,2))
plot(early.int1$age0,early.int1$cog,type="n",xlim=c(0,1),ylim=c(40,160),main="No intervention",xlab="Age-1 (in years)",ylab="IQ",axes=F)
axis(side=1,at=c(0,0.5,1),labels=c(0,0.5,1))
axis(side=2,at=seq(40,160,20))
box()
for (i in 1:103)
{if (reg.coef[i,3]==0) 
{curve(cbind(1,x)%*%reg.coef[i,1:2],add=T,col="gray")}}
curve(cbind(1,x)%*%c(mean.int[1],mean.slope[1]),add=T,lwd=2)

plot(early.int1$age0,early.int1$cog,type="n",xlim=c(0,1),ylim=c(40,160),main="Intervention",xlab="Age-1 (in years)",ylab="IQ",axes=F)
axis(side=1,at=c(0,0.5,1),labels=c(0,0.5,1))
axis(side=2,at=seq(40,160,20))
box()
for (i in 1:103)
{if (reg.coef[i,3]==1) 
{curve(cbind(1,x)%*%reg.coef[i,1:2],add=T,col="gray")}}
curve(cbind(1,x)%*%c(mean.int[2],mean.slope[2]),add=T,lwd=2)


## Fitting the model with ML

## Different random intercept and slope

early.lmer1<-lmer(cog~1+age0*program+(1 + age0|id), REML = FALSE, data=early.int1)
mcp.fnc(early.lmer1)
help(mcp.fnc)

summary(early.lmer1)
display(early.lmer1)
anova(early.lmer1)

## Estimating the fixed effects via bootstrap

fixed.boot=bootMer(early.lmer1,  fixef, use.u = TRUE, nsim = 250)
fixed.boot
summary(fixed.boot)

#help(pvalues)

## Calculating confidence intervals for the fixed effects via bootstrap

# Another way to calculate profile likelihood CI
early.prof=profile(early.lmer1,signames=FALSE)
confint(early.prof)
warnings()

confint(early.lmer1, level = 0.95,method="profile",oldNames = FALSE)
confint(early.lmer1,par=5:8,method="Wald",oldNames = FALSE)
confint(early.lmer1,method="boot",boot.type ="perc",oldNames = FALSE,nsim=1000)
confint(early.lmer1,method="boot",boot.type ="basic",oldNames = FALSE,nsim=1000)

## Get the KR-approximated degrees of freedom

early.lmer1.df.KR <- get_Lb_ddf(early.lmer1, fixef(early.lmer1))

## Get p-values from the t-distribution using the t-values and approximated
## degrees of freedom

early.lmer1.coef=coef(summary(early.lmer1))
early.lmer1.p.KR <- cbind(early.lmer1.coef,df=early.lmer1.df.KR,2 * (1 - pt(abs(early.lmer1.coef[,3]), early.lmer1.df.KR)))
early.lmer1.p.KR

## Another way to get the p-values require(lmerTest) and refit the model

require(lmerTest)
early.lmer1<-lmer(cog~1+age0*program+(1 + age0|id), REML = FALSE, data=early.int1)
summary(early.lmer1)
anova(early.lmer1)
anova(early.lmer1, type=1)
help(anova)

## Fitting the model with lme

early.lme1<-lme(cog~1+age0*program,random=~1+age0|id, method = "ML", data=early.int1)
summary(early.lme1)
VarCorr(early.lme1)
getVarCov(early.lme1)


## Likelihood ratio tests

early.lmer1.noprog<-lmer(cog~1+age0+(1 + age0|id), REML = FALSE, data=early.int1)
early.lmer1.intprog<-lmer(cog~1+age0+program+(1 + age0|id), REML = FALSE, data=early.int1)
anova(early.lmer1.noprog,early.lmer1.intprog,early.lmer1)

Final model

early.lmer1.slopeprog<-lmer(cog~1+age0+age0:program+(1 + age0|id), REML = FALSE, data=early.int1)

summary(early.lmer1.slopeprog)

## Random effects covariance matrix

D.early=unclass(VarCorr(early.lmer1))$id
D.early

## CI lme 

intervals(early.lme1)
warnings()

## Plotting the profiles

xyplot(cog~age0|Sx, groups=program, data=early.int1, fit=early.lmer1,
strip = FALSE, aspect="xy", pch=16, grid=TRUE,
panel = function(x, y, fit, subscripts) {
panel.xyplot(x, y)
ypred <- predict(early.lmer1)[subscripts]
panel.lines(x, ypred, col = "black")
},
xlab = "Centered age", ylab = "Cognitive Ability")

## The fitted and predict function produce the same result

cbind(early.int1,predict(early.lmer1),fitted(early.lmer1))


## Predicted random effects

early.lmer1.re=ranef(early.lmer1)$id
head(early.lmer1.re,10)
str(early.lmer1.re)
table(early.int1$program)

## The residual from lmer and lme are basically the same

cbind(ranef(early.lmer1)$id,ranef(early.lme1))

## Ploting the randon intercept and slope for the entire data set

dotplot(ranef(early.lmer1, condVar = T))

## Ploting the randon intercept and slope for a subset of the data
## In this case the first 20 individuals given in s

r.int <- ranef(early.lmer1, condVar=TRUE)

# We create a vector of the desired subset row numbers.

s <-1:20

# Inside a lapply function we subset both in the list, the data.frame and importantly the 
# attributes, where the variances are hidden in an array.

r.int <- lapply(r.int, function(x) {
  s2 <- which(rownames(x) %in% s)
  x <- x[s2, ]
  attributes(x)$postVar <- attributes(x)$postVar[, , s2]
  return(x)
})

# We hack the required class label and make the plot

class(r.int) <- "ranef.mer"
dotplot(r.int)

plot(early.lmer1.re[1:58,], main=bquote("Random intercept"~b[i0]~"versus random slope"~ b[i1]~"Program=1"), xlab=expression(b[i0]), ylab="")
mtext(expression(b[i1]),  side = 2, line = 3, las = 1)
identify(early.lmer1.re,n=3) 

Random intercept (b0i) versus random slope (b1i)


## Creating the subject specific intercepts and slopes
## Here we have to use the model fitted with lme

ind.coef=coef(early.lme1)
head(ind.coef)
prog=early.int1[early.int1$age0==0,]$program

int.subject=ind.coef[,1]+ind.coef[,3]*prog
slope.subject=ind.coef[,2]+ind.coef[,4]*prog
plot(int.subject[1:58],slope.subject[1:58], xlab=expression(pi[i0]), ylab="",
main="Random intercept versus random slope \n(Including the fixed effects) \n Program=1")
identify(int.subject,slope.subject,n=3)
mtext(expression(pi[i1]),  side = 2, line = 3, las = 1) 


## Comparing the intercepts and slopes obtained from the subject specific regressions and
## the mixed effects models

## Creating a data based with the  coefficients from both analysis

subj.model.lmer1=data.frame(intsub=reg.coef[,1],slopesub=reg.coef[,2],intmodel=int.subject,slopemodel=slope.subject,program=reg.coef[,3])
head(subj.model.lmer1)

## Selecting program==0

subj.model.lmer10=subj.model.lmer1[subj.model.lmer1$program==0,]

## Making the plot

plot(subj.model.lmer10$intsub,subj.model.lmer10$slopesub,type="n",main="Per subject OLS versus LMM estimates \n of the slope and intercepts (program=0)",
     xlab="Intercept",ylab="Slope")
points(subj.model.lmer10$intsub,subj.model.lmer10$slopesub, col = "red", pch=16)
points(subj.model.lmer10$intmodel,subj.model.lmer10$slopemodel, col = "blue", pch=15)
for (i in 1:103)
{
arrows(subj.model.lmer10$intsub[i],subj.model.lmer10$slopesub[i], 
       subj.model.lmer10$intmodel[i],subj.model.lmer10$slopemodel[i],
       angle = 10,length = 0.15,code=2)
}

## Adding the population coefficients

population.int0=fixef(early.lmer1)[1]
population.slop0=fixef(early.lmer1)[2]

points(population.int0,population.slop0, col = "red", pch=15,cex=2)

## Selecting program==1. MUCH BETTER!

subj.model.lmer11=subj.model.lmer1[subj.model.lmer1$program==1,]

plot(subj.model.lmer11$intsub,subj.model.lmer11$slopesub, xlab=expression(pi[i0]), ylab="",
main="OLS intercept versus OLS slope \n Program=1")
identify(subj.model.lmer11$intsub,subj.model.lmer11$slopesub,n=3)
mtext(expression(pi[i1]),  side = 2, line = 3, las = 1) 

## Making the plot

plot(subj.model.lmer11$intsub,subj.model.lmer11$slopesub,type="n",main="Per subject OLS versus LMM estimates \n of the slope and intercepts (program=1)",
     xlab="Intercept",ylab="Slope")
points(subj.model.lmer11$intsub,subj.model.lmer11$slopesub, col = "red", pch=16)
points(subj.model.lmer11$intmodel,subj.model.lmer11$slopemodel, col = "blue", pch=15)
for (i in 1:103)
{
arrows(subj.model.lmer11$intsub[i],subj.model.lmer11$slopesub[i], 
       subj.model.lmer11$intmodel[i],subj.model.lmer11$slopemodel[i],
       angle = 10,length = 0.15,code=2)
}

## Adding the population coefficients

population.int1=fixef(early.lmer1)[1]+fixef(early.lmer1)[3]
population.slop1=fixef(early.lmer1)[2]+fixef(early.lmer1)[4]

points(population.int1,population.slop1, col = "red", pch=18,cex=2)

legend("topright", inset=.05, ,legend = c("OLS","LMM","Population"),horiz=TRUE
         , text.col = c("red", "blue","red")
         , pt.bg = c("red","blue","red")
         , pch = c(16,15,18),col = c("red", "blue","red"))

identify(subj.model.lmer11$intsub,subj.model.lmer11$slopesub,n=3)
identify(subj.model.lmer11$intmodel,subj.model.lmer11$slopemodel,n=3)

## Making the plot. Both groups together

plot(subj.model.lmer1$intsub,subj.model.lmer1$slopesub,type="n",main="Per subject OLS versus LMM estimates \n of the slope and intercepts",
     xlab="Intercept",ylab="Slope")
points(subj.model.lmer1$intsub,subj.model.lmer1$slopesub, col = "red", pch=16)
points(subj.model.lmer1$intmodel,subj.model.lmer1$slopemodel, col = "blue", pch=15)
for (i in 1:103)
{
arrows(subj.model.lmer1$intsub[i],subj.model.lmer1$slopesub[i], 
       subj.model.lmer1$intmodel[i],subj.model.lmer1$slopemodel[i],
       angle = 10,length = 0.15,code=2)
}

## Adding the population coefficients

population.int1=fixef(early.lme1)[1]+fixef(early.lme1)[3]
population.slop1=fixef(early.lme1)[2]+fixef(early.lme1)[4]

points(population.int1,population.slop1, col = "red", pch=18,cex=2)

legend("topright", inset=.05, ,legend = c("OLS","LMM","Population"),horiz=TRUE
         , text.col = c("red", "blue","red")
         , pt.bg = c("red","blue","red")
         , pch = c(16,15,18),col = c("red", "blue","red"))



## Displaying the linear regression per person:

cf<-sapply(early.int1$id, function(x) 
    coef(lm(cog~age0, data=subset(early.int1, id==x))))

Sx2<-reorder(subj.model.lmer1, subj.model.lmer1$intsub)

xyplot(cog ~ age0|Sx,groups=program,data=early.int1,
 type=c('p','r'),auto.key=T,aspect="xy",
 par.settings=list(axis.text=list(cex=0.6),
 fontsize=list(text=8, points=10)),
 scales=list(
    x=list(
      at=c(0,0.5,1),
      labels=c("0","0.5","1")))
panel = function(x, y, ...) {
  panel.abline(a = 0, b = 2)}
)

pred.re=ranef(early.lmer1,condVar = TRUE)
dotplot(pred.re, lattice.options=list(layout=c(1,2)))
qqmath(pred.re, lattice.options=list(layout=c(1,2)))


RLRsim(early.lmer1)

detach(early.int1)


#############################################################################################
##################################  The Rat Pup Example #####################################
#############################################################################################

## The data come from a study in which 30 female rats were randomly assigned to receive one of 
## three doses (high, low or control) of an experimental compound. Although 10 female rats 
## were initially assigned to receive each treatment dose, three of the female rats in the 
## high-dose group died, so there are no data for their litters. In addition, litter sizes 
## varied widely, ranging from 2 to 18 pups.
## Objective of the study: To compare the birth weights of pups from litters born to female 
## rats that received the high- and low-dose treatments to the birth weights of pups from 
## litters that received the control treatment.

## Jose Pinheiro and Doug Bates, (2000) Mixed-Effects Models in S and S-PLUS.


## Reading the data

ratpup <- read.table("rat_pup.dat", h = T)
ratpup$sex1[ratpup$sex == "Female"] <- 1
ratpup$sex1[ratpup$sex == "Male"] <- 0
head(ratpup)

attach(ratpup)

## Table describing the  data

g <- function(x)c(N=length(x),MEAN=mean(x,na.rm=TRUE),SD=sd(x,na.rm=TRUE),MIN=min(x,na.rm=TRUE),MAX=max(x,na.rm=TRUE))
summarize(ratpup$weight,by=llist(ratpup$treatment,ratpup$sex),g)
xtable(summarize(ratpup$weight,by=llist(ratpup$treatment,ratpup$sex),g))

## Comparing the distributions of birth weights for each treatment by sex combination

library(lattice)  # trellis graphics
library(grid)

bwplot(weight ~ sex|treatment, data=ratpup,aspect = 2, ylab="Birth Weights", 
    xlab="SEX",main = "Boxplots of birth weights for levels of treatment by sex")

## require(lattice)

## Comparing the distributions of birth weights for each treatment

dotplot(litterid ~ weight,group=treatment, data =ratpup, xlab="Weight", ylab="Litter",
auto.key=list(space="top", column=3, cex=.8, title="", 
                      cex.title=1, lines=FALSE, points=TRUE) )

dotplot(litterid ~ weight,group=interaction(treatment,sex), data =ratpup, xlab="Weight", ylab="Litter",
auto.key=list(space="top", column=3, cex=.8, title="", 
                      cex.title=1, lines=FALSE, points=TRUE) )

with(ratpup, interaction.plot(treatment,sex,weight))

## Fitting a homocedastic model

library(nlme)

meanfull.hom <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1, 
                    random = ~1 | litterid, ratpup, method = "REML")

meanfull.hom$apVar

summary(meanfull.hom)
anova(meanfull.hom)

VarCorr(meanfull.hom)
getVarCov(meanfull.hom)


### Extracting the re with the sd in lmer. It is not possible to get the conditional 
## variance of the re in lme
hom.lmer <- lmer(weight ~ treatment + sex1 + litsize + treatment:sex1+(1 | litterid), ratpup)
dotplot(random.effects(hom.lmer, condVar = T))

str(random.effects(hom.lmer, condVar = T))
unlist(random.effects(hom.lmer, condVar = T))
###################################################################

## Display the random effects (EBLUPs) from the model.

random.effects(meanfull.hom)

## Fitting a heterocedastic model

meanfull.het <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1, 
                    random = ~1 | litterid, ratpup, method = "REML",
                    weights = varIdent(form = ~1 | treatment))

summary(meanfull.het)
anova(meanfull.het)

## Heterocedastic versus homocedastic model

anova(meanfull.hom,meanfull.het)
anova(meanfull.het,meanfull.hom)

## High-low dose: Equal residual variance

ratpup$trtgrp[treatment=="Control"] <- 1
ratpup$trtgrp[treatment == "Low" | treatment == "High"] <- 2

meanfull.hilo <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1, 
                    random = ~1 | litterid, ratpup, method = "REML",
                    weights = varIdent(form = ~1 | trtgrp))

summary(meanfull.hilo)
anova(meanfull.hilo)

anova(meanfull.het,meanfull.hilo)

## Is there a litter effect?

meanfull.hilo.nolitter <- gls(weight ~ treatment + sex1 + litsize +
        treatment:sex1, data = ratpup, weights = varIdent(form = ~1 | trtgrp))

summary(meanfull.hilo.nolitter)
anova(meanfull.hilo.nolitter)
anova(meanfull.hilo.nolitter,meanfull.hilo)

## Modeling the mean structure

meanfull.hilo.ml <- lme(weight ~ treatment + sex1 + litsize + treatment:sex1, 
                    random = ~1 | litterid, ratpup, method = "ML",
                    weights = varIdent(form = ~1 | trtgrp))

summary(meanfull.hilo.ml)
anova(meanfull.hilo.ml)


## Simulation based (exact) restricted likelihood ratio test based on 
## simulated values from the finite sample distribution for testing 
## whether the variance of a random effect is 0 in a linear mixed model 
## with  known correlation structure of the tested random 
## effect and i.i.d. errors.

require(RLRsim)
exactRLRT(meanfull.hilo)



require(varTestnlme)
varCompTest(meanfull.hilo,meanfull.hilo.nolitter)

detach(ratpup)




###############################################################################################
###############################################################################################
###############################################################################################
###############################################################################################


