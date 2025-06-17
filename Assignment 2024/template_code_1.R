## Multilevel Models: Longitudinal

## Girls = 2, Boys = 1

# Home


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

install.packages("merTools")
install.packages("pbkrtest")
 
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
#library(Rcmdr) # This one opens a windows

library(lme4)
library(lattice)
library("lsmeans")
library(arm)
library(car)
library(pbkrtest)
library(LMERConvenienceFunctions)
library(coda)
library(languAGER)
library(Hmisc)

require(tigerstats)
require(merTools)
require(pbkrtest)

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
############################## Growth Sizes ###################################
#############################################################################################


#########################
############################## Analysis of the data.
#########################

## Reading the data

growth.int1 <- read.table(file="growth.txt", header=TRUE)
rownames(growth.int1)=NULL
growth.int1$AGE0<-growth.int1$AGE-8
growth.int1$SEX<-growth.int1$SEX-1
growth.int1 <- growth.int1[, -which(names(growth.int1) %in% c("INDIV"))]
head(growth.int1)

## Attach data to the search path

attach(growth.int1)

## Spaghettiplot

n=length(unique(IDNR))

interaction.plot(AGE,IDNR,MEASURE, xlab="AGE in years", ylab="Measure", legend=F) 

# Plot individual profiles + mean with function

Spaghetti.Plot.new(Dataset=growth.int1, Outcome=MEASURE, Time=AGE, Id=IDNR, tto=SEX,
			  Add.Profiles = TRUE, Add.Mean = TRUE, Add.Median = FALSE, Col = 8, 
			  colme=c(2,1), ltyme=c(1,1), Lwd.Me = 3, xlim=c(1,2), ylim=c(60,150), bytto=TRUE) 

legend("topright", inset=.05, ,legend = c("Control","Intervertion"),horiz=TRUE
         , text.col = c("black","red")
         , pt.bg = c("black","red")
         , pch = c(16,16),col = c("black","red"))

## Adding a loess curve per group. There is a problem probably due to small n
with (growth.int1, {
lines (loess.smooth (AGE[SEX==0], MEASURE[SEX==0], family = "gaussian"),lty = 3,col=4,lwd = 4)
lines (loess.smooth (AGE[SEX==1], MEASURE[SEX==1], family = "gaussian"),lty = 4,col=5,lwd = 4)
})

# Plot individual profiles + mean  by hand

plot(AGE, MEASURE, type = "n")

for (i in 1:length(unique(IDNR)))
{ lines(y =growth.int1$MEASURE[growth.int1$IDNR == unique(growth.int1$IDNR)[i]], 
                x =growth.int1$AGE[growth.int1$IDNR == unique(growth.int1$IDNR)[i]], 
                col = 8)
}

mean.prog0=tapply(MEASURE[SEX==0], INDEX =AGE[SEX==0], FUN = mean, na.rm = TRUE)
lines(mean.prog0, x = unique(AGE), lwd = 3)

mean.prog1=tapply(MEASURE[SEX==1], INDEX =AGE[SEX==1], FUN = mean, na.rm = TRUE)
lines(mean.prog1, x = unique(AGE), lwd = 4,col=2)


## Descriptives

## Mean:
growth.mean=tapply(MEASURE,list(AGE,SEX),mean)

## Standard deviation:
growth.sd=tapply(MEASURE,list(AGE,SEX),sd)

## Variance:
growth.var=tapply(MEASURE,list(AGE,SEX),var)

## Frequency:
growth.n=table(AGE,SEX)

## Boxplots:

boxplot(MEASURE~AGE,xlab="AGE (in years)",ylab="Measure", col = "gray")

## Boxplots per SEX

par(mfrow=c(2,1))
boxplot(MEASURE[SEX==0]~AGE[SEX==0],main="No intervention",xlab="AGE (in years)",ylab="Measure", col = "gray")
boxplot(MEASURE[SEX==1]~AGE[SEX==1],main="Intervention",xlab="AGE (in years)",ylab="Measure", col = "gray")

## Plotting mean evolutions

par(mfrow=c(1,1))
plot(AGE[IDNR==1],growth.mean[,1],type="b",xlim=c(8,14), ## SEX=1
ylim=c(10,30),xlab="AGE (months?)",ylab="Distance (mm?)",axes=F,
main="Mean evolution (with 1 SE intervals)")
axis(side=1,at=c(8,10,12,14),labels=c(8,10,12,14))
axis(side=2,at=seq(10,30,5))

box()
points(AGE[IDNR==1],growth.mean[,2],type="b",col="red") ## SEX=2
errbar(AGE[IDNR==1]-.005,growth.mean[,1],growth.sd[,1],.1)
errbar(AGE[IDNR==1]+.005,growth.mean[,2],growth.sd[,2],.1,col="red")

legend("topright", inset=.05, ,legend = c("Boys","Girls"),horiz=TRUE
         , text.col = c("black","red")
         , pt.bg = c("black","red")
         , pch = c(16,16),col = c("black","red"))

## Correlations

## Reshaping the data into a wide form

growth.int2 <- reshape(growth.int1[,c(1,2,3,4)], 
  timevar = "AGE", idvar = c("IDNR", "SEX"), direction = "wide")
tail(growth.int2)

## Correlation between the Measure scores  at  different AGEs

cor(growth.int2[,3:6])
cor(growth.int2[growth.int2$SEX==1,3:6])
cor(growth.int2[growth.int2$SEX==0,3:6])

## Linear regression per person + histograms

## Displaying the linear regression per person:

cf<-sapply(growth.int1$IDNR, function(x) 
    coef(lm(MEASURE~AGE0, data=subset(growth.int1, IDNR==x))))

plot(cf[1,],cf[2,],xlab="Individual regression intercept",
     ylab="Individual regression slope", main="Individual regression intercept versus slope")
identify(cf[1,],cf[2,],n=1)

Sx<-reorder(growth.int1$IDNR, cf[1,])

xyplot(MEASURE ~ AGE0|Sx,groups=SEX,data=growth.int1,
 type=c('p','r'),auto.key=T,aspect="xy",
 par.settings=list(axis.text=list(cex=0.6),
 fontsize=list(text=8, points=10)),
 scales=list(
    x=list(
      at=c(0,2,4,6),
      labels=c("0","2","4","6")))
)


## Linear regression per participant of MEASURE on AGE

## Coefficients

lin.reg.coef <- by(growth.int1, growth.int1$IDNR, 
          function(data) coef(lm(MEASURE ~ AGE0, data=data)))
lin.reg.coef1 <- unlist(lin.reg.coef)
names(lin.reg.coef1) <- NULL 
lin.reg.coef2=matrix(lin.reg.coef1,length(lin.reg.coef1)/2,2,byrow = TRUE)

## R squared

lin.reg.r.squared <- by(growth.int1, growth.int1$IDNR, 
          function(data) summary(lm(MEASURE ~ AGE, data=data))$r.squared )
lin.reg.r.squared1<- as.vector(unlist(lin.reg.r.squared))

## Histograms

par(mfrow=c(3,1))
hist(lin.reg.coef2[,1],xlab="Intercept",col="lightblue",main="Histogram of individual intercepts")
hist(lin.reg.coef2[,2],xlab="Slope",col="lightblue",main="Histogram of individual slopes")
hist(lin.reg.r.squared1,xlab="R squared",col="lightblue",main="Histogram of individual R squared")

##   Correlations

par(mfrow=c(1,1))
int.slope.corr=cor(lin.reg.coef2[,1],lin.reg.coef2[,2])
plot(lin.reg.coef2[,1],lin.reg.coef2[,2],xlab="Intercept", ylab="Slope", main="Intercept versus Slope")

## Plotting individual regression lines per group

reg.coef=cbind(lin.reg.coef2, growth.int1[growth.int1$AGE==8,]$SEX)

mean.int<-tapply(reg.coef[,1],reg.coef[,3],mean)
mean.slope<-tapply(reg.coef[,2],reg.coef[,3],mean)

x_lim <- c(0,6)
y_lim <- c(15,35)
par(mfrow=c(1,2))
plot(growth.int1$AGE0,growth.int1$MEASURE,type="n",xlim=x_lim,ylim=y_lim,main="No intervention",xlab="AGE-1 (in years)",ylab="Measure",axes=F)
axis(side=1,at=c(0,2,4,6),labels=c(0,2,4,6))
axis(side=2,at=seq(10,30,5))
box()
for (i in 1:n)
{if (reg.coef[i,3]==0) 
{curve(cbind(1,x)%*%reg.coef[i,1:2],add=T,col="gray")}}
curve(cbind(1,x)%*%c(mean.int[1],mean.slope[1]),add=T,lwd=2)

plot(growth.int1$AGE0,growth.int1$MEASURE,type="n",xlim=x_lim,ylim=y_lim,main="Intervention",xlab="AGE-1 (in years)",ylab="Measure",axes=F)
axis(side=1,at=c(0,2,4,6),labels=c(0,2,4,6))
axis(side=2,at=seq(10,30,5))
box()
for (i in 1:n)
{if (reg.coef[i,3]==1) 
{curve(cbind(1,x)%*%reg.coef[i,1:2],add=T,col="gray")}}
curve(cbind(1,x)%*%c(mean.int[2],mean.slope[2]),add=T,lwd=2)

## ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML 
## Fitting the model with ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML 
## ML but with LMER instead of LME

## Different random intercept and slope

growth.lmer1<-lmer(MEASURE~1+AGE0*SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1)
mcp.fnc(growth.lmer1)
help(mcp.fnc)

summary(growth.lmer1)
display(growth.lmer1)
anova(growth.lmer1)

##################### FIXED EFFETCS?

## Estimating the fixed effects via bootstrap

fixed.boot=bootMer(growth.lmer1,  fixef, use.u = TRUE, nsim = 250)
fixed.boot
help(bootMer)
help(fixef)
summary(fixed.boot)

#help(pvalues)

## Calculating confidence intervals for the fixed effects via bootstrap

# Another way to calculate profile likelihood CI
growth.prof=profile(growth.lmer1,signames=FALSE)
confint(growth.prof)
warnings()

confint(growth.lmer1, level = 0.95,method="profile",oldNames = FALSE)
confint(growth.lmer1,par=5:8,method="Wald",oldNames = FALSE)
confint(growth.lmer1,method="boot",boot.type ="perc",oldNames = FALSE,nsim=1000)
confint(growth.lmer1,method="boot",boot.type ="basic",oldNames = FALSE,nsim=1000)

## Get the KR-approximated degrees of freedom

growth.lmer1.df.KR <- get_Lb_ddf(growth.lmer1, fixef(growth.lmer1))

## Get p-values from the t-distribution using the t-values and approximated
## degrees of freedom

growth.lmer1.coef=coef(summary(growth.lmer1))
growth.lmer1.p.KR <- cbind(growth.lmer1.coef,df=growth.lmer1.df.KR,2 * (1 - pt(abs(growth.lmer1.coef[,3]), growth.lmer1.df.KR)))
growth.lmer1.p.KR

## Another way to get the p-values require(lmerTest) and refit the model

require(lmerTest)
growth.lmer1<-lmer(MEASURE~1+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1)
summary(growth.lmer1)
anova(growth.lmer1)
anova(growth.lmer1, type=1)
help(anova)

anova(lmer(MEASURE~1+SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))
anova(lmer(MEASURE~1+SEX+AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))
anova(lmer(MEASURE~1+SEX+AGE0+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))
anova(lmer(MEASURE~1+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))
anova(lmer(MEASURE~1+AGE0*SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))

anova(lmer(MEASURE~1+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1), type=1)
anova(lmer(MEASURE~1+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1), type=2)
anova(lmer(MEASURE~1+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1), type=3)

summary(lmer(MEASURE~1+AGE0*SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))

anova(lmer(MEASURE~1+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1),
      lmer(MEASURE~1+SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1), 
      lmer(MEASURE~1+SEX+AGE+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1),
      lmer(MEASURE~1+SEX+AGE+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))
anova(lmer(MEASURE~1+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1),
      lmer(MEASURE~1+AGE+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1), 
      lmer(MEASURE~1+AGE+SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1),
      lmer(MEASURE~1+AGE+SEX+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))
anova(lmer(MEASURE~1+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1),
      lmer(MEASURE~1+SEX*AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1),
      lmer(MEASURE~1+SEX*AGE0+SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1),
      lmer(MEASURE~1+SEX*AGE0+AGE0+SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1))

## ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML 
## Fitting the model with lme  ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML ML 
## ML but with LME instead of LMER (p-values)

help(lme)
growth.lme1<-lme(MEASURE~1+AGE0*SEX,random=~1+AGE0|IDNR, method = "ML", data=growth.int1)
summary(growth.lme1)
VarCorr(growth.lme1)
getVarCov(growth.lme1)


## Likelihood ratio tests

growth.lmer1.boy<-lmer(MEASURE~1+AGE0+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1)
growth.lmer1.girl<-lmer(MEASURE~1+AGE0+SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1)
anova(growth.lmer1.boy,growth.lmer1.girl,growth.lmer1)

## Final model

growth.lmer1.slopesex <- lmer(MEASURE~1+AGE0+AGE0:SEX+(1 + AGE0|IDNR), REML = FALSE, data=growth.int1)

summary(growth.lmer1.slopesex)

## Random effects covariance matrix

D.growth=unclass(VarCorr(growth.lmer1))$IDNR
D.growth

## CI lme 

intervals(growth.lme1)
warnings()

## Plotting the profiles

xyplot(MEASURE~AGE0|Sx, groups=SEX, data=growth.int1, fit=growth.lmer1,
strip = FALSE, aspect="xy", pch=16, grid=TRUE,
panel = function(x, y, fit, subscripts) {
panel.xyplot(x, y)
ypred <- predict(growth.lmer1)[subscripts]
panel.lines(x, ypred, col = "black")
},
xlab = "AGE - 8", ylab = "Measure")

## The fitted and predict function produce the same result

cbind(growth.int1,predict(growth.lmer1),fitted(growth.lmer1))


## Predicted random effects

growth.lmer1.re=ranef(growth.lmer1.slopesex)$IDNR
head(growth.lmer1.re,10)
str(growth.lmer1.re)
table(growth.int1$SEX)

## The residual from lmer and lme are basically the same

cbind(ranef(growth.lmer1)$IDNR,ranef(growth.lme1))

## Ploting the randon intercept and slope for the entire data set

dotplot(ranef(growth.lmer1, condVar = T))

## Ploting the randon intercept and slope for a subset of the data
## In this case the first 20 individuals given in s

r.int <- ranef(growth.lmer1, condVar=TRUE)

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

plot(growth.lmer1.re[1:58,], main=bquote("Random intercept"~b[i0]~"versus random slope"~ b[i1]~"SEX=1"), xlab=expression(b[i0]), ylab="")
mtext(expression(b[i1]),  side = 2, line = 3, las = 1)
identify(growth.lmer1.re,n=3) 

# Random intercept (b0i) versus random slope (b1i)


## Creating the subject specific intercepts and slopes
## Here we have to use the model fitted with lme

ind.coef=coef(growth.lme1)
head(ind.coef)
prog=growth.int1[growth.int1$AGE0==0,]$SEX

int.subject=ind.coef[,1]+ind.coef[,3]*prog
slope.subject=ind.coef[,2]+ind.coef[,4]*prog
plot(int.subject[1:58],slope.subject[1:58], xlab=expression(pi[i0]), ylab="",
main="Random intercept versus random slope \n(Including the fixed effects) \n SEX=1")
identify(int.subject,slope.subject,n=3)
mtext(expression(pi[i1]),  side = 2, line = 3, las = 1) 


## Comparing the intercepts and slopes obtained from the subject specific regressions and
## the mixed effects models

## Creating a data based with the  coefficients from both analysis

subj.model.lmer1=data.frame(intsub=reg.coef[,1],slopesub=reg.coef[,2],intmodel=int.subject,slopemodel=slope.subject,SEX=reg.coef[,3])
head(subj.model.lmer1)

## Selecting SEX==0

subj.model.lmer10=subj.model.lmer1[subj.model.lmer1$SEX==0,]

# Making the plot

par(mfrow=c(1,1))
plot(subj.model.lmer10$intsub,subj.model.lmer10$slopesub,type="n",main="Per subject OLS versus LMM estimates \n of the slope and intercepts (SEX=0)",
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

population.int0=fixef(growth.lmer1)[1]
population.slop0=fixef(growth.lmer1)[2]

points(population.int0,population.slop0, col = "red", pch=15,cex=2)

## Selecting SEX==1. MUCH BETTER!

subj.model.lmer11=subj.model.lmer1[subj.model.lmer1$SEX==1,]

plot(subj.model.lmer11$intsub,subj.model.lmer11$slopesub, xlab=expression(pi[i0]), ylab="",
main="OLS intercept versus OLS slope \n SEX=1")
identify(subj.model.lmer11$intsub,subj.model.lmer11$slopesub,n=3)
mtext(expression(pi[i1]),  side = 2, line = 3, las = 1) 

## Making the plot

plot(subj.model.lmer11$intsub,subj.model.lmer11$slopesub,type="n",main="Per subject OLS versus LMM estimates \n of the slope and intercepts (SEX=1)",
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

population.int1=fixef(growth.lmer1)[1]+fixef(growth.lmer1)[3]
population.slop1=fixef(growth.lmer1)[2]+fixef(growth.lmer1)[4]

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

population.int1=fixef(growth.lme1)[1]+fixef(growth.lme1)[3]
population.slop1=fixef(growth.lme1)[2]+fixef(growth.lme1)[4]

points(population.int1,population.slop1, col = "red", pch=18,cex=2)

legend("topright", inset=.05, ,legend = c("OLS","LMM","Population"),horiz=TRUE
         , text.col = c("red", "blue","red")
         , pt.bg = c("red","blue","red")
         , pch = c(16,15,18),col = c("red", "blue","red"))



## Displaying the linear regression per person:

cf<-sapply(growth.int1$IDNR, function(x) 
    coef(lm(MEASURE~AGE0, data=subset(growth.int1, IDNR==x))))

Sx2<-reorder(subj.model.lmer1, subj.model.lmer1$intsub)

xyplot(MEASURE ~ AGE0|Sx,groups=SEX,data=growth.int1,
 type=c('p','r'),auto.key=T,aspect="xy",
 par.settings=list(axis.text=list(cex=0.6),
 fontsize=list(text=8, points=10)),
 scales=list(
    x=list(
      at=c(0,2,4,6),
      labels=c("0","2","4","6"))),
panel = function(x, y, ...) {
  panel.abline(a = 0, b = 2)}
)

pred.re=ranef(growth.lmer1,condVar = TRUE)
dotplot(pred.re, lattice.options=list(layout=c(1,2)))
qqmath(pred.re, lattice.options=list(layout=c(1,2)))


RLRsim(growth.lmer1)

detach(growth.int1)

