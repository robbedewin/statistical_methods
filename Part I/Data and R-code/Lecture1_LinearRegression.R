
# Work
setwd("C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Bioinformatics\\Linear-Regression\\R-code-data")

install.packages("Hmisc")
install.packages("pastecs")
install.packages("doBy")
install.packages("stats")
install.packages("raster")
install.packages("XLConnect")

library(Hmisc)
library(xtable)
library(pastecs)
library(doBy)
library(stats)
library(raster)
require(plotrix)
require(car)

library(TeachingDemos)
library(mvtnorm)
library(SDMTools)
library(animation)
library(XLConnect)
library(gdata) 

##############################################################################################
################################# Important Functions ########################################
##############################################################################################


#### Function to plot a t confidence interval in a simulation context

conf.int.ariel<-function (level = 0.95, size = 50, cl = c("red", "gray"), ...) 
{
    n = ani.options("nmax")
    d = replicate(n, rnorm(size))
    m = colMeans(d)
    z = qnorm(1 - (1 - level)/2)
    y0 = m - z/sqrt(size)
    y1 = m + z/sqrt(size)
    rg = range(c(y0, y1))
    cvr = y0 < 0 & y1 > 0
    xax = pretty(1:n)
    for (i in 1:n) {
        dev.hold()
        plot(1:n, ylim = rg, type = "n", xlab = "Samples", ylab = expression("CI: [" ~ 
            bar(x) - T[1-alpha/2] * S/sqrt(n) ~ ", " ~ bar(x) + 
            T[1-alpha/2] * S/sqrt(n) ~ "]"), xaxt = "n", ...)
        abline(h = 0, lty = 2)
        axis(1, xax[xax <= i])
        arrows(1:i, y0[1:i], 1:i, y1[1:i], length = par("din")[1]/n * 
            0.5, angle = 90, code = 3, col = cl[cvr[1:i] + 1])
        points(1:i, m[1:i], col = cl[cvr[1:i] + 1])
        legend("topright", legend = format(c(i - sum(cvr[1:i]), 
            sum(cvr[1:i])), width = nchar(n)), fill = cl, bty = "n", 
            ncol = 2)
        legend("topleft", legend = paste("coverage rate:", format(round(mean(cvr[1:i]), 
            3), nsmall = 3)), bty = "n")
        ani.pause()
    }
    CI = cbind(y0, y1)
    colnames(CI) = paste(round(c((1 - level)/2, 1 - (1 - level)/2), 
        2) * 100, "%")
    rownames(CI) = 1:n
    invisible(list(level = level, size = size, CI = CI, CR = mean(cvr)))
}

#### Function to plot a histogram with a normal curve
#### g=variable, xlable=label for the x axes often the name of the variable, bins=number of intervals

histnorm = function(g, title="", xlable="", bins=10){
    m<-mean(g)
    std<-sqrt(var(g))
    myhist <- hist(g, density=10, breaks=10, prob=TRUE,plot=FALSE, 
         xlab="x-variable", ylim=c(0, 2), 
         main="Normal curve and histogram")
    hist(g, density=10, breaks=10, prob=TRUE, 
         xlab=xlable, ylim=c(0, 1.08*max(myhist$density)), 
         main=title)
    curve(dnorm(x, mean=m, sd=std), 
          col="darkblue", lwd=2, add=TRUE, yaxt="n")}


#### Function to plot a histogram with a normal curve
#### g=variable, xlable=label for the x axes often the name of the variable, bins=number of intervals

histnorm = function(g, xlable="", bins=10){
    m<-mean(g)
    std<-sqrt(var(g))
    myhist <- hist(g, density=10, breaks=10, prob=TRUE,plot=FALSE, 
         xlab="x-variable", ylim=c(0, 2), 
         main="Normal curve and histogram")
    hist(g, density=10, breaks=10, prob=TRUE, 
         xlab=xlable, ylim=c(0, 1.08*max(myhist$density)), 
         main="Normal curve and histogram")
    curve(dnorm(x, mean=m, sd=std), 
          col="darkblue", lwd=2, add=TRUE, yaxt="n")}

#### Function to obtained a CI for  the mean based on a t distribution
#### Argument t is the variable

myci <- function(t) {
  n <- length(t) # n is the sample size
  se <- sd(t)/sqrt(n); # Find the standard error of the sample
  m <- mean(t); # Find the sample mean
  cv <- qt(0.975,df=n-1) # cv is a critical value for the t distribution. P( t > cv ) = 0.025 = P( t < -cv )
  c(m-cv*se,m+cv*se) # Return the 95% confidence interval
}

#### Function to plot curly brackets
# N determines how many points in each curve
# Tilt is the ratio between the axis in the ellipse 
#  defining the curliness of each curve
# Long is the length of the straight line in the curly brackets 
#  in units of the projection of the curly brackets in this dimension
# 2*scale is the absolute size of the projection of the curly brackets 
#  in the y dimension (when theta=0)
# xcent is the location center of the x axis of the curly brackets
# ycent is the location center of the y axis of the curly brackets
# theta is the angle (in radians) of the curly brackets orientation
# col and lwd are passed to points/grid.lines

curly <- function(N = 100, Tilt = 1, Long = 2, scale = 0.1, xcent = 0.5,
                  ycent = 0.5, theta = 0, col = 1, lwd = 1, grid = FALSE){

           ymin <- scale / Tilt
           y2 <- ymin * Long
           i <- seq(0, pi/2, length.out = N)

           x <- c(ymin * Tilt * (sin(i)-1),
                  seq(0,0, length.out = 2),
                  ymin * (Tilt * (1 - sin(rev(i)))),
                  ymin * (Tilt * (1 - sin(i))),
                  seq(0,0, length.out = 2),
                  ymin * Tilt * (sin(rev(i)) - 1))

           y <- c(-cos(i) * ymin,
                  c(0,y2),
                  y2 + (cos(rev(i))) * ymin,
                  y2 + (2 - cos(i)) * ymin,
                  c(y2 + 2 * ymin, 2 * y2 + 2 * ymin),
                  2 * y2 + 2 * ymin + cos(rev(i)) * ymin)

           x <- x + xcent
           y <- y + ycent - ymin - y2

           x1 <- cos(theta) * (x - xcent) - sin(theta) * (y - ycent) + xcent
           y1 <- cos(theta) * (y - ycent) + sin(theta) * (x - xcent) + ycent

           ##For grid library:
           if(grid){
              grid.lines(unit(x1,"npc"), unit(y1,"npc"),gp=gpar(col=col,lwd=lwd))
           }

           ##Uncomment for base graphics
           else{
              par(xpd=TRUE)
              points(x1,y1,type='l',col=col,lwd=lwd)
              par(xpd=FALSE)
           }
}


##############################################################################################
##############################################################################################
##############################################################################################

######### Kalama Study

## In het kader van een onderzoek naar de fysieke ontwikkeling van kinderen heeft een 
## gezondheidswetenschapper de leeftijd (maanden) en de lichaamslengte (cm) van 12 kinderen 
## van Kalama in Egypte gemeten.

## The data is generated one time and save it does not have to be generated again

## age<-18:29
## height<-c(76.1,77,78.1,78.2,78.8,79.7,79.9,81.1,81.2,81.8,82.8,83.5)
## kalama<-data.frame(age,height)
## write.table(kalama, file = "kalama.txt", sep = "\t")


## Reading the data

kalama=read.table("kalama.txt", header=T)


## Kalama data

kalama.table.t<- xtable(t(kalama))
print(kalama.table.t)

kalama.table<- xtable(kalama)
print(kalama.table)

sum(height)
sum((height-mean(height))^2)/11


## Descriptive Statistics

options(digits=2)
descrip.kalama<-stat.desc(kalama[,c("age","height")],basic=TRUE, desc=TRUE)
descript.kalama.table<- xtable(descrip.kalama)

## Calculating the covariance and correlatio

cov.age.height<-cov(kalama$age,kalama$height)
corr.age.height<-cor(kalama$age,kalama$height)

## Testing if the population correlation is zero

corr.age.height.test= cor.test(kalama$age, kalama$height, alternative="two.sided", 
                      method = "pearson")
cor.test(~ age+height, data=kalama, alternative="two.sided", method = "pearson")


## Box plots for heigth and age

par(mfrow = c(1, 2)) 
boxplot(kalama$age,main="Age",ylab="age", col = "gray")
boxplot(kalama$height,main="Height",ylab="height", col = "gray")

## Histogram for heigth

colors = c("yellow", "red", "green", "violet") 
hist(kalama$height, right=FALSE, freq = TRUE, labels = TRUE, col=colors, main="Histogram for height", xlab="Height in cm")                 # x-axis label 
histnorm(kalama$height , xlable="Height")

stem(kalama$height,scale=2)

## Creating a stem&leaf plot and a corresponding latex table

sink("stem.out")
stem(kalama$height,scale=2)
sink()
latex(readLines("stem.out"))


## Fitting the model

res<-lm(height~age, data=kalama)
kalama.anova<-anova(res)
kalama.summary<-summary(res)

kalama.anova.table<- xtable(kalama.anova)
print(kalama.anova.table)

kalama.summary.table<- xtable(kalama.summary)
print(kalama.summary.table)


plot(growth$age,growth$height, xlab="Age", ylab="Height", main="Height versus Age")

abline(a=65.8,b=0.61,col ="red",lwd=2)
abline(a=63,b=0.71,col ="red",lwd=2)
abline(res,col ="blue",lwd=2)

XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

######### Patien satisfaction

## A hospital administrator wished to study the relation between patient satisfaction ($Y$) 
## and patient's age ($X_1$, in years), severity of illness ($X_2$, an index), and anxiety 
## level ($X_3$, an index).\\[2mm] The administrator randomly selected 46 patients and 
## collected data on the previous variables. Larger values of $Y$, $X_2$, and $X_3$ are, 
## respectively, associated with more satisfaction, increased severity of illness, and more 
## anxiety.

## Reading the data

satisfaction=read.table("satisfaction.txt", header=T)
head(satisfaction,10)

## Exploring the data

cor(satisfaction) 
plot(satisfaction)

## Descriptive Statistics

options(digits=2)
descrip.satisfaction<-stat.desc(satisfaction,basic=TRUE, desc=TRUE)

## Fitting the model

satisfaction.lm<-lm(satis~age+severity+anxiety, data=satisfaction)
satisfaction.summary<-summary(satisfaction.lm)

## Likelihood ratio test null model versus full model

satisfaction.lm.int<-lm(satis~1, data=satisfaction) # Null model
anova(satisfaction.lm.int,satisfaction.lm)          # Null versus full

## Sequential building og the model

satisfaction.anova<-anova(satisfaction.lm)

satisfaction.lm2<-lm(satis~age+anxiety+severity, data=satisfaction)
satisfaction.anova2<-anova(satisfaction.lm2)

## Final model

satisfaction.lm.final<-lm(satis~age+anxiety, data=satisfaction)
satisfaction.final.summary<-summary(satisfaction.lm.final)

## Predicting a new observation

newdata = data.frame(age=43, anxiety=2.7) 
pred.w.plim <- predict(satisfaction.lm.final, newdata, interval="predict") 
pred.w.clim <- predict(satisfaction.lm.final, newdata, interval = "confidence")

help(anova)


######### Smoking and Cancer

### The data are per capita numbers of cigarettes smoked (sold) by 43 states and the 
### District of Columbia in 1960 together with death rates per thouusand population from
### various forms of cancer.

### Number of cases: 44
### Variable Names:

###     CIG = Number of cigarettes smoked (hds per capita)
###     BLAD = Deaths per 100K population from bladder cancer
###     LUNG = Deathes per 100K population from lung cancer
###     KID = Deaths per 100K population from bladder cancer
###     LEUK = Deaths per 100 K population from leukemia

### Reference: J.F. Fraumeni, "Cigarette Smoking and Cancers of the Urinary Tract: Geographic Variations in the United States," Journal of the National Cancer Institute, 41, 1205-1211. 

## Work

FPath.work.smoking.cancer <- "C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\SmokingCancer.txt" 

## Home

FPath.work.smoking.cancer <- "C:\\Users\\ariel alonso\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\SmokingCancer.txt" 

Smoking.cancer = read.table(FPath.work.smoking.cancer, header=T)

## Showing the data

Smoking.cancer.table<- xtable(Smoking.cancer)
print(Smoking.cancer.table)

## Boxplot for cigarrestes

boxplot(Smoking.cancer$CIG  ,main="Number of cigarettes",ylab="Cigarettes", col = "gray")
onevec <- rep(1,length(Smoking.cancer$CIG))
identify(onevec,Smoking.cancer$CIG, plot = T)
points(onevec[7], Smoking.cancer$CIG[7], col = "red", pch=16, cex=1.5)
points(onevec[25], Smoking.cancer$CIG[25], col = "red", pch=16, cex=1.5)

## Boxplot for mortality versus cancer type

boxplot(Smoking.cancer[,3:6]   ,main="Deaths per 100 K population",ylab="Deaths", col = "gray")

## Scatterplot smoking versus lung cancer mortality

plot(Smoking.cancer$CIG,Smoking.cancer$LUNG , xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Lung Cancer")
abline(lm(Smoking.cancer$LUNG~Smoking.cancer$CIG), lwd=2, col="red")
lines(lowess(Smoking.cancer$CIG,Smoking.cancer$LUNG,f=0.9), lwd=2, col="blue")
abline(lm(Smoking.cancer$LUNG~Smoking.cancer$CIG+Smoking.cancer$CIG^3), lwd=2, col="blue")

# pnt <- identify(Smoking.cancer$CIG,Smoking.cancer$LUNG, plot = T)
# it is al identified and pnt is equal to 5

points(Smoking.cancer$CIG[7], Smoking.cancer$LUNG[7], col = "red", pch=16, cex=1.5)
points(Smoking.cancer$CIG[25], Smoking.cancer$LUNG[25], col = "red", pch=16, cex=1.5)


## General descriptive statistics

options(digits=2)
descript.Smoking.cancer<-stat.desc(Smoking.cancer[,c("CIG","BLAD","LUNG","KID", "LEUK")],basic=TRUE, desc=TRUE)
descript.table<- xtable(descript.Smoking.cancer)

## Histogram with Normal Curve

par(mfrow = c(2, 2))
histnorm(Smoking.cancer$LUNG, xlable="Lung cancer", bins=10)
histnorm(Smoking.cancer$BLAD, xlable="Bladder cancer", bins=10)
histnorm(Smoking.cancer$KID, xlable="Kidney cancer", bins=10)
histnorm(Smoking.cancer$LEUK, xlable="Leukemia", bins=10)

## Incidence of smoking versus cancer mortality

par(mfrow = c(2, 2))
plot(Smoking.cancer$CIG,Smoking.cancer$LUNG , xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Lung Cancer")
abline(lm(Smoking.cancer$LUNG~Smoking.cancer$CIG), lwd=2, col="red")
plot(Smoking.cancer$CIG,Smoking.cancer$BLAD, xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Bladder cancer")
abline(lm(Smoking.cancer$BLAD~Smoking.cancer$CIG), lwd=2, col="red")
plot(Smoking.cancer$CIG,Smoking.cancer$KID, xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Kidney cancer")
abline(lm(Smoking.cancer$KID~Smoking.cancer$CIG), lwd=2, col="red")
plot(Smoking.cancer$CIG,Smoking.cancer$LEUK, xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Leukemia")
abline(lm(Smoking.cancer$LEUK~Smoking.cancer$CIG), lwd=2, col="red")

## Calculating the z-scores and correlation by hand between smoking and lung cancer

Lung.Data.Summary<-data.frame(smoking=Smoking.cancer$CIG,devsmok=Smoking.cancer$CIG-mean(Smoking.cancer$CIG),
              smoking.std=scale(Smoking.cancer$CIG), lung=Smoking.cancer$LUNG, devlung=Smoking.cancer$LUNG-mean(Smoking.cancer$LUNG),
              lung.std=scale(Smoking.cancer$LUNG))

Lung.Data.Summary.table<- xtable(Lung.Data.Summary)
print(Lung.Data.Summary.table)

mean(Smoking.cancer$CIG)
sd(Smoking.cancer$CIG)

mean(Smoking.cancer$LUNG)
sd(Smoking.cancer$LUNG)

## Correlation between smoking and cancer

corr.smoking.lung<-cor(Smoking.cancer$CIG,Smoking.cancer$LUNG)
corr.smoking.BLAD<-cor(Smoking.cancer$CIG,Smoking.cancer$BLAD)
corr.smoking.KID<-cor(Smoking.cancer$CIG,Smoking.cancer$KID)
corr.smoking.LEUK<-cor(Smoking.cancer$CIG,Smoking.cancer$LEUK)

## Testing if the population correlation is zero

corr.smoking.lung.test= cor.test(Smoking.cancer$CIG,Smoking.cancer$LUNG, alternative="two.sided", 
                      method = "pearson")

cor.test.LUNG=cor.test(~ CIG+LUNG, data=Smoking.cancer, alternative="two.sided", method = "pearson")
cor.test.BLAD=cor.test(~ CIG+BLAD, data=Smoking.cancer, alternative="two.sided", method = "pearson")
cor.test.KID=cor.test(~ CIG+KID, data=Smoking.cancer, alternative="two.sided", method = "pearson")
cor.test.LEUK=cor.test(~ CIG+LEUK, data=Smoking.cancer, alternative="two.sided", method = "pearson")


## Plotting confidence interval for all groups

ci.Smoking.cancer.new=matrix(0,4,3,byrow = TRUE)
ci.Smoking.cancer.new[1,]=c(corr.smoking.lung,cor.test.LUNG$conf.int)
ci.Smoking.cancer.new[2,]=c(corr.smoking.BLAD,cor.test.BLAD$conf.int)
ci.Smoking.cancer.new[3,]=c(corr.smoking.KID,cor.test.KID$conf.int)
ci.Smoking.cancer.new[4,]=c(corr.smoking.LEUK,cor.test.LEUK$conf.int)
rownames(ci.Smoking.cancer.new)=c("LUNG","BLAD","KID","LEUK")

## Plot the confidence intervals. The first row of ci.Smoking.cancer.new, namely 
## ci.Smoking.cancer.new[2,], has the left end points. 
## The second row of ci.new, namely ci.Smoking.cancer.new[3,], has the right 
## end points.

jpeg('smoking-cancer-corre-ci.jpg')
plotCI(x = 1:4, y = ci.Smoking.cancer.new[,1], xaxt="n", xlab="" , ylab="Correlation with smoking incidence", li = ci.Smoking.cancer.new[,2], ui = ci.Smoking.cancer.new[,3], col = 2, lwd = 3)
axis(1, at=1:4,labels=rownames(ci.Smoking.cancer.new), las=1) 
abline(h=0)
dev.off() 

## Incidence of smoking versus cancer mortality with correlation

par(mfrow = c(2, 2))
plot(Smoking.cancer$CIG,Smoking.cancer$LUNG , xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Lung Cancer")
mtext(paste("r(x,y)=",signif(corr.smoking.lung, digits = 2)))
abline(lm(Smoking.cancer$LUNG~Smoking.cancer$CIG), lwd=2, col="red")

plot(Smoking.cancer$CIG,Smoking.cancer$BLAD, xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Bladder cancer")
mtext(paste("r(x,y)=",signif(corr.smoking.BLAD, digits = 2)))
abline(lm(Smoking.cancer$BLAD~Smoking.cancer$CIG), lwd=2, col="red")

plot(Smoking.cancer$CIG,Smoking.cancer$KID, xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Kidney cancer")
mtext(paste("r(x,y)=",signif(corr.smoking.KID, digits = 2)))
abline(lm(Smoking.cancer$KID~Smoking.cancer$CIG), lwd=2, col="red")

plot(Smoking.cancer$CIG,Smoking.cancer$LEUK, xlab="Cigarettes per capita", ylab="Deaths per 100 K population", pch=16, cex=1.5, main="Smoking and Leukemia")
mtext(paste("r(x,y)=",signif(corr.smoking.LEUK, digits = 2)))
abline(lm(Smoking.cancer$LEUK~Smoking.cancer$CIG), lwd=2, col="red")


######### Muscle mass study

### A nutritionist would like to study the relationship between age and muscle mass. 
### To that effect he randomly selected 60 women from different age groups and recorded the age 
### and muscle mass of each women.  (Kutner et al. pg. 36)
### Number of cases: 60
### Variable Names:

###     muscle: Muscle mass
###     age: Age in years

## Work

FPath.work.muscle.mass  <- "C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\muscle.txt" 

## Home
 
FPath.work.muscle.mass <- "C:\\Users\\ariel alonso\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\muscle.txt" 

muscle<-read.table(FPath.work.muscle.mass, col.names=c("muscle", "age"))

## Boxplot for age and muscle mass

par(mfrow = c(1, 2))
boxplot(muscle$age, main="Age",ylab="", col = "gray")
boxplot(muscle$muscle ,main="Muscle mass",ylab="", col = "gray")

## General descriptive statistics

options(digits=2)
descript.muscle.mass<-stat.desc(muscle[,c("age","muscle")],basic=TRUE, desc=TRUE)
descript.muscle.mass.table<- xtable(descript.muscle.mass)

## Histogram with Normal Curve

setEPS()
postscript("Muscle-histogram.eps",horizontal=FALSE)
par(mfrow = c(1, 2))
histnorm(muscle$age, xlable="Age", bins=10)
histnorm(muscle$muscle, xlable="Muscle mass", bins=10)
dev.off()

pdf("Muscle-histogram.pdf")
par(mfrow = c(1, 2))
histnorm(muscle$age, xlable="Age", bins=10)
histnorm(muscle$muscle, xlable="Muscle mass", bins=10)
dev.off()

png("Muscle-histogram.png")
par(mfrow = c(1, 2))
histnorm(muscle$age, xlable="Age", bins=10)
histnorm(muscle$muscle, xlable="Muscle mass", bins=10)
dev.off()

## Scaterplot

plot(muscle$age,muscle$muscle, ylab="muscle mass", xlab="age", main="Muscle mass versus age")

## Correlation between muscle and age

corr.muscle<-cor(muscle$age,muscle$muscle)

mean(muscle$age)
sd(muscle$age)

mean(muscle$muscle)
sd(muscle$muscle)

## Linear regression model

res.muscle<-lm(muscle$muscle~muscle$age)
summary(res.muscle)

res.muscle.anova<-anova(res.muscle)

res.muscle.anova.table<- xtable(res.muscle.anova)
print(res.muscle.anova.table)


## Incidence of smoking versus cancer mortality with correlation

plot(muscle$age,muscle$muscle, xlab="Age", ylab="Muscle mass", pch=16, cex=1.5, main="Muscle mass versus")
mtext(paste("r(x,y)=",signif(corr.muscle, digits = 2)))
abline(lm(muscle$muscle~muscle$age), lwd=2, col="red")



######### Smoking-Lung-Cancer 2

### Below is a table showing the incidence rates of smoking and lung cancer per 100,000 people, 
### as tracked by the Centers for Disease Control and Prevention, from 1999 - 2007.
### Source smoking: National program of cancer registries, CDC
### http://apps.nccd.cdc.gov/uscs/cantcerbystateandregion.aspx
### Source Lung cancer: National program of cancer registries, CDC
### http://www.cdc.gov/nchs/data/nhis/earlyrelease/earlyrelease201006.pdf#page=52

## The data is generated one time and save it does not have to be generated again

## town<-c(1,2,3,4,5,6,7,8,9)
## smoking<-c(23.3,23.1,22.6,22.3,21.5,20.8,20.8,20.8,19.7)
## lung<-c(93.5,91.5,91.0,89.7,89.3,87.8,86.6,84.2,80.5)

## lung.data<-data.frame(town,smoking,lung)
## write.table(lung.data, file = "LungCancer.txt", sep = "\t")

## Work

FPath.work.lung  <- "C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\LungCancer.txt" 

## Home
 
FPath.work.lung <- "C:\\Users\\ariel alonso\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\LungCancer.txt" 

lung.data<-read.table(FPath.work.lung)


## Scaterplot

plot(lung.data$smokin,lung.data$lung, ylab="Lung Cancer", xlab="Smoking", main="Smoking versus Lung Cancer")


## Calculating the z-scores and correlation by hand between smoking and lung cancer

Lung.Data.Summary2<-data.frame(smoking=lung.data$smoking,devsmok=lung.data$smoking-mean(lung.data$smoking),
              smoking.std=scale(lung.data$smoking), lung=lung.data$lung, devlung=lung.data$lung-mean(lung.data$lung),
              lung.std=scale(lung.data$lung))

Lung.Data.Summary2.table<- xtable(Lung.Data.Summary2)
print(Lung.Data.Summary2.table)

mean(lung.data$smoking)
sd(lung.data$smoking)

mean(lung.data$lung)
sd(lung.data$lung)

## Correlation between smoking and cancer

corr.smoking.lung.cancer<-cor(lung.data$smoking,lung.data$lung)

## Testing if the population correlation is zero

corr.smoking.lung.cancer.test= cor.test(lung.data$smoking,lung.data$lung, alternative="two.sided", 
                      method = "pearson")

## Linear regression model

res.lung<-lm(lung.data$lung~lung.data$smoking)
summary(res.lung)

res.lung.anova<-anova(res.lung)

res.lung.anova.table<- xtable(res.lung.anova)
print(res.lung.anova.table)


## Incidence of smoking versus cancer mortality with correlation

plot(lung.data$smoking,lung.data$lung, xlab="Smoking", ylab="Lung Cancer", pch=16, cex=1.5, main="Smoking and Lung Cancer")
mtext(paste("r(x,y)=",signif(corr.smoking.lung.cancer, digits = 2)))
abline(lm(lung.data$lung~lung.data$smoking), lwd=2, col="red")


######### Larsen and Marx (1986) write

### Since Word War II, plutonium for use in atomic weapons has been produced at an Atomic Energy Commission facility in Hanford, Washington. 
### One of the major safety problems encountered there has been the storage of radioactive wastes. Over the years, significant quantities of these 
### substances - including strontium 90 and cesium 137 - have leaked from their open-pit storage areas into the nearby Columbia River, which flows 
### along the Washington-Oregon border, and eventually empties into the Pacific Ocean.

### To measure the health consequences of this contamination, an index of exposure was calculated for each of the nine Oregon counties having frontage 
### on either the Columbia River or the Pacific Ocean. This particular index was based on several factors, including the county's stream distance from 
### Hanford and the average distance of its population from any water frontage. As a covariate, the cancer mortality rate was determined for each of these same counties.

### The data give the index of exposure and the cancer mortality rate during 1959-1964 for the nine Oregon counties affected. Higher index values represent higher levels of contamination.

### Variable 		Description
### County 		      Name of county
### Exposure 		Index of exposure
### Mortality 		Cancer mortality per 100,000 man-years
	
### Fadeley, R. C. (1965). Oregon malignancy pattern physiographically related to Hanford, Washington, Radioisotope Storage. Journal of Environmental Health 27, 883-897.
### Larsen, R.J., and Marx, M.L., (1986). An Introduction to Mathematical Statistics and Its Applications 2nd Edition. Prentice-Hall, Englewood Cliffs, New Jersey. Case Study 1.2.4.

## Work

FPath.work.hanford<- "C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\CancerHanfordReactor.txt" 


## Home

FPath.work.hanford <- "C:\\Users\\ariel alonso\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\CancerHanfordReactor.txt" 

Cancer.Hanford = read.table(FPath.work.hanford, header=T, sep=",")

plot(Cancer.Hanford$Exposure,Cancer.Hanford$Mortality, xlab="Exposure", ylab="Mortality", pch=16, cex=1.5, main="Kankersterfte in de buurt van Hanford Reactor")

## Create Line Chart

ord1<-order(Cancer.Hanford$Mortality)

Cancer.Hanford1 <- Cancer.Hanford[ord1,]

Cancer.Hanford1$Countynum<-1:9

print(xtable(Cancer.Hanford1))

jpeg('Hanford-Reactor-plot2.jpg')
par(mfrow = c(2, 1))

plot(Cancer.Hanford1$Countynum, Cancer.Hanford1$Exposure , xaxt="n", xlab="", ylab="Exposure" )
lines(Cancer.Hanford1$Countynum, Cancer.Hanford1$Exposure ,  type="b", lwd=2.0, lty=1, col="red")
axis(1, at=Cancer.Hanford1$Countynum,labels=Cancer.Hanford1$County, col.axis="red", las=2)

plot(Cancer.Hanford1$Countynum, Cancer.Hanford1$Mortality , xaxt="n", xlab="", ylab="Mortality")
lines(Cancer.Hanford1$Countynum, Cancer.Hanford1$Mortality , type="b", lwd=2.0, lty=1, col="blue")
axis(1, at=Cancer.Hanford1$Countynum,labels=Cancer.Hanford1$County, col.axis="red", las=2)
dev.off()

## Create Barchart

print(xtable(Cancer.Hanford))

ord<-order(Cancer.Hanford$Mortality)

Cancer.Hanford.table <- xtable(Cancer.Hanford[ord,])
print(Cancer.Hanford.table)

barX.Cancer.Hanford<-barplot(Cancer.Hanford$Mortality[ord],names.arg=1:9,ylim=c(0, max(Cancer.Hanford$Mortality)+10),col=rainbow(12)[ord])
text(x=barX.Cancer.Hanford, y=Cancer.Hanford$Mortality[ord]+3.8, label=Cancer.Hanford$Mortality[ord])
legend("topleft",legend =Cancer.Hanford$County[ord],fill=rainbow(12)[ord],cex=.75)

barX.Cancer.Hanford<-barplot(Cancer.Hanford$Mortality,names.arg=1:9,ylim=c(0, max(Cancer.Hanford$Mortality)+10),col=rainbow(12))
text(x=barX.Cancer.Hanford, y=Cancer.Hanford$Mortality+3.8, label=Cancer.Hanford$Mortality)
legend("topleft",legend =Cancer.Hanford$County,fill=rainbow(12),cex=.75)

## Calculating the z-scores and correlation by hand between smoking and lung cancer

Hanford.Data.Summary<-data.frame(Exposure=Cancer.Hanford$Exposure,devExposure=Cancer.Hanford$Exposure-mean(Cancer.Hanford$Exposure),
              Exposure.std=scale(Cancer.Hanford$Exposure), Mortality=Cancer.Hanford$Mortality, devMortality=Cancer.Hanford$Mortality-mean(Cancer.Hanford$Mortality),
              Mortality.std=scale(Cancer.Hanford$Mortality))

Hanford.Data.Summary.table<- xtable(Hanford.Data.Summary)
print(Hanford.Data.Summary.table)

mean(Cancer.Hanford$Exposure)
sd(Cancer.Hanford$Exposure)

mean(Cancer.Hanford$Mortality)
sd(Cancer.Hanford$Mortality)

## Correlation between expoture and cancer

corr.Hanford<-cor(Cancer.Hanford$Exposure,Cancer.Hanford$Mortality)

## Testing if the population correlation is zero

corr.Hanford.test=cor.test(~ Exposure+Mortality, data=Cancer.Hanford, alternative="two.sided", method = "pearson")


# Regression lines and plots

hanford.res<-lm(Mortality~Exposure, data=Cancer.Hanford)
hanford.anova<-anova(hanford.res)
hanford.summary<-summary(hanford.res)

hanford.anova.table<- xtable(hanford.anova)
print(hanford.anova.table)

hanford.summary.table<- xtable(hanford.summary)
print(hanford.summary.table)

## Expoture versus cancer mortality with correlation

plot(Cancer.Hanford$Exposure,Cancer.Hanford$Mortality, xlab="Exposure", ylab="Mortality", pch=16, cex=1.5, main="Kankersterfte in de buurt van Hanford Reactor")
mtext(paste("r(x,y)=",signif(corr.Hanford, digits = 2)))
abline(lm(Cancer.Hanford$Mortality~Cancer.Hanford$Exposure), lwd=2, col="red")


######### Cancer Survival

### A one-way ANOVA with Organ as the discrete factor and Survival as the dependent variable is appropriate for these data. 
### A normal probability plot of the residuals shows that they are skewed to the right. We transform the data to alleviate this problem. 
### A square root transformation appears to make the residuals more symmetric, so we use the square root of Survival (Survival) as 
### our new dependent variable.

### The ANOVA table for this analysis indicates that the means of Survival are significantly different at the 5% level of significance with 
### a p-value of 0.0002. Figure 1 is a boxplot of ÃSurvival by Organ. We can see from this plot that breast cancer patients seem to survive 
### much longer after receiving treatments of ascorbate than patients with cancer of the bronchus, colon, or stomach. 
### A post-hoc test using Scheffe's method confirms the results of the boxplot. 

### Reference: Cameron, E. and Pauling, L. (1978) Supplemental ascorbate in the supportive treatment of cancer: re-evaluation of prolongation 
### of survival times in terminal human cancer. Proceedings of the National Academy of Science USA, 75, 4538Ð4542. 
### Also found in: Manly, B.F.J. (1986) Multivariate Statistical Methods: A Primer, New York: Chapman & Hall, 11. Also found in: Hand, D.J., et al. (1994) 
### A Handbook of Small Data Sets, London: Chapman & Hall, 255. 

## Work

FPath.survival  <- "C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\CancerSurvival.txt" 

## Home

FPath.survival <- "C:\\Users\\ariel alonso\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\CancerSurvival.txt" 

Cancer.Survival = read.table(FPath.survival, header=T)

head(Cancer.Survival, n=15)#

Cancer.Survival.table <- xtable(head(Cancer.Survival, n=15))
print(Cancer.Survival.table)


## General descriptive statistics

Cancer.Survival1<-Cancer.Survival
Cancer.Survival1$Sur<-Cancer.Survival$Survival
options(digits=1)
descript.cancer.survival<-stat.desc(Cancer.Survival1[,c("Survival","Sur")],basic=TRUE, desc=TRUE)
descript.cancer.survival.table<- xtable(descript.cancer.survival)

Cancer.Survival.organ<-cbind(Freq=table(Cancer.Survival$Organ),
                             relative=100*prop.table(table(Cancer.Survival$Organ)))


Cancer.Survival.organ.table <- xtable(Cancer.Survival.organ)
print(Cancer.Survival.organ.table)

## General descriptive statistics per Organ

co.var <- function(x,na.rm=TRUE) (sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm))
descript.organ<-aggregate(Survival~Organ, data =Cancer.Survival, 
                  FUN = function(x) { c(mean = mean(x), median=median(x),std = sd(x), cv(x))})
descript.organ.table <- xtable(t(descript.organ))

## Histogram with Normal Curve

histnorm(Cancer.Survival$Survival, xlable="Survival", bins=10)

Breast<-Cancer.Survival$Organ=="Breast"
Stomach<-Cancer.Survival$Organ=="Stomach"
Bronchus<-Cancer.Survival$Organ=="Bronchus"
Colon<-Cancer.Survival$Organ=="Colon"
Ovary<-Cancer.Survival$Organ=="Ovary"

par(mfrow = c(2, 3))
histnorm(Cancer.Survival$Survival[Breast], xlable="Breast cancer", bins=10)
histnorm(Cancer.Survival$Survival[Stomach], xlable="Stomach cancer", bins=10)
histnorm(Cancer.Survival$Survival[Bronchus], xlable="Bronchus cancer", bins=10)
histnorm(Cancer.Survival$Survival[Colon], xlable="Colon cancer", bins=10)
histnorm(Cancer.Survival$Survival[Ovary], xlable="Ovary cancer", bins=10)

## Boxplot

boxplot(Survival~Organ,data=Cancer.Survival, main="Survival versus organ",
   xlab="Organ", ylab="Survival", col = "gray")

## Calculating and plotting confidence interval for all groups

w=list(Breast=Cancer.Survival$Survival[Breast],Stomach=Cancer.Survival$Survival[Stomach],
	 Bronchus=Cancer.Survival$Survival[Bronchus],Colon=Cancer.Survival$Survival[Colon],
	 Ovary=Cancer.Survival$Survival[Ovary])

## Apply the myci function to each column of w. Each column of ci 
## has the endpoints of a conficence interval. The first row has the 
## left end points, the second row has the right end points.

ci <- lapply(w, myci)  
average <- lapply(w, mean)

ci.temp=matrix(unlist(ci),5,2,byrow = TRUE)
ci.new=cbind(unlist(average),ci.temp)
ci.new.table <- xtable(ci.new)

## Plot the confidence intervals. The first row of ci.new, namely ci.new[1,], has the left end 
## points. The second row of ci.new, namely ci.new[2,], has the right end points.

jpeg('survival-CI.jpg')
plotCI(x = 1:5, y = ci.new[,1], xaxt="n", xlab="" , ylab="Survival", li = ci.new[,2], ui = ci.new[,3], col = 2, lwd = 3)
axis(1, at=1:5,labels=rownames(ci.new), las=2) 
dev.off() 

## Selecting Bronchus and Colon cancers

Cancer.Survival.BC <- Cancer.Survival[Cancer.Survival$Organ %in% c("Bronchus", "Colon"), ]
Cancer.Survival.BC$Organ<-factor(Cancer.Survival.BC$Organ,levels = c("Bronchus", "Colon"),
					   labels = c("Bronchus", "Colon"))

## Boxplot Bronchus and Colon cancers

boxplot(Survival~Organ,data=Cancer.Survival.BC, main="Survival versus organ",
   xlab="Organ", ylab="Survival", col = "gray")

## Histogram with Normal Curve for  Bronchus and Colon cancers

par(mfrow = c(1, 2))
histnorm(Cancer.Survival.BC$Survival[Cancer.Survival.BC$Organ=="Bronchus"], xlable="Bronchus cancer", bins=10)
histnorm(Cancer.Survival.BC$Survival[Cancer.Survival.BC$Organ=="Colon"], xlable="Colon cancer", bins=10)

## Difference in means

mean(Cancer.Survival.BC$Survival[Cancer.Survival.BC$Organ=="Bronchus"])-mean(Cancer.Survival.BC$Survival[Cancer.Survival.BC$Organ=="Colon"])

## Testing if the variances are equal

var.equal.survival<-leveneTest(Survival~Organ,data=Cancer.Survival.BC)

var.equal.table.survival<-xtable(var.equal.survival, digits=2)

## t-test for equual variance

ind.test.var.equal.survival<-t.test(Survival~Organ, alternative = "two.sided", mu=0, data=Cancer.Survival.BC, var.equal = T)

## t-test for unequual variance

ind.test.var.unequal.survival<-t.test(Survival~Organ, alternative = "two.sided", mu=0, data=Cancer.Survival.BC, var.equal = F)

## Wilcox test to compare survival times for bronchus and colon cancer

wilcox.test(Survival~Organ, data = Cancer.Survival.BC)


######### Normal human body temperature

## What is normal human body temperature (taken orally)? We've all been taught since grade school
## that it's 98.6 degrees Fahrenheit, and never mind that what's normal for one person may not be "normal" for another! 
## So from a statistical point of view, we should abandon the word "normal" and confine ourselves to talking about mean 
## human body temperature. We hypothesize that mean human body temperature is 98.6 degrees, because that's what we were 
## told in the third grade. The data set "Normal Body Temperature, Gender, and Heart Rate" bears on this hypothesis. It is 
## not built-in, so we will have to enter the data. The data are from a random sample (supposedly) of 130 cases and has been 
## posted at the Journal of Statistical Education's data archive. 

## Data structure
## Body temperature (degrees Fahrenheit)
## Gender (1 = male, 2 = female)
## Heart rate (beats per minute)

## Mackowiak, P. A., Wasserman, S. S., and Levine, M. M. (1992). A Critical Appraisal of 98.6 Degrees F, the Upper Limit of the 
## Normal Body Temperature, and Other Legacies of Carl Reinhold August Wunderlich. Journal of the American Medical Association, 
## 268, 1578-1580. 

## Work

FPath.temperature  <- "C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\BodyTemperature.txt" 

## Home

FPath.temperature <- "C:\\Users\\ariel alonso\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\BodyTemperature.txt" 

Body.temperature = read.table(FPath.temperature, header=T)

## Temperature in Celsius

reference.temp<-((98.6 - 32) * 5)/9

Body.temperature$TempC = ((Body.temperature$Temperature - 32) * 5)/9
Body.temperature$GenderF<-factor(Body.temperature$Gender,levels = c(1,2),
					   labels = c("Male", "Female"))


## Data

head(Body.temperature[,c(2,3,4)], n=15)#

Body.temperature.table <- xtable(head(Body.temperature[,c(2,3,4)], n=15))
print(Body.temperature.table)


## General descriptive statistics

options(digits=4)
descript.body.temperature<-stat.desc(Body.temperature[,c("Heart","TempC")],basic=TRUE, desc=TRUE)
descript.body.temperature.table<- xtable(descript.body.temperature)

## Percentile for the  CI

qt(.975, df=129) 

table(Body.temperature$Gender)

## General descriptive statistics per Gender

co.var <- function(x,na.rm=TRUE) (sd(x,na.rm=na.rm)/mean(x,na.rm=na.rm))
Body.temperature.gender<-aggregate(TempC~Gender, data =Body.temperature, 
                  FUN = function(x) { c(mean = mean(x), median=median(x),std = sd(x), cv(x))})
Body.temperature.gender.table <- xtable(t(Body.temperature.gender))


## Histogram with Normal Curve

histnorm(Body.temperature$TempC , xlable="Temperature", bins=10)

## Boxplot

boxplot(Body.temperature$TempC  ,main="Body temperature",ylab="Temperature", col = "gray")
onevec <- rep(1,length(Body.temperature$TempC))
identify(onevec,Body.temperature$TempC, plot = T)
points(onevec[130], Body.temperature$TempC[130], col = "red", pch=16, cex=1.0)
points(onevec[66], Body.temperature$TempC[66], col = "red", pch=16, cex=1.0)
points(onevec[1], Body.temperature$TempC[1], col = "red", pch=16, cex=1.0)

## Boxplot by gender

boxplot(TempC~GenderF, data=Body.temperature  ,main="Body temperature by  gender",ylab="Temperature", col = "gray")

## Plot temperaure versus heart rate

plot(Body.temperature$TempC,Body.temperature$Heart , xlab="Temperature", ylab="Heart rate", pch=16, cex=1.5, main="Heart rate versus temperature")
abline(lm(Body.temperature$Heart~Body.temperature$TempC), lwd=2, col="red")


## One sample t-test

t.test(Body.temperature$TempC, mu=37, alternative="two.sided")

mean(Body.temperature$TempC)

## Comparing the average temperature for boys and girs

## Testing if the variances are equal

var.equal.temp<-leveneTest(TempC~as.factor(Gender),data=Body.temperature)

var.equal.table.temp<-xtable(var.equal.temp, digits=2)

## t-test for equual variance

ind.test.var.equal.temp<-t.test(TempC ~ Gender, alternative = "two.sided", mu=0, data=Body.temperature, var.equal = T)

## t-test for unequual variance

ind.test.var.unequal.temp<-t.test(TempC ~ Gender, alternative = "two.sided", mu=0, data=Body.temperature, var.equal = F)

## 95% CI by gender

Male<-Body.temperature$GenderF=="Male"
Female<-Body.temperature$GenderF=="Female"
w.temp=list(Male=Body.temperature$TempC[Male],Female=Body.temperature$TempC[Female])


ci.temp <- lapply(w.temp, myci)  
average.temp <- lapply(w.temp, mean)

ci.temp=matrix(unlist(ci.temp),2,2,byrow = TRUE)
ci.new.temp=cbind(unlist(average.temp),ci.temp)
ci.new.table.temp <- xtable(ci.new.temp)

## Plot the confidence intervals. The first row of ci.new, namely ci.new[1,], has the left end 
## points. The second row of ci.new, namely ci.new[2,], has the right end points.


plotCI(x = 1:2, y = ci.new.temp[,1], xaxt="n", xlab="" , ylab="Temperature", 
       xlim=c(0.5,2.5), li = ci.new.temp[,2], ui = ci.new.temp[,3], col = 2, lwd = 3,
	main="95% confidence intervals for males and females")
axis(1, at=1:2,labels=rownames(ci.new.temp), las=1) 

## Calculating the z-scores and correlation by hand

Data.Summary.Body<-data.frame(Temperature=Body.temperature$TempC,devTemperature=Body.temperature$TempC-mean(Body.temperature$TempC),
              Temperature.std=scale(Body.temperature$TempC), Heart=Body.temperature$Heart, devHeart=Body.temperature$Heart-mean(Body.temperature$Heart),
              Heart.std=scale(Body.temperature$Heart))

Data.Summary.Body.table<- xtable(Data.Summary.Body)
print(Data.Summary.Body.table)

cov.temperature.heart<-cov(Body.temperature$TempC,Body.temperature$Heart)
corr.temperature.heart<-cor(Body.temperature$TempC,Body.temperature$Heart)

## Testing if the population correlation is zero

corr.temperature.heart.test= cor.test(~ TempC+Heart, data=Body.temperature, alternative="two.sided", method = "pearson")

# Regression lines and plots

temperature.res<-lm(Heart~TempC, data=Body.temperature)
temperature.anova<-anova(temperature.res)
temperature.summary<-summary(temperature.res)

temperature.anova.table<- xtable(temperature.anova)
print(temperature.anova.table)

temperature.summary.table<- xtable(temperature.summary)
print(temperature.summary.table)


######### High Fiber Diet Plan

## A manufacturer was considering marketing crackers high in a certain kind of edible fiber as a dieting aid. 
## Dieters would consume some crackers before a meal, filling their stomachs so that they would feel less hungry 
## and eat less. A laboratory studied whether people would in fact eat less in this way.

## Overweight female subjects ate crackers with different types of fiber (bran fiber, gum fiber, both, and 
## a control cracker) and were then allowed to eat as much as they wished from a prepared menu. The amount of 
## food they consumed and their weight were monitored, along with any side effects they reported.  

## Effects of dietary fiber: 12 female subjects were fed a controlled diet. Before each meal they ate crackers 
## containing either bran fiber, gum fiber, a combination of both, or no fiber (control). Their caloric intake 
## was monitored. Subjects reported any gastric or other problems. 

## Unfortunately, some subjects developed uncomfortable bloating and gastric upset from some of the fiber crackers. A contingency 
## table of "Cracker" versus "Bloat" shows the relationship between the four different types of cracker and the 
## four levels of severity of bloating as reported by the subjects.

## Number of cases: 12
## Variable Names:

## 1. Cracker: Type of fiber in the cracker
## 2. Diet: One of four diets (type of cracker)
## 3. Subject: An identification for each of the 12 subjects
## 4. Digested: Digested calories. Difference between caloric intake and calories passed through system
## 5. Bloat: Degree of bloating and flatulence reported by the subjects 

## Work

FPath.fiber <- "C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\FiberDiet.txt" 

## Home

FPath.fiber <- "C:\\Users\\ariel alonso\\OneDrive\\Equizo\\Courses\\KULeuven\\Verpleegwetenschappen\\Data-Codes\\FiberDiet.txt" 

FiberDiet = read.table(FPath.fiber, header=T)

FiberDiet.table <- xtable(head(FiberDiet, n=15)#)
print(FiberDiet.table)

################################# Graphs used in the slides
##################################################################
##################################################################

## Scatter plot

sigma.s<-matrix(c(1,0.8,0.8,1),2,2)

set.seed(8)
data.temp.s<-rmvnorm(20, mean =c(1,0), sigma =sigma.s)
data.s<-data.frame(x=data.temp.s[,1], y=data.temp.s[,2])
plot(data.s$x,data.s$y, main="Scatterplot", xlab="x", ylab="y", col=1, pch=16, cex=1.5)
abline(lm(data.s$y~data.s$x), lwd=2, col="blue")

## Explaining least squares

sigma.l<-matrix(c(1,0.6,0.6,1),2,2)

set.seed(1)
data.temp.l<-rmvnorm(10, mean =c(1,0), sigma =sigma.l)
data.l<-data.frame(x=data.temp.l[,1], y=data.temp.l[,2])

## Plot an empty chart with tight axis boundaries, and axis lines on bottom and left

plot(data.l$x,data.l$y, type="n", xaxs="i", yaxs="i", xlim=c(min(data.l$x)-0.5, 2.5), ylim=c(min(data.l$y)-0.5, 2),
     bty="l", xlab="x", ylab="y", main="Scatterplot")
points(data.l$x,data.l$y,col=1, pch=16, cex=1.5)

## Fitting a regression model

fit.l<-lm(data.l$y~data.l$x)

## Drawing different lines

abline(-0.5,0.5, lwd=2, col="black")
abline(-0.5,1.0, lwd=2, col="black")
abline(-1,1.0, lwd=2, col="black")
abline(fit.l, lwd=2, col="blue")

## Selecting a point under the line

# This lets you click on the points you want to change the color of.  Right click and select 
"stop" when # you have clicked all the points you want
# pnt <- identify(data.l$x,data.l$y, plot = F)
# it is al identified and pnt is equal to 4

points(data.l$x[7], data.l$y[7], col = "red", pch=16, cex=1.5)
segments(data.l$x[7],min(data.l$y)-0.5,  data.l$x[7], fitted(fit.l)[7], col= 'red', lwd=2, lty="dashed")
segments(data.l$x[7],data.l$y[7],  min(data.l$x)-0.5, data.l$y[7], col= 'red', lwd=2, lty="dashed")
segments(data.l$x[7], fitted(fit.l)[7],  min(data.l$x)-0.5, fitted(fit.l)[7], col= 'red', lwd=2, lty="dashed")
axis(1, at=data.l$x[7], labels=expression(x[1]), cex=2.5)
axis(2, fitted(fit.l)[7]-0.02, labels=expression(hat(y)[1]), cex=3.5, las = 2)
axis(2, data.l$y[7], labels=expression(y[1]), cex=3.5, las = 2)

longa<-abs(fitted(fit.l)[7]-data.l$y[7])
xcenta<-data.l$x[7]+0.1
ycenta<-(fitted(fit.l)[7]+data.l$y[7])/2
curly(N=100,Tilt=0.65,Long=longa+0.52,scale=0.1,xcent=xcenta, ycent=ycenta,theta=2*pi,col="red",lwd=2,grid=FALSE)

text(xcenta+0.2,fitted(fit.l)[7]-(fitted(fit.l)[7]-data.l$y[7])/2, expression(e[1]),col="red", cex=2)


## Selecting a point above the line

# This lets you click on the points you want to change the color of.  Right click and select 
"stop" when # you have clicked all the points you want
# pnt <- identify(data.l$x,data.l$y, plot = F)
# it is al identified and pnt is equal to 5

points(data.l$x[2], data.l$y[2], col = "red", pch=16, cex=1.5)
segments(data.l$x[2],min(data.l$y)-0.5,  data.l$x[2], data.l$y[2], col= 'red', lwd=2, lty="dashed")
segments(data.l$x[2],data.l$y[2],  min(data.l$x)-0.5, data.l$y[2], col= 'red', lwd=2, lty="dashed")
segments(data.l$x[2], fitted(fit.l)[2],  min(data.l$x)-0.5, fitted(fit.l)[2], col= 'red', lwd=2, lty="dashed")
axis(1, at=data.l$x[2], labels=expression(x[2]), cex=2.5)
axis(2, fitted(fit.l)[2], labels=expression(hat(y)[2]), cex=3.5, las = 2)
axis(2, data.l$y[2], labels=expression(y[2]), cex=3.5, las = 2)

longa<-abs(fitted(fit.l)[2]-data.l$y[2])
xcenta<-data.l$x[2]-0.1
ycenta<-(fitted(fit.l)[2]+data.l$y[2])/2
curly(N=100,Tilt=0.4,Long=4.1,scale=0.05,xcent=xcenta, ycent=ycenta,theta=pi,col="red",lwd=2,grid=FALSE)

text(xcenta-0.2,data.l$y[2]-(-fitted(fit.l)[2]+data.l$y[2])/2, expression(e[2]),col="red", cex=2)

## Explaining different sources of variation

sigma.g<-matrix(c(1,0.8,0.8,1),2,2)

set.seed(3)
data.temp.g<-rmvnorm(10, mean =c(25,40), sigma =sigma.g)
data.g<-data.frame(x=data.temp.g[,1], y=data.temp.g[,2])

plot(data.g$x,data.g$y, type="n", xaxs="i", yaxs="i", xlim=c(23, 27), ylim=c(38, 42),
     bty="l", xlab="", ylab="Weight", main="Scatterplot")
points(data.g$x,data.g$y,col=1, pch=16, cex=1.5)

meang.x<-mean(data.g$x)
meang.y<-mean(data.g$y)
abline(h=meang.y, lwd=2)

axis(2, at=meang.y, labels=expression(bar(y)), cex=3.5, las = 2)

## Total variation

abline(h=meang.y, lwd=2)
axis(2, at=meang.y, labels=expression(bar(y)), cex=3.5, las = 2)
for(i in 1:10)
{
arrows(data.g$x[i],data.g$y[i], data.g$x[i], meang.y, lwd=2, length=0.1,angle=20, col="black")
}

## Variation due to regression

plot(data.g$x,data.g$y, type="n", xaxs="i", yaxs="i", xlim=c(23, 27), ylim=c(38, 42),
     bty="l", xlab="Size", ylab="Weight", main="Scatterplot")
points(data.g$x,data.g$y,col=1, pch=16, cex=1.5)

meang.x<-mean(data.g$x)
meang.y<-mean(data.g$y)
abline(v=meang.x, lwd=2)
abline(h=meang.y, lwd=2)

axis(1, at=meang.x, labels=expression(bar(x)), cex=2.5)
axis(2, at=meang.y, labels=expression(bar(y)), cex=3.5, las = 2)

fit.g<-lm(data.g$y~data.g$x)
abline(fit.g, lwd=2, col="red")

text(26.2,41.67, expression(paste(hat(y),"=", b[0]+b[1]*x)),col="red", cex=2)

for(i in 1:10)
{
arrows(data.g$x[i],fitted(fit.g)[i], data.g$x[i], meang.y, lwd=2, length=0.1,angle=20, col="red")
}


## Variation from regression line

plot(data.g$x,data.g$y, type="n", xaxs="i", yaxs="i", xlim=c(23, 27), ylim=c(38, 42),
     bty="l", xlab="Size", ylab="Weight", main="Scatterplot")
points(data.g$x,data.g$y,col=1, pch=16, cex=1.5)

meang.x<-mean(data.g$x)
meang.y<-mean(data.g$y)
abline(v=meang.x, lwd=2)
abline(h=meang.y, lwd=2)

axis(1, at=meang.x, labels=expression(bar(x)), cex=2.5)
axis(2, at=meang.y, labels=expression(bar(y)), cex=3.5, las = 2)

fit.g<-lm(data.g$y~data.g$x)
abline(fit.g, lwd=2, col="red")

text(26.2,41.67, expression(paste(hat(y),"=", b[0]+b[1]*x)),col="red", cex=2)

for(i in 1:10)
{
arrows(data.g$x[i],data.g$y[i], data.g$x[i], fitted(fit.g)[i], lwd=2, length=0.1,angle=20, col="blue")
}

## Example of total variation

plot(data.g$x,data.g$y, type="n", xaxs="i", yaxs="i", xlim=c(23, 27), ylim=c(38, 42),
     bty="l", xlab="Size", ylab="Weight", main="Scatterplot")
points(data.g$x,data.g$y,col=1, pch=16, cex=1.5)

meang.x<-mean(data.g$x)
meang.y<-mean(data.g$y)
abline(v=meang.x, lwd=2)
abline(h=meang.y, lwd=2)

axis(1, at=meang.x, labels=expression(bar(x)), cex=2.5)
axis(2, at=meang.y, labels=expression(bar(y)), cex=3.5, las = 2)

fit.g<-lm(data.g$y~data.g$x)
abline(fit.g, lwd=2, col="red")

text(26.2,41.67, expression(paste(hat(y),"=", b[0]+b[1]*x)),col="red", cex=2)

arrows(data.g$x[6],data.g$y[6], data.g$x[6], fitted(fit.g)[6], lwd=2, length=0.1,angle=20, col="blue")
arrows(data.g$x[6],fitted(fit.g)[6], data.g$x[6], meang.y, lwd=2, length=0.1,angle=20, col="red")

## Model explains the total variation

sigma.t<-matrix(c(1,1,1,1),2,2)

set.seed(4)
data.temp.t<-rmvnorm(10, mean =c(25,40), sigma =sigma.t)
data.t<-data.frame(x=data.temp.t[,1], y=data.temp.t[,2])

plot(data.t$x,data.t$y, type="n", xaxs="i", yaxs="i", xlim=c(23, 27), ylim=c(38, 42),
     bty="l", xlab="Size", ylab="Weight", main="Scatterplot")
points(data.t$x,data.t$y,col=1, pch=16, cex=1.5)

meant.x<-mean(data.t$x)
meant.y<-mean(data.t$y)
abline(v=meant.x, lwd=2)
abline(h=meant.y, lwd=2)

axis(1, at=meant.x, labels=expression(bar(x)), cex=2.5)
axis(2, at=meant.y, labels=expression(bar(y)), cex=3.5, las = 2)

fit.t<-lm(data.t$y~data.t$x)
abline(fit.t, lwd=2, col="red")

text(26.1,41.8, expression(paste(hat(y),"=", b[0]+b[1]*x)),col="red", cex=2)

Plotting parallel lines

plot(1:40,1:40,type="n",, xlim=c(15, 45), ylim=c(10, 40),xlab="Age",ylab="Hospital stay")
abline(a=10,b=0.8, lwd=2,col="red")
abline(a=1,b=0.8, lwd=2,col="blue")
abline(a=1,b=1, lwd=2,col="blue") #interaction model
legend("topleft", inset=.05, ,legend = c("Males","Females"),horiz=TRUE
         , text.col = c("red", "blue")
         , pt.bg = c("red","blue")
         , pch = c(16,16),col = c("red", "blue"))


############################################################################################### 
###############################################################################################
###############################################################################################
###############################################################################################



## Momentarily no more analysis for this data set



# Example 2: Smoking and lung cancer

# Below is a table showing the incidence rates of smoking and lung cancer per 100,000 people, 
# as tracked by the Centers for Disease Control and Prevention, from 1999 - 2007.
# Source smoking: National program of cancer registries, CDC
# http://apps.nccd.cdc.gov/uscs/cantcerbystateandregion.aspx
# Source Lung cancer: National program of cancer registries, CDC
# http://www.cdc.gov/nchs/data/nhis/earlyrelease/earlyrelease201006.pdf#page=52

#year<-c(1999,2000,2001,2002,2003,2004,2005,2006,2007)
town<-c(1,2,3,4,5,6,7,8,9)

smoking<-c(23.3,23.1,22.6,22.3,21.5,20.8,20.8,20.8,19.7)
lung<-c(93.5,91.5,91.0,89.7,89.3,87.8,86.6,84.2,80.5)
lung.data<-data.frame(town,smoking,lung)

lung.data.table<- xtable(lung.data)
print(lung.data.table)

# Create Line Chart

par(mfrow = c(2, 1))

plot(lung.data$town, lung.data$lung, type="n", xlab="Town", ylab="Incidence Cancer" )
lines(lung.data$town, lung.data$lung, type="b", lwd=2.0, lty=1, col="red")

plot(lung.data$town, lung.data$smoking, type="n", xlab="Town", ylab="Incidence Smoking" )
lines(lung.data$town, lung.data$smoking, type="b", lwd=2.0, lty=1, col="blue")

# For some misterious reason when you save a boxplot in R into a eps latex later rotates that figure
# and the entire slide. The solution is to save the boxplot into a pdf like the code below does and then
# to open it with Acroba and save it as an eps.

pdf(file="lungboxplot.pdf")
par(mfrow = c(1, 2))
boxplot(lung.data$lung,main="Incidence of Cancer",ylab="Incidence", col = "gray")
boxplot(lung.data$smoking,main="Incidence of Smoking",ylab="Incidence", col = "gray")
dev.off()

# Histograms for heartbeat and temperature

par(mfrow = c(1, 2))
colors <- rainbow(6) 
hist(lung.data$lung, right=FALSE, freq = TRUE, labels = TRUE,  col=colors, main="Histogram for lung cancer", xlab="Incidence")                  
hist(lung.data$smoking, right=FALSE, freq = TRUE, labels = TRUE, col=colors, main="Histogram for smoking", xlab="Incidence")


# Summary statistics and constructing table with them

Data.lung.Summary<-data.frame(lung=lung.data$lung,devlung=lung.data$lung-mean(lung.data$lung),
              lung.std=scale(lung.data$lung), smoking=lung.data$smoking, devsmoking=lung.data$smoking-mean(lung.data$smoking),smoking.std=scale(lung.data$smoking))

Data.lung.Summary.table<- xtable(Data.lung.Summary)
digits(Data.lung.Summary) <- 2 
print(Data.lung.Summary.table)

mean(lung.data$lung)
sd(lung.data$lung)

mean(lung.data$smoking)
sd(lung.data$smoking)

corr.lung.smoking<-cor(lung.data$lung,lung.data$smoking)

(sd(lung.data$lung)/sd(lung.data$smoking))*corr.lung.smoking

res.lung<-lm(lung.data$lung~lung.data$smoking)
summary(res.lung)

plot(lung.data$smoking,lung.data$lung, ylab="Lung cancer", xlab="Smoking", main="Lung cancer versus smoking")

abline(res.lung,col ="blue",lwd=2.0)



# Example 3

temperature<-c(35.7, 35.9, 36.1, 36.1, 36.2, 36.2, 36.2, 36.2, 36.3, 36.3, 36.3, 36.3, 36.3, 36.4)
heart<-c(70, 71, 74, 80, 73, 75, 82, 64, 69, 70, 68, 72, 78, 70)
heart.data<-data.frame(temperature, heart)

heart.data.table<- xtable(heart.data)
print(heart.data.table)

# Summary statistics and constructing table with them

Data.heart.Summary<-data.frame(heart=heart.data$heart,devheart=heart.data$heart-mean(heart.data$heart),
              heart.std=scale(heart.data$heart), temperature=heart.data$temperature, devtemperature=heart.data$temperature-mean(heart.data$temperature),temperature.std=scale(heart.data$temperature))

Data.heart.Summary.table<- xtable(Data.heart.Summary)
print(Data.heart.Summary.table)

mean(heart.data$heart)
sd(heart.data$heart)

mean(heart.data$temperature)
sd(heart.data$temperature)


res.heart<-lm(heart.data$heart~heart.data$temperature)
summary(res.heart)

plot(heart.data$temperature,heart.data$heart, ylab="Heartbeat", xlab="Temperature", main="Heartbeat versus temperature")

abline(res.heart,col ="blue")

mean(heart.data$heart)
sd(heart.data$heart)
heart.center<-(heart.data$heart-mean(heart.data$heart))
heart.std<-scale(heart.data$heart)

mean(heart.data$temperature)
sd(heart.data$temperature)
temperature.center<-(heart.data$temperature-mean(heart.data$temperature))
temperature.std<-scale(heart.data$temperature)

# For some misterious reason when you save a boxplot in R into a eps latex later rotates that figure
# and the entire slide. The solution is to save the boxplot into a pdf like the code below does and then
# to open it with Acroba and save it as an eps.

pdf(file="heartboxplot.pdf")
par(mfrow = c(1, 2))
boxplot(heart.data$temperature,main="Temperature",ylab="Temperature", col = "gray")
boxplot(heart.data$heart,main="Heartbeat",ylab="Heartbeat", col = "gray")
dev.off()

# Histograms for heartbeat and temperature

par(mfrow = c(1, 2))
colors = c("red", "yellow", "green", "violet", "orange")
hist(heart.data$heart, right=FALSE, freq = TRUE, labels = TRUE, col=colors, main="Histogram for heartbeat", xlab="Heartbeat")                  
colors = c("red", "yellow", "green", "violet", "orange", "pink")
hist(heart.data$temperature, right=FALSE, freq = TRUE, labels = TRUE, col=colors, main="Histogram for temperature", xlab="Temperature")



stem(heart.data$temperature)
stem(heart.data$heart)


cor(heart.data$temperature,heart.data$heart)






