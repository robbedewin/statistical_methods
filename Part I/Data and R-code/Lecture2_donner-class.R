## Home

setwd("C:\\Users\\u0065129\\OneDrive\\Equizo\\Courses\\KULeuven\\Bioinformatics\\Logistic-Regression\\Doner-Data-Programs")

## Reading the data

donner<-read.table("donner-class.txt", row.names = 1, header=TRUE) 

head(donner,10)

install.packages("qpcR")
library(qpcR)

install.packages("epicalc")


## Keeping only the variables of interest

donner.na<-na.omit(subset(donner,select=c('Age','Outcome','Sex'))) 
donner.na$fem = as.numeric(donner.na$Sex=="Female")
head(donner.na,10)

## Graph plotting Survival versus gender

cols=ifelse(donner.na$fem==1,"red","black")
plot(donner.na$Age,donner.na$Outcome,col=cols,pch=20,cex=1.2,xlab="Age",ylab="Survival")
legend("topright",pch=20,col=c("red","black"),c("women","men"))

## Graph plotting Survival versus gender but jitteriing the Survival variable

plot(donner.na$Age,jitter(donner$Outcoome,.2),col=cols,pch=20,cex=1.2,xlab="Age",ylab="Status (jittered)")
legend("topright",pch=20,col=c("red","black"),c("women","men"))

## Fitting a linear model and adding the linear graph to the previous plot 

donner.lm<-lm(Outcome ~ Age*fem,data=donner.na)
summary(donner.lm)
c<-coef(donner.lm)
cols=ifelse(donner.na$fem==1,"red","black")
plot(donner.na$Age,donner.na$Outcome,col=cols,pch=20,cex=1.2,xlab="Age",ylab="Survival")
curve(c[1]+c[2]*x+c[3]*0+c[4]*x*0, 0, 70, add=T)
curve(c[1]+c[2]*x+c[3]*1+c[4]*x*1, 0, 70, add=T,col="red")
legend("topright",pch=20,col=c("red","black"),c("women","men"))

## Plotting the logit curve

logit<-function(x){log(x/(1-x))}
ilogit<-function(x,a,b){exp(a+b*x)/(1+exp(a+b*x))}
x=seq(-5,5,.01)
y=ilogit(x,0,1)
plot(x,y,type="l",xlab=expression(X),ylab="Probability")
#plot(x,y,type="l",xlab=expression(eta(X)),ylab="Probability")
#legend("bottomright","pi = logit(eta)")

## Fitting a logistic regression

xtabs(~Outcome+fem+cut(Age,breaks=quantile(Age,probs=seq(0,1,0.25))), data = donner.na)

donner.log<-glm(Outcome ~ Age + fem,data=donner.na,family=binomial(link="logit"))
summary(donner.log)


donner.log.age=glm(Outcome ~ Age,data=donner.na,family=binomial(link="logit"))
anova(donner.log.age,donner.log,test='Chisq')

library(aod)
wald.test(b=coef(donner.log), Sigma=vcov(donner.log), Terms=3)

## Odds ratios

exp(donner.log$coefficients)
exp(confint(donner.log))
exp(cbind(OR =donner.log$coefficients, confint(donner.log)))

## Odd ratio for Survival after 10 years increased

exp(donner.log$coefficients*10)
exp(c(OR =donner.log$coefficients[2]*10, confint(donner.log)[2,]*10))

## Plotting survival for men versus women

cl=coef(donner.log)
plot(donner.na$Age,jitter(donner.na$Outcome,.2),col=cols,pch=20,cex=1.2,xlab="Age",ylab="Status (jittered)")
curve(ilogit(cl[1]+cl[2]*x+cl[3]*0,0,1),add=T)
curve(ilogit(cl[1]+cl[2]*x+cl[3]*1,0,1),add=T,col="red")
legend("topright",pch=20,lty="solid",col=c("red","black"),c("women","men"))


## Predicted probabilities of survival

newdata2<-data.frame(fem=1, Age=mean(donner.na$Age))
newdata2$greP<-predict(donner.log,newdata=newdata2,type="response")

newdata3<-data.frame(fem=0, Age=mean(donner.na$Age))
newdata3$greP<-predict(donner.log,newdata=newdata3,type="response")

newdata4<-data.frame(fem=c(0,1),Age=mean(donner.na$Age))
newdata4$greP<-predict(donner.log,newdata=newdata4,type="response")

## Interaction model

m4<-glm(Outcome ~ Age*fem,data=donner.na,family=binomial(link="logit"))
summary(m4)

## Akaike weights

#install.packages("qpcR") ## old version
#library(qpcR) ## old version

#install.packages("AICcmodavg")
library(AICcmodavg)

## Fitting the models

donner.list=list()

donner.list[[1]]=glm(Outcome ~ Age,data=donner.na,family=binomial(link="logit"))
donner.list[[2]]=glm(Outcome ~ fem,data=donner.na,family=binomial(link="logit"))
donner.list[[3]]=glm(Outcome ~ Age + fem,data=donner.na,family=binomial(link="logit"))
donner.list[[4]]=glm(Outcome ~ Age*fem,data=donner.na,family=binomial(link="logit"))

donner.modnames <- c("Age", "Sex", "Age+Sex", "Age+Sex+Age:Sex")

## Preparing a data set for the AIC (by hand)

aics<-data.frame(paste("m",1:4,sep=""),c(m1$aic,m2$aic,m3$aic,m4$aic),row.names=NULL)
colnames(aics)<-c("model","AIC")
aics<-aics[order(aics$AIC),]
for(i in 1:dim(aics)[1]){
aics$diff[i]<-aics$AIC[i]-aics$AIC[1]}
aics$wi<-exp(-0.5*aics$diff)
aics$aic.weights<-aics$wi/sum(aics$wi)
aics

## Akaike weights with AICcmodavg

donner.aictab=aictab(cand.set = donner.list, modnames = donner.modnames)
donner.aictab

## Model average results

modavg(cand.set= donner.list, parm="Age", second.ord=TRUE, modnames = donner.modnames, 
uncond.se="revised", exclude = list("Age:fem"), conf.level=0.95, warn = TRUE)

modavg(cand.set= donner.list, parm="fem", second.ord=TRUE, modnames = donner.modnames, 
uncond.se="revised", exclude = list("Age:fem"), conf.level=0.95, warn = TRUE)

?modavg

## Akaike weights with qpcR ## old version

# modList <- list(m1,m2,m3,m4) ## old version
# aics2 <- sapply(modList, function(x){AIC(x)}) ## old version
# akaike.weights(aics)## old version

## Gender odd ratio in the interaction model

x=seq(1,70,0.01)
y=exp(coef(m4)[3]+coef(m4)[4]*x)

plot(x,y, type = "n", ylim=c(0.7, 4.5), xlab = "Age", ylab = "Odds ratio", main="Odds ratio for gender")
lines(x,y, lty = 1, col="red")
abline(h=1)


############################################################################################
############################################################################################
############################################################################################
############################################################################################



select.forward<-step (glm(Outcome ~ 1 , data=donner.na),scope=list(lower=~1,upper=~Age*fem), direction='forward')
select.backward<-step (glm(Outcome ~ Age*fem , data=donner.na),scope=list(lower=~1,upper=~Age*fem), direction='backward')
add1(m1, ~.^2,test="Chisq")


# Some graphical output for the final model. Not possible with binary data

donner.log<-glm(Outcome ~ Age + fem,data=donner.na,family=binomial(link="logit"))

donner.log.fitted1 <- predict(donner.log, type = "response")

myplot <- function(role.fitted) 
{ f <- donner.na$fem == 1
  plot(donner.na$Age, donner.log.fitted1, type = "n", ylab = "Probability of surviving", xlab = "Age", ylim = c(0,1))
  lines(donner.na$Age[!f], donner.log.fitted1[!f], lty = 1)
  lines(donner.na$Age[f], donner.log.fitted1[f], lty = 2)
  lgtxt <- c("Fitted (Males)", "Fitted (Females)")
  legend("topright", lgtxt, lty = 1:2, bty = "n")
}
  y <- womensrole$agree / (womensrole$agree +
+ womensrole$disagree)
+ text(womensrole$education, y, ifelse(f, "\\VE", "\\MA"),
+ family = "HersheySerif", cex = 1.25)
+ }

myplot(donner.log.fitted1)

# Plotting the probability curves in exercise 4 in the second practical class

ilogit<-function(x,a,b){exp(a+b*x)/(1+exp(a+b*x))}
x=seq(2,70,.01)

## For women 

y=ilogit(x,1.876,-0.048)
plot(x,y,type="l",xlab="age",ylab="Probability",col="red")
legend("topright",lty="solid",col="red","Vrouwen")


## For men 

y=ilogit(x,0.398,-0.028)
plot(x,y,type="l",xlab="age",ylab="Probability")
legend("topright",lty="solid","Mannen")


## Men versus women 

x=seq(2,70,.01)
ywomen=ilogit(x,1.876,-0.048)
ymen=ilogit(x,0.398,-0.028)

plot(x,ywomen,type="n",xlab="Age",ylab="Probability")
curve(ilogit(x,0.398,-0.028),add=T)
curve(ilogit(x,1.876,-0.048),add=T,col="red")
legend("topright",pch=20,lty="solid",col=c("red","black"),c("Vrouwen","Mannen"))

plot(x,y,type="l",xlab="age",ylab="Probability")
legend("topright","Vrouwen")


cl=coef(donner.log)
plot(donner.na$Age,jitter(donner.na$Outcome,.2),col=cols,pch=20,cex=1.2,xlab="Age",ylab="Status (jittered)")
curve(ilogit(cl[1]+cl[2]*x+cl[3]*0,0,1),add=T)
curve(ilogit(cl[1]+cl[2]*x+cl[3]*1,0,1),add=T,col="red")
legend("topright",pch=20,lty="solid",col=c("red","black"),c("women","men"))

#SPSS class

Flushots<-read.table("Flushots.txt", header=TRUE,sep=",") 

Flushots.int<-glm(Shot ~1 ,data=Flushots,family=binomial(link="logit"))
Flushots.Age<-glm(Shot ~ Age,data=Flushots,family=binomial(link="logit"))
Flushots.Aware<-glm(Shot ~ Awareness,data=Flushots,family=binomial(link="logit"))
Flushots.Age.Aware<-glm(Shot ~ Age+Awareness,data=Flushots,family=binomial(link="logit"))

AIC(Flushots.Age)


