# TODO: Add comment
# 
# Author: jelier
###############################################################################

library(glmnet)
library(pROC)
library(doBy)
library(pls)
load('/Users/robjelier/Projects/Courses/LM&GLM/2021/Lecture_2-3/VIJVER.Rdata')
dim(data) # 188 4949 ==> 4948 variables for only 188 data points
#nearly 1*10^6 datapoints!
colnames(data)[1:100]
npts = nrow(data) # number of data points
nvars = ncol(data) - 1 # number of variables (each var is a gene expr)

geneExpr=as.matrix(data[,-1])#no first column
meta= data[,1]#response
summary(meta)

#3. For a couple of genes evaluate association with the phenotype. 
#Do you see proof for some predictive potential? 
#Test your intuition with a formal statistical test.

index=9
plot(geneExpr[,index]~meta)
model=glm(meta~geneExpr[,index],family="binomial")
summary(model)

#4. Demonstrate if collinearity occurs between genes in this dataset
correlationMatrix=(cor(geneExpr))
MT1=apply(correlationMatrix>0.9,2,sum)
max(MT1)
length(which(MT1>1))/length(MT1)
#7% of genes have at least one gene for which correlation >0.9! 
MT2=apply(correlationMatrix>0.5,2,sum)
length(which(MT2>1))/length(MT2)
#82%

#illustrate some of the weaker correlations would already be identified as significant

model=lm(geneExpr[,3]~geneExpr[,4])
summary(model)
cor(geneExpr[,3],geneExpr[,4])
plot(geneExpr[,3]~geneExpr[,4])
abline(model)



#5. Use lasso, ridge and PCR methodology and make a predictor based on the gene expression values. 
#How many genes are used for an optimal predictor? 
#Evaluate the performance of the predictors, and comment on what you find. 

#make test and training set
set.seed(1)
train=sample(1:nrow(geneExpr), nrow(geneExpr)*2/3)
test=(-train)

#do a lasso!
grid=10^seq(10,-2,length=100)

lasso.mod=glmnet(y=meta[train],x=(geneExpr[train,]),alpha=1,family="binomial")#note error family
plot(lasso.mod)
cv.lasso=cv.glmnet(geneExpr[train ,],meta[train],alpha=1,family="binomial")
plot(cv.lasso)
lasso.pred=predict(lasso.mod,s=cv.lasso$lambda.min,newx=geneExpr[test,],type="response")
plot(lasso.pred~meta[test])
pred=rep("DM",length(meta[test]))
pred[lasso.pred>0.5]="NODM" #set threshold of 0.5
table(meta[test],pred)
performanceLasso=length(which(pred==meta[test]))/length(meta[test])
performanceLasso #accuracy, 0.71
#Receiver operating statistics!
rocLasso=roc((meta[test]), lasso.pred[,1]) #AUC 0.7357
plot(rocLasso)
#selected genes
vals=predict(lasso.mod,s=cv.lasso$lambda.min,type="coefficients")
selected=colnames(geneExpr)[vals@i]

#do a ridge!
ridge.mod=glmnet(y=meta[train],x=(geneExpr[train,]),alpha=0,family="binomial")
plot(ridge.mod)
ridge.cv=cv.glmnet(geneExpr[train ,],meta[train],alpha=0,family="binomial")
plot(ridge.cv)
ridge.pred=predict(ridge.mod,s=ridge.cv$lambda.min,newx=geneExpr[test,],type="response")
pred=rep("DM",length(meta[test]))
pred[ridge.pred>0.5]="NODM"
table(meta[test],pred)
performanceRidge=length(which(pred==meta[test]))/length(meta[test])
performanceRidge #accuracy, 0.63

#Receiver operating statistics!
rocRidge=roc((meta[test]), ridge.pred[,1]) #AUC 0.7173
plot(rocRidge)

#do a pcr
metaRef= meta=="NODM"
DF=data.frame(meta=metaRef,geneExpr)
pcr.mod <- pcr(meta ~ .,family=binomial(link=logit), data=DF[train,], subset=train, scale=TRUE,validation="CV")
summary(pcr.mod)
validationplot(pcr.mod,val.type="RMSEP")#14 components
pcr.pred <- predict(pcr.mod,data[test,], ncomp=14,type="response")
pred=rep("DM",length(meta[test]))
hist(pcr.pred)#Oops need pls-glm plsRglm for logistic regression!
library(Compositional)
pcr.mod <- glm.pcr(metaRef[train],geneExpr[train,],k=14,xnew=geneExpr[test,])
pcr.pred <- pcr.mod$est
pred=rep("DM",length(meta[test]))
hist(pcr.pred)#Oops need pls-glm plsRglm for logistic regression!
pred[pcr.pred>0.5]="NODM"
table(meta[test],pred)
performancePCR=length(which(pred==meta[test]))/length(meta[test])
performancePCR #accuracy 0.7

#Receiver operating statistics!
rocPCR=roc((meta[test]), pcr.pred) #AUC 0.74
plot(rocPCR)

