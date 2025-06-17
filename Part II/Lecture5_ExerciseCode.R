library(tree)
library(randomForest)
library(gbm)


heart = read.csv('SAheart.data', header=T,row.names=1)
chdOriginal=heart$chd
heart$chd=as.factor(heart$chd)
set.seed(1)
train=sample(1:nrow(heart), nrow(heart)*2/3)
test=(-train)

tr=tree(chd~.,data=heart, subset=train)
summary(tr)
plot(tr)
text(tr,pretty=0,cex=0.5)

cv.tr=cv.tree(tr,FUN=prune.misclass)
plot(cv.tr)

prune.tree=prune.misclass(tr,best=6)
plot(prune.tree)
text(prune.tree,pretty=0,cex=0.65)
summary(prune.tree)
predtree = predict(tr,newdata=heart[test,],type="class")
table(heart$chd[test],predtree) 
performanceTree1=length(which(predtree==heart$chd[test]))/length(heart$chd[test])

predtree2 = predict(prune.tree,newdata=heart[test,],type="class")
table(heart$chd[test],predtree2) 
performanceTree2=length(which(predtree2==heart$chd[test]))/length(heart$chd[test])

#do a randomForest
set.seed(1)
rf=randomForest(chd~.,mtry=(ncol(heart[train,])/3),data=heart,subset=train,importance=TRUE,ntree=1000,nodesize=20) 
predRF = predict(rf,newdata=heart[test ,])
table(heart$chd[test],predRF) 
performanceRF=length(which(predRF==heart$chd[test]))/length(heart$chd[test])
performanceRF

#do a bagged tree
set.seed(1)
bag=randomForest(chd~.,mtry=ncol(heart[,])-1,data=heart,subset=train,importance=TRUE,ntree=1000,nodesize=20) 
predbag = predict(bag,newdata=heart[test ,])
table(heart$chd[test],predbag) 
performanceBag=length(which(predbag==heart$chd[test]))/length(heart$chd[test])
performanceBag

importance(rf)
varImpPlot (rf)
importance(bag)
varImpPlot (bag)

#BOOSTING
set.seed(1)
heart$chd=chdOriginal
gbm1 <-gbm(chd~.,data=heart[train,],				
				distribution="bernoulli",     
				n.trees=1000,                # number of trees
				shrinkage=0.01,              # shrinkage or learning rate, 0.001 to 0.1 usually work
				interaction.depth=3,         # 1: additive model, 2: two-way interactions, etc.
				n.minobsinnode = 10,         # minimum total weight needed in each node
				cv.folds = 5,                # do 5-fold cross-validation
				keep.data=TRUE,              # keep a copy of the dataset with the object
				verbose=FALSE,               # don't print out progress
				n.cores=1)                   # use only a single core (detecting #cores is error-prone, so avoided here)

# check performance using cross-validation
best.iter <- gbm.perf(gbm1,method="cv")
print(best.iter)
predBO=predict(gbm1,newdata=heart[test,],n.trees=best.iter,type="response")
pred=rep("1",length(heart$chd[test]))
pred[predBO<0.5]="0"
table(heart$chd[test],pred)
performanceBO=length(which(pred==heart$chd[test]))/length(heart$chd[test])
summary(gbm1,n.trees=best.iter)
print(pretty.gbm.tree(gbm1,1))
print(pretty.gbm.tree(gbm1,gbm1$n.trees))
plot(gbm1,c(7,9),best.iter)
plot(gbm1,1,best.iter)
plot(gbm1,2:3,best.iter)


