library(SIS)
library(glmnet)
library(e1071)
library(class)
library(tree)
library(randomForest)
library(gbm)
library(nnet)
library(klaR)
library(caret)
### training data
filename <- paste("http://pubs.broadinstitute.org/mpr/projects/",
                  "Global_Cancer_Map/GCM_Training.res", sep="")
dat0 <- read.delim(filename, sep="\t", header=FALSE, skip=3, quote="")
tmp <- dat0[,1:290]
tmp <- tmp[, -seq(4, 290, by=2)]
tmp <- tmp[, -(1:2)]
train <- t(tmp)
filename <- paste("http://pubs.broadinstitute.org/mpr/projects/",
                  "Global_Cancer_Map/GCM_Training.cls", sep="")
train.classes <- read.table(filename, skip=2)+1
train.classes <- unlist(train.classes)
### test data
filename <- paste("http://pubs.broadinstitute.org/mpr/projects/",
                  "Global_Cancer_Map/GCM_Test.res", sep="")
dat1 <- read.delim(filename, sep="\t", header=FALSE, skip=3, quote="")
tmp <- dat1[,1:110]
tmp <- tmp[, -seq(4, 110, by=2)]
tmp <- tmp[, -(1:2)]
test <- t(tmp)[1:46,]
filename <- paste("http://pubs.broadinstitute.org/mpr/projects/",
                  "Global_Cancer_Map/GCM_Test.cls", sep="")
test.classes <- read.table(filename, skip=2)+1
test.classes <- test.classes[test.classes!=15]
test.classes <- unlist(test.classes)
### pre-processing data
train[train < 20] <- 20
train[train > 16000] <- 16000
filter <- apply(train, 2, function(x) if((max(x)/min(x) > 5) && (max(x)-min(x) > 500))
                  return(1)
                else return(0))
train <- train[, filter==1]
train <- log10(train)
x <- train
meanx <- colMeans(x)
one <- rep(1,nrow(x))
normx <- sqrt(drop(one %*% (x^2)))
train <- scale(train, meanx, normx)
tmp <- cbind(train, train.classes)
tmp <- cbind(train, train.classes)
nx <- dim(tmp)[2]
a0 <- b0 <- 0
for(k in 1:length(table(train.classes))){
  tmp1 <- subset(tmp, tmp[,nx]==k)
  xc.bar <- colMeans(tmp1[,-nx]) ###average of gene j across class k
  xa.bar <- colMeans(tmp[,-nx]) ### average of gene j across all samples
  a0 <- a0 + dim(tmp1)[1] * ((xc.bar - xa.bar)^2)
  b0 <- b0 + colSums((tmp[,-nx] - xc.bar)^2)
}
bw <- a0/b0
npre <- nx - 1 ### use all genes and ignore bw values
bw1 <- order(bw, decreasing=TRUE)[1:npre]
train <- train[,bw1]
test[test < 20] <- 20
test[test > 16000] <- 16000
test <- test[, filter==1]
test <- log10(test)
test <- scale(test, meanx, normx)[, bw1]
colnames(train) <- paste("x", 1:dim(train)[2], sep="")
colnames(test) <- paste("x", 1:dim(test)[2], sep="")
###heat map
library(gplots)
# Refresh the plot window
graphics.off()
# Heatmap
heatmap.2(train, dendrogram = "none", col = redgreen, 
          key = TRUE, trace = "none", cexRow = 0.4, cexCol = 0.4, labCol = FALSE)
### 
Xtrain <- train
Xtest <- test
Ytrain <- train.classes
Ytest <- test.classes
### ISIS binomial classification
xtrain <- Xtrain[Ytrain == 9 | Ytrain == 5, ]
ytrain <- Ytrain[Ytrain == 9 | Ytrain == 5]
ytrain <- ifelse(ytrain == 9, 0, 1)
model0 <- SIS(xtrain, ytrain, family = "binomial", iter = TRUE,
              standardize = FALSE, type.measure = "deviance", 
              penalty = "lasso", tune = "cv", varISIS = "cons")
xtest <- Xtest[Ytest == 9 | Ytest == 5, ]
ypred0 <- predict(model0, xtest, type='class')
ytest <- Ytest[Ytest == 9 | Ytest == 5]
ytest <- ifelse(ytest == 9, 0, 1)
table(ypred0, ytest)
### SIS before LASSO,one vs all
coef0=NULL
for(i in 1:14){
  coef0=union(coef0,SIS(Xtrain, Ytrain==i, family = "binomial", iter = F, nsis = 50,
                  standardize = FALSE, type.measure = "deviance", 
                  penalty = "lasso", tune = "cv", varISIS = "cons")$sis.ix0
  )
}
cv.out0 <- cv.glmnet(Xtrain[,coef0], as.factor(Ytrain), type.measure = "class", 
                     lambda = 10^(seq(-5, -2, length.out = 50)), standardize = FALSE, 
                     family = "multinomial", alpha = 1, nfolds = 4)
ypred0 = predict(cv.out0, s = "lambda.min", newx = Xtest[,coef0], type = "class")
table(as.integer(ypred0),Ytest)
coefpred1 <- predict(cv.out0, s = "lambda.min", newx = Xtest, type = "coefficient")
coef1 <- unlist(sapply(coefpred1, function(u) u@i))
coef1 <- unique(coef1)
coef0 = coef0[coef1]
### one vs all SIS before LASSO
coef0=coef1=NULL
for(i in 1:14){
  coefpred0 = SIS(Xtrain, Ytrain==i, family = "binomial", iter = F, nsis = 50,
                        standardize = FALSE, type.measure = "deviance", 
                        penalty = "lasso", tune = "cv", varISIS = "cons")$sis.ix0
  coef0=union(coef0, coefpred0)
  cv.out1 = cv.glmnet(Xtrain[,coefpred0], Ytrain==i, type.measure = "class", 
            lambda = 10^(seq(-5, -2, length.out = 50)), standardize = FALSE, 
            family = "binomial", alpha = 1, nfolds = 4)
  coefpred1 <- predict(cv.out1, s = "lambda.min", type = "coefficient")
  coef1 <- union(coef1, coefpred0[coefpred1@i])
}
coef0 = coef1
### SIS one vs all + elastic.net
coef0=NULL
for(i in 1:14){
  coef0=union(coef0,SIS(Xtrain, Ytrain==i, family = "binomial", iter = F, nsis = 50,
                        standardize = FALSE, type.measure = "deviance", 
                        penalty = "lasso", tune = "cv", varISIS = "cons")$sis.ix0
  )
}
Data = data.frame(X = Xtrain[, coef0], Y = as.factor(Ytrain))
tuneGrid = data.frame(alpha = rep(seq(0.05, 0.95, length.out = 10), 50), 
        lambda = as.vector(sapply(10^(seq(-5, -2, length.out = 50)), rep, 10)))
cv_4 = trainControl(method = "cv", number = 4, search = "random")
hit_elnet = train(
  Y ~ ., data = Data, 
  method = "glmnet",
  trControl = cv_4, 
  tuneGrid = tuneGrid, 
  standardize = FALSE
)
cv.out0 <- glmnet(Xtrain[, coef0], as.factor(Ytrain), type.measure = "class", 
                     lambda = 10^(seq(-5, -2, length.out = 50)), standardize = FALSE, 
                     family = "multinomial", alpha = hit_elnet$bestTune$alpha)
ypred0 = predict(cv.out0, s = hit_elnet$bestTune$lambda, 
                 newx = Xtest[,coef0], type = "class")
table(as.integer(ypred0),Ytest)
coefpred1 <- predict(cv.out0, s = hit_elnet$bestTune$lambda, 
                     newx = Xtest, type = "coefficient")
coef1 <- unlist(sapply(coefpred1, function(u) u@i))
coef1 <- unique(coef1)
coef0 = coef0[coef1]
### SIS before LASSO, one vs one
ix=NULL
for(i in 1:13)
  for(j in (i+1):14){
    ix=union(ix,SIS(Xtrain[Ytrain == i | Ytrain == j,], 
                    Ytrain[Ytrain == i | Ytrain == j], 
                    family = "binomial", iter = F, nsis =4,
                    standardize = FALSE, type.measure = "deviance", 
                    penalty = "lasso", tune = "cv", varISIS = "cons")$sis.ix0
    )
  }
coef0 = ix
cv.out0 <- cv.glmnet(Xtrain[,coef0], as.factor(Ytrain), type.measure = "class", 
                     lambda = 10^(seq(-5, -2, length.out = 50)), standardize = FALSE, 
                     family = "multinomial", alpha = 1, nfolds = 4)
ypred0 = predict(cv.out0, s = "lambda.min", newx = Xtest[,coef0], type = "class")
table(as.integer(ypred0),Ytest)
coefpred1 <- predict(cv.out0, s = "lambda.min", newx = Xtest, type = "coefficient")
coef1 <- unlist(sapply(coefpred1, function(u) u@i))
coef1 <- unique(coef1)
coef0 = coef0[coef1]
### LASSO multinomial classification
cv.out <- cv.glmnet(Xtrain, as.factor(Ytrain), type.measure = "deviance", 
                    lambda = 10^(seq(-5, -2, length.out = 50)), standardize = FALSE, 
                    family = "multinomial", alpha = 1, nfolds = 4, 
                    type.multinomial = "ungrouped")
ypred1 = predict(cv.out, s = "lambda.min", newx = Xtest, type = "class")
table(as.integer(ypred1), Ytest)
coefpred1 <- predict(cv.out, s = "lambda.min", newx = Xtest, type = "coefficient")
coef1 <- unlist(sapply(coefpred1, function(u) u@i))
coef1 <- unique(coef1)
### elastic net
Data = data.frame(X = Xtrain, Y = as.factor(Ytrain))
tuneGrid = data.frame(alpha = rep(seq(0.05, 0.95, length.out = 10), 10), 
                      lambda = as.vector(sapply(10^(seq(-5, -2, length.out = 10)), 
                                                rep, 10)))
cv_4 = trainControl(method = "cv", number = 4, search = "random")
hit_elnet = train(
  Y ~ ., data = Data, 
  method = "glmnet",
  trControl = cv_4, 
  tuneGrid = tuneGrid, 
  standardize = FALSE
)
cv.out0 <- glmnet(Xtrain, as.factor(Ytrain), type.measure = "class", 
                  lambda = 10^(seq(-5, -2, length.out = 10)), standardize = FALSE, 
                  family = "multinomial", alpha = hit_elnet$bestTune$alpha)
ypred0 = predict(cv.out0, s = hit_elnet$bestTune$lambda, 
                 newx = Xtest[,coef0], type = "class")
table(as.integer(ypred0),Ytest)
coefpred1 <- predict(cv.out0, s = hit_elnet$bestTune$lambda, 
                     newx = Xtest, type = "coefficient")
coef1 <- unlist(sapply(coefpred1, function(u) u@i))
coef1 <- unique(coef1)
coef0 = coef1
### SVM
Xtrain.red0 <- Xtrain[, coef0]
Xtest.red0 <- Xtest[, coef0]

Data <- data.frame(X = rbind(Xtrain.red0, Xtest.red0), 
                   Y = c(Ytrain, Ytest))

svm.model0 <- svm(Xtrain.red0, as.factor(Ytrain), cost = 1e4, cross = 4)
svm.pred0 <- predict(svm.model0, Xtest.red0)
table(svm.pred0, Ytest)

tune.svm0 <- tune.svm(Xtrain.red0, as.factor(Ytrain), scale = FALSE, 
                      cost = 10^seq(-3, 5, 1), cross = 4)
svm.pred0 <- predict(tune.svm0$best.model, Xtest.red0)
table(svm.pred0, Ytest)

Xtrain.red <- Xtrain[, coefpred1$`1`@i]
Xtest.red <- Xtest[, coefpred1$`1`@i]
svm.model <- svm(Xtrain.red, as.factor(Ytrain))
svm.pred <- predict(svm.model, Xtest.red)
table(svm.pred, Ytest)
### Logistic
multi0 <- multinom(Y~., data = Data, 
                   subset = 1:length(Ytrain), MaxNWts = 10000)
multi.pred0 <- predict(multi0, newdata = Data[-(1:length(Ytrain)), ], type = "class")
table(as.integer(multi.pred0), Ytest)
### KNN
knn.pred0 <- knn(Xtrain.red0, Xtest.red0, Ytrain, k = 3)
table(knn.pred0, Ytest)
### Tree methods
### CART
Xtrain.red0 <- Xtrain[, coef0]
Xtest.red0 <- Xtest[, coef0]

tree0 <- tree(as.factor(Y) ~ ., data = Data, subset = 1:length(Ytrain))
tree.pred0 <- predict(tree0, newdata = Data[-(1:length(Ytrain)), ], type = "class")
table(tree.pred0, Ytest)

cv.tree0 <- cv.tree(tree0, FUN = prune.misclass, K = 4)
prune.tree0 <- prune.misclass(tree0, best = 14)
tree.cv.pred0 <- predict(prune.tree0, newdata = Data[-(1:length(Ytrain)), ], 
                         type = "class")
table(tree.cv.pred0, Ytest)
### Bagging and Random Forest
bag0 <- randomForest(Xtrain.red0, as.factor(Ytrain), mtry = ncol(Xtrain.red0))
bag.pred0 <- predict(bag0, newdata = Xtest.red0, type = "class")
table(bag.pred0, Ytest)

rf0 <- randomForest(Xtrain.red0, as.factor(Ytrain), 
                     mtry = floor(sqrt(ncol(Xtrain.red0))))
rf.pred0 <- predict(rf0, newdata = Xtest.red0, type = "class")
table(rf.pred0, Ytest)
### Boosting method
gbm0 <- gbm(Y ~ ., data = Data[1:length(Ytrain), ], distribution = "multinomial", 
            n.trees = 1000, shrinkage = 0.1, interaction.depth = 5, 
            bag.fraction = 0.5, train.fraction = 0.5, n.minobsinnode = 10, 
            cv.folds = 0, keep.data = TRUE, verbose = TRUE)
#best.iter <- gbm.perf(gbm0, method = "cv")
gbm.pred0 <- predict(gbm0, newdata = Data[-(1:length(Ytrain)), ], type = "response")
table(apply(gbm.pred0, 1, which.max), Ytest)
### adaboost(dont work)
ada0 <- boosting(as.factor(Y) ~ ., data = Data[1:length(Ytrain), ],
                 boos = TRUE, mfinal = 10, maxdepth = 5)
### RDA(dont work)
rda0 <- rda(Y ~ ., data = Data[1:length(Ytrain), ], gamma = 0.05, lambda = 0.2)
rda.pred0 <- predict(rda0, Data[-(1:length(Ytrain)), ], type = "class")
table(rda.pred0$class, Ytest)

### plot for different classes
# Heatmap
heatmap.2(Xtrain[Ytrain == 6 | Ytrain == 1, coef0], 
          dendrogram = "col", col = redgreen, key = TRUE, trace = "none", 
          cexRow = 0.4, cexCol = 0.4, labCol = FALSE)
# gene centroid plot
par(mfrow = c(1, 5))
for(k in 1:14) 
{
  cv.out1 = cv.glmnet(Xtrain[,coef0], Ytrain==k, type.measure = "class", 
                      lambda = 10^(seq(-5, -2, length.out = 50)), 
                      standardize = FALSE, 
                      family = "binomial", alpha = 1, nfolds = 4)
  coefpred1 <- predict(cv.out1, s = "lambda.min", type = "coefficient")
  coef01 <- coef0[coefpred1@i]
  cent <- colMeans(Xtrain[Ytrain == k, coef0])
  plot(0, 700, 'n', ylim = c(0, 700), ylab = "selected variables", 
       xlab = paste("class", k, sep = ""))
  for(i in 1:length(coef0)) {
    lines(c(0,cent[i]*4e1), c(i,i))
  }
  for(i in 2:length(coef01)) {
    lines(c(0,cent[coef01[i]]*4e1), c(coef01[i],coef01[i]), 
          lwd = 2, col = 'blue')
  }
}

