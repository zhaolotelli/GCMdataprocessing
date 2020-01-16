library(plsgenomics)
data(SRBCT)
Xtrain <- SRBCT$X[1:63, ]
ytrain <- as.matrix(SRBCT$Y[1:63])
Xtest <- SRBCT$X[64:83, ]
ytest <- as.matrix(SRBCT$Y[64:83])

library(SIS)
meanx <- colMeans(Xtrain)
one <- rep(1, nrow(Xtrain))
normx <- sqrt(drop(one %*% (Xtrain^2)))
xtrain <- scale(Xtrain, meanx, normx)
model0 <- SIS(xtrain, ytrain, family = "gaussian", iter = TRUE, nsis = 15, 
              penalty = "lasso", nfolds = 4, tune = "cv")
xtest <- scale(Xtest, meanx, normx)
ypred0 <- predict(model0, xtest, type='response')
ypred0 <- round(ypred0)
ypred0[ypred0 >= 5] = 4
ypred0[ypred0 <= 0] = 1
table(ypred0, ytest)

Xtrain <- Xtrain[ytrain == 3 | ytrain == 4, ]
meanx <- colMeans(Xtrain)
one <- rep(1, nrow(Xtrain))
normx <- sqrt(drop(one %*% (Xtrain^2)))
xtrain <- scale(Xtrain, meanx, normx)
ytrain <- ytrain[ytrain == 3 | ytrain == 4]
ytrain <- ifelse(ytrain == 3, 0, 1)
model0 <- SIS(xtrain, ytrain, family = "binomial", iter = TRUE, nsis = 15, 
              penalty = "lasso", nfolds = 4, tune = "cv")
Xtest <- Xtest[ytest == 3 | ytest == 4, ]
xtest <- scale(Xtest, meanx, normx)
ypred0 <- predict(model0, xtest, type='class')
ytest <- ytest[ytest == 3 | ytest == 4]
ytest <- ifelse(ytest == 3, 0, 1)
table(ypred0, ytest)

library(glmnet)
model1 <- tune.fit(xtrain, ytrain, family = "gaussian", 
                   penalty = "lasso", nfolds = 4, tune = "cv")
bestlam <- model1$lambda
out <- glmnet(xtrain, ytrain, alpha = 1, lambda = bestlam)
ypred1 = predict(out, s = bestlam, newx = xtest)
ypred1 <- round(ypred1)
ypred1[ypred1 >= 5] = 4
ypred1[ypred1 <= 0] = 1
table(ypred1, ytest)
