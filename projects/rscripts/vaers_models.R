# Daniel Wang
# Predictive Modelling Analysis
# 4/28/26

library(tidyverse)
library(caret) # for QDA and RF
library(glmnet) # for Ridge and Logistic regression
library(gam) # for GAMS

setwd("projects/data")
load("vaers.2018-2021.RData")

# ______________________________________________________________________________
# Transforming predictors / converting column types / Imputation
vaers_data <- vaers_dat %>%
  mutate(across(c(3,5,7:18), as.factor)) %>%
  mutate(AGE_YRS = ifelse(is.na(AGE_YRS), mean(AGE_YRS, na.rm=T), AGE_YRS)) %>%
  mutate(NUMDAYS = ifelse(is.na(NUMDAYS), median(NUMDAYS, na.rm=T), NUMDAYS)) %>%
  select(-c(19,20)) # dropping two extremely rare vaccine types
# note: Imputing before adding polynomial term, so that the distribution of the polynomial term is consistent
vaers_data$AGE_YRS2 <- vaers_data$AGE_YRS**2
vaers_data <- vaers_data %>%
  relocate(AGE_YRS2, .after = AGE_YRS)
# ______________________________________________________________________________


# ______________________________________________________________________________
# Sample for analysis
set.seed(208)

samp <- sample(1:nrow(vaers_data), size = 10000, replace = FALSE)
vaer.samp <- vaers_data[samp, -c(1,2)]

# training and test sets
index <- sample(1:nrow(vaer.samp), round(0.7*nrow(vaer.samp)), replace = F)

trainer <- vaer.samp[index,]
test <- vaer.samp[-index,]
# ______________________________________________________________________________


# ______________________________________________________________________________
# Logistic Regression

# logistic fit on training data
log.fit <- glm(high_severity ~ ., data=trainer, family = "binomial")

# predicting on test set
log.pred <- predict(log.fit, newdata = test, type = "response")
log.pred <- ifelse(log.pred > 0.5, 1, 0)

# evaluating model
pred.table <- table(Predicted = log.pred, True = test$high_severity)
confusionMatrix(pred.table)
# high accuracy, low specificity rate
# ______________________________________________________________________________


# ______________________________________________________________________________
# Linear/quadratic Discriminant Analysis

# training method
tr <- trainControl(method = "repeatedcv", number = 10, repeats = 20)

lda.fit <- train(high_severity ~ .,
                 data = trainer, method="lda",
                 trControl = tr)

lda.pred <- predict(lda.fit, newdata = test, type = "raw")

pred.table3 <- table(Predicted = lda.pred, True = test$high_severity)
confusionMatrix(pred.table3)
# high accuracy, better specificity rate than logistic regression
# so far best balance of accuracy and specificity


qda.fit <- train(high_severity ~ .,
                 data = trainer, method="qda",
                 trControl = tr)

qda.pred <- predict(qda.fit, newdata = test, type = "raw")

pred.table4 <- table(Predicted = qda.pred, True = test$high_severity)
confusionMatrix(pred.table4)
# moderate accuracy, higher specificity rate
# ______________________________________________________________________________


# ______________________________________________________________________________
# Random Forest
tr2 <- trainControl(method = "oob", number = 3000)
tune.m <- data.frame(mtry = 1:(ncol(trainer)-1))

# Bagged random forest, might take a minute to tun
rf.fit <- train(high_severity ~ ., method="rf",
                data = trainer, tuneGrid=tune.m,
                trControl=tr2)

rf.pred <- predict(rf.fit, newdata = test, type = "raw")

pred.table5 <- table(Predicted = rf.pred, True = test$high_severity)
confusionMatrix(pred.table5)
# highest accuracy, though only slightly
# slightly higher specificity than logistic regression, still lower than LDA
# ______________________________________________________________________________


# ______________________________________________________________________________
# Regularization with Ridge regression
# Not included in analysis

# response vector
Y <- vaer.samp$high_severity
# model matrix
model.mat <- makeX(vaer.samp[,-1])

# grid of possible lambda values
grid.1 <- 10**seq(-2,10,length=100)

# 10 fold CV to optimize lambda
cv.out.1 <- cv.glmnet(x=model.mat, y=Y, family="binomial",
                      alpha=0, lambda = grid.1)

bestlam <- data.frame(min = cv.out.1$lambda.min,
                      one_se = cv.out.1$lambda.1se)

# refit ridge regression
# 1se means 1-SE
ridge.min <- predict(cv.out.1$glmnet.fit,
                     type = "coefficients",
                     s = bestlam$min)
ridge_1se <- predict(cv.out.1$glmnet.fit,
                     type = "coefficients",
                     s = bestlam$one_se)

# parameter estimates
betahat <- cbind(ridge.min, ridge_1se)
colnames(betahat) <- c("min", "1se")
rownames(betahat) <- attributes(betahat)$Dimnames[[1]]
betahat
# ______________________________________________________________________________