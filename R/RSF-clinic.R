# ================================
library(survival)
library(survminer)
library(randomForestSRC)
library(timeROC)
library(data.table)
library(tidyr)
library(dplyr)
library(tibble)
library(caret)
library(ranger)
library(ggRandomForests)
library(rsample)
library(purrr)
# ================================
setwd("D:/desktop/train/model/RSF/clinic")
lung <- read.csv("clinical+data.csv", header = TRUE, row.names = 1, check.names = FALSE)
lung$status <- ifelse(lung$status == 'Alive', 0, 1)
lung$stage <- as.numeric(as.factor(lung$stage))
lung <- lung[, c("time", "status", "stage", "age")]
# ================================
set.seed(376)  
train_indices <- sample(1:nrow(lung), 0.7 * nrow(lung))
train_data <- lung[train_indices, ]
validation_data <- lung[-train_indices, ]
# ================================
rf_model <- rfsrc(Surv(time, status) ~ age + stage, data = train_data)
# ================================
train_pred <- predict(rf_model, train_data)
cindex_train <- survConcordance(Surv(train_data$time, train_data$status) ~ train_pred$predicted)
cat("Training set C-index:", round(cindex_train$concordance, 4), "\n")
# ================================
validation_pred <- predict(rf_model, validation_data)
cindex_validation <- survConcordance(Surv(validation_data$time, validation_data$status) ~ validation_pred$predicted)
cat("Validation set C-index:", round(cindex_validation$concordance, 4), "\n")
# ================================
test_data <- read.csv("test+data.csv", header = TRUE, row.names = 1, check.names = FALSE)
test_data$status <- ifelse(test_data$status == 'Alive', 0, 1)
test_data$stage <- as.numeric(as.factor(test_data$stage))
test_data <- test_data[, c("time", "status", "stage", "age")]
test_pred <- predict(rf_model, test_data)
cindex_test <- survConcordance(Surv(test_data$time, test_data$status) ~ test_pred$predicted)
cat("Test set C-index:", round(cindex_test$concordance, 4), "\n")
print(cindex_test)