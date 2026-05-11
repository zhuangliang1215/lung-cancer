setwd("D:/desktop/train/model/XGBoost/xgboost_model")
library(mlr3learners)
library(mlr3verse)
library(mlr3proba)
library(dplyr)
lrn("surv.xgboost")
# ================================
data <- read.csv("D:/desktop/train/model/XGBoost/xgboost_model/gene+data.csv", 
                 sep = ",", header = TRUE)
data <- data[, c("time", "status", "riskScore")]
task <- as_task_surv(data, 
                     time = "time",
                     event = "status", 
                     type = "right")
task
# ================================
learner <- lrn("surv.xgboost")
learner

set.seed(376)
split <- partition(task, ratio = 0.7)
learner$train(task, row_ids = split$train)
fit <- learner$model
fit
predictions_train <- learner$predict(task, row_ids = split$train)
predictions_train
cindex_train <- predictions_train$score(msr("surv.cindex"))
cat("Train C-index:", cindex_train, "\n")

predictions_test <- learner$predict(task, row_ids = split$test)
predictions_test

cindex_test <- predictions_test$score(msr("surv.cindex"))
cat("Validation C-index:", cindex_test, "\n")
# ================================
test_data <- read.csv("D:/desktop/train/model/XGBoost/xgboost_model/test+data.csv", 
                      sep = ",", header = TRUE)
# Alive -> 0, Dead -> 1
if ("status" %in% colnames(test_data)) {
  test_data$status <- ifelse(test_data$status == 'Alive', 0, 1)
}
test_data <- test_data[, c("time", "status", "riskScore")]
task_test <- as_task_surv(test_data, 
                          time = "time", 
                          event = "status", 
                          type = "right")
predictions_external <- learner$predict_newdata(test_data)
predictions_external
cindex_external <- predictions_external$score(msr("surv.cindex"))
cat("Test C-index:", cindex_external, "\n")
# ================================
cat("\n========== Data Sum ==========\n")
cat("Train C-index:", round(cindex_train, 4), "\n")
cat("Validation C-index:", round(cindex_test, 4), "\n")
cat("Test C-index:", round(cindex_external, 4), "\n")
cat("==================================\n")