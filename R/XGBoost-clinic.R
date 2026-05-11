# ================================
setwd("D:/desktop/train/model/XGBoost/clinic")
library(mlr3learners)
library(mlr3verse)
library(mlr3proba)
library(dplyr)
lrn("surv.xgboost")
# ================================
data <- read.csv("clinical+data.csv", sep = ",", header = TRUE)
data$status <- ifelse(data$status == 'Alive', 0, 1)
data$stage <- as.numeric(as.factor(data$stage))
data <- data[, c("time", "status", "stage", "age")]
str(data)
# ================================
task <- as_task_surv(data, 
                     time = "time",
                     event = "status", 
                     type = "right")
task
# ================================
learner <- lrn("surv.xgboost")
learner
# ================================
set.seed(376)
split <- partition(task, ratio = 0.7)
# ================================
learner$train(task, row_ids = split$train)
fit <- learner$model
fit
# ================================
predictions_train <- learner$predict(task, row_ids = split$train)
predictions_train
cindex_train <- predictions_train$score(msr("surv.cindex"))
cat("Train C-index:", cindex_train, "\n")
# ================================
predictions_valid <- learner$predict(task, row_ids = split$test)
predictions_valid
cindex_valid <- predictions_valid$score(msr("surv.cindex"))
cat("Validation C-index:", cindex_valid, "\n")
# ================================
test_data <- read.csv("test+data.csv", sep = ",", header = TRUE)
test_data$status <- ifelse(test_data$status == 'Alive', 0, 1)
test_data$stage <- as.numeric(as.factor(test_data$stage))
test_data <- test_data[, c("time", "status", "stage", "age")]
task_test <- as_task_surv(test_data, 
                          time = "time", 
                          event = "status", 
                          type = "right")
# ================================
predictions_test <- learner$predict_newdata(test_data)
predictions_test
cindex_test <- predictions_test$score(msr("surv.cindex"))
cat("Test C-index:", cindex_test, "\n")
# ================================
result_summary <- data.frame(
  Dataset = c("Training", "Validation", "Test"),
  C_index = c(cindex_train, cindex_valid, cindex_test)
)
print(result_summary)
write.csv(result_summary, "xgboost_cindex_results.csv", row.names = FALSE)