# ================================
setwd("D:/desktop/train/model/XGBoost/pet")

library(mlr3learners)
library(mlr3verse)
library(mlr3proba)
library(dplyr)

lrn("surv.xgboost")

# ================================
data <- read.csv("D:/desktop/train/model/XGBoost/pet/pet+data.csv", 
                 sep = ',', header = TRUE)

data$status <- ifelse(data$status == 'Alive', 0, 1)
data <- data[, c("time", "status", "SUVmax", "SUVmean")]
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
train_predictions <- learner$predict(task, row_ids = split$train)
train_predictions
train_cindex <- train_predictions$score(msr("surv.cindex"))
cat("Training set C-index:", train_cindex, "\n")
# ================================
test_predictions <- learner$predict(task, row_ids = split$test)
test_predictions
test_cindex <- test_predictions$score(msr("surv.cindex"))
cat("Validation set C-index:", test_cindex, "\n")
# ================================
test1 <- read.csv("D:/desktop/train/model/XGBoost/pet/test+data.csv", 
                  sep = ',', header = TRUE)
test1$status <- ifelse(test1$status == 'Alive', 0, 1)
test1 <- test1[, c("time", "status", "SUVmax", "SUVmean")]
task_test1 <- as_task_surv(test1, 
                           time = "time", 
                           event = "status", 
                           type = "right")

external_predictions <- learner$predict(task_test1)
external_predictions

external_cindex <- external_predictions$score(msr("surv.cindex"))
cat("External test set C-index:", external_cindex, "\n")
# ================================
performance_summary <- data.frame(
  Dataset = c("Training", "Validation", "External Test"),
  C_index = c(train_cindex, test_cindex, external_cindex)
)
print(performance_summary)
write.csv(performance_summary, "xgboost_survival_performance.csv", row.names = FALSE)