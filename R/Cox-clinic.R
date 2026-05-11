# ================================
setwd("D:/desktop/train/model/Cox/clinic")
library(rms)
library(survival)
library(survminer)
# ================================
bc <- read.csv("clinical+data.csv", sep = ",", header = TRUE)
bc <- bc[, c("time", "status", "stage", "age")]
bc$status <- ifelse(bc$status == 'Alive', 0, 1)
bc$stage <- as.numeric(as.factor(bc$stage))
# ================================
set.seed(376) 
train_indices <- sample(1:nrow(bc), 0.7 * nrow(bc), replace = FALSE)
train_data <- bc[train_indices, ]
validation_data <- bc[-train_indices, ]
# ================================
cox_model <- coxph(Surv(time, status) ~ age + stage, data = train_data)
summary(cox_model)
# ================================
validation_data$predicted_risk <- predict(cox_model, newdata = validation_data, type = "risk")
c_index_val <- survConcordance(Surv(time, status) ~ predicted_risk, data = validation_data)
print(paste("Validation C-index:", round(c_index_val$concordance, 4)))
# ================================
test_data <- read.csv("test+data.csv", header = TRUE, row.names = 1, check.names = FALSE)
test_data$status <- ifelse(test_data$status == 'Alive', 0, 1)
test_data$stage <- as.numeric(as.factor(test_data$stage))
test_data <- test_data[, c("time", "status", "stage", "age")]
test_data$predicted_risk <- predict(cox_model, newdata = test_data, type = "risk")
c_index_test <- survConcordance(Surv(time, status) ~ predicted_risk, data = test_data)
print(paste("Test C-index:", round(c_index_test$concordance, 4)))