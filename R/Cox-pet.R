setwd("D:/desktop/train/model/Cox/pet")
library(rms)
library(survival)
library(survminer)
# ================================
bc <- read.csv("D:/desktop/train/model/Cox/pet/pet+data.csv", sep = ",", header = TRUE)
bc <- bc[, c("time", "status", "SUVmean", "SUVmax")]
bc$status <- ifelse(bc$status == 'Alive', 0, 1)
bc <- na.omit(bc)
# ================================
set.seed(376)  
train_indices <- sample(1:nrow(bc), 0.6 * nrow(bc), replace = FALSE)
train_data <- bc[train_indices, ]
validation_data <- bc[-train_indices, ]
# ================================
cox_model <- coxph(Surv(time, status) ~ SUVmean + SUVmax, data = train_data)
summary(cox_model)
# ================================
validation_data$predicted_risk <- predict(cox_model, newdata = validation_data, type = "risk")
c_index_val <- survConcordance(Surv(time, status) ~ predicted_risk, data = validation_data)
print(paste("Validation set C-index:", round(c_index_val$concordance, 4)))
# ================================
test_data <- read.csv("test+data.csv", header = TRUE, row.names = 1, check.names = FALSE)
test_data$status <- ifelse(test_data$status == 'Alive', 0, 1)
test_data <- test_data[, c("time", "status", "SUVmax", "SUVmean")]
test_data <- na.omit(test_data)
test_data$predicted_risk <- predict(cox_model, newdata = test_data, type = "risk")

c_index_test <- survConcordance(Surv(time, status) ~ predicted_risk, data = test_data)
print(paste("Test set C-index:", round(c_index_test$concordance, 4)))