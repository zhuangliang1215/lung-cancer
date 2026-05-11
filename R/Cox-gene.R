# ================================
setwd("D:/desktop/train/model/Cox/gene")
library(rms)
library(survival)
library(survminer)
# ================================
bc <- read.csv("D:/desktop/train/model/Cox/gene/gene+data.csv", 
               sep = ",", header = TRUE)
bc <- bc[, c("time", "status", "riskScore")]
bc <- na.omit(bc)
# ================================
set.seed(376) 
train_indices <- sample(1:nrow(bc), size = 0.6 * nrow(bc), replace = FALSE)
train_data <- bc[train_indices, ]
validation_data <- bc[-train_indices, ]
# ================================
cox_model <- coxph(Surv(time, status) ~ riskScore, data = train_data)
cat("\n========== Cox ==========\n")

print(summary(cox_model))
# ================================
validation_data$predicted_risk <- predict(cox_model, newdata = validation_data, type = "risk")
c_index_val <- survConcordance(Surv(time, status) ~ predicted_risk, data = validation_data)
cat("\n========== Validation data） ==========\n")
cat(paste("C-index:", round(c_index_val$concordance, 4), "\n"))
# ================================
test_data <- read.csv("test+data.csv", header = TRUE, row.names = 1, check.names = FALSE)
test_data$status <- ifelse(test_data$status == 'Alive', 0, 1)
test_data <- test_data[, c("time", "status", "riskScore")]
# ================================
test_data$predicted_risk <- predict(cox_model, newdata = test_data, type = "risk")
c_index_test <- survConcordance(Surv(time, status) ~ predicted_risk, data = test_data)
cat("\n========== Test data ==========\n")
cat(paste("C-index:", round(c_index_test$concordance, 4), "\n"))
# ================================
cat("\n========== Data Sum ==========\n")
cat(paste("Train C-index:", round(summary(cox_model)$concordance[1], 4), "\n"))
cat(paste("Validation C-index:", round(c_index_val$concordance, 4), "\n"))
cat(paste("Test C-index:", round(c_index_test$concordance, 4), "\n"))