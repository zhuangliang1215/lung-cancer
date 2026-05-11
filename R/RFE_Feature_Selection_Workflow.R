# ============================================================================
# RFE_Feature_Selection_Workflow.R
# 
# Description: Recursive Feature Elimination (RFE) for multiple algorithms
#              Includes Linear Regression, Random Forest, and SVM-RFE
# 
# Algorithms: 
#   1. Linear Regression RFE (lmFuncs)
#   2. Random Forest RFE (rfFuncs) 
#   3. SVM-RFE (caretFuncs with svmRadial)
# 
# Output: Selected features, variable importance plots, performance metrics
# ============================================================================

# ============================================================
# Part 1: Load required packages
# ============================================================
library(caret)
library(psych)
library(mlbench)
library(rio)
library(ggplot2)
library(future)

# ============================================================
# Part 2: BloodBrain data example (demonstration)
# ============================================================

# Load and explore data
data(BloodBrain)
psych::headTail(bbbDescr)
head(logBBB)

# Remove near-zero variance predictors
x <- scale(bbbDescr[, -nearZeroVar(bbbDescr)])

# Remove highly correlated predictors
x <- x[, -findCorrelation(cor(x), 0.8)]
x <- as.data.frame(x, stringsAsFactors = TRUE)
psych::headTail(x)

# ============================================================
# Part 3: Linear Regression RFE (lmFuncs)
# ============================================================

# Setup parallel processing
plan("multisession", workers = 8)

# RFE control settings
set.seed(1)
rfeControl_lm <- rfeControl(
  functions = lmFuncs,
  method = "cv",           # Cross-validation
  saveDetails = TRUE,
  number = 5,              # 5-fold CV
  allowParallel = TRUE
)

# Run RFE
set.seed(1)
lmProfile <- rfe(
  x, logBBB,
  sizes = c(2:25, 30, 35, 40, 45, 50, 55, 60, 65),
  rfeControl = rfeControl_lm
)

# View results
lmProfile
predictors(lmProfile)

# Plot performance metrics
ggplot(data = lmProfile, metric = "RMSE") + theme_bw()
ggplot(data = lmProfile, metric = "MAE") + theme_bw()

# Variable importance
varImp(lmProfile)

# Create importance plot
varimp_data_lm <- data.frame(
  feature = row.names(varImp(lmProfile))[1:22],
  importance = varImp(lmProfile)[1:22, 1]
)

ggplot(data = varimp_data_lm, aes(x = reorder(feature, -importance), 
                                  y = importance, fill = feature)) +
  geom_bar(stat = "identity") + 
  labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust = 1.6, color = "white", size = 4) + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# Part 4: Random Forest RFE (rfFuncs)
# ============================================================

# RFE control for Random Forest
set.seed(1)
rfeControl_rf <- rfeControl(
  functions = rfFuncs,
  method = "cv",
  number = 5
)

# Setup parallel processing
plan("multisession", workers = 8)

# Run RFE with Random Forest
set.seed(1)
rfProfile <- rfe(
  x, logBBB,
  sizes = c(2:25, 30, 35, 40, 45, 50, 55, 60, 65),
  rfeControl = rfeControl_rf,
  allowParallel = TRUE
)

# View results
rfProfile

# ============================================================
# Part 5: SVM-RFE (caretFuncs with svmRadial)
# ============================================================

# RFE control for SVM
set.seed(1)
rfeControl_svm <- rfeControl(
  functions = caretFuncs,
  method = "cv",
  number = 5
)

# Setup parallel processing
plan("multisession", workers = 8)

# Run SVM-RFE
set.seed(1)
svmProfile <- rfe(
  x, logBBB,
  sizes = c(2:25, 30, 35, 40, 45, 50, 55, 60, 65),
  method = "svmRadial",
  rfeControl = rfeControl_svm,
  allowParallel = TRUE
)

# View results
svmProfile

# ============================================================
# Part 6: Python-style SVM-RFE (sklearn demonstration)
# ============================================================

# Note: This section is Python code, not R
# from sklearn.datasets import make_friedman1
# from sklearn.feature_selection import RFE
# from sklearn.svm import SVR
# 
# X, y = make_friedman1(n_samples=50, n_features=95, random_state=0)
# estimator = SVR(kernel="linear")
# selector = RFE(estimator, n_features_to_select=50, step=1)
# selector = selector.fit(X, y)
# 
# print(selector.support_)  # Selected features (True = selected)
# print(selector.ranking_)  # Feature ranking (lower = better)
# print(selector.n_features_)  # Number of selected features

# ============================================================
# Part 7: SVM-RFE for classification (binary outcome)
# ============================================================

set.seed(7)

# Load data from clipboard (or use your own data)
# predcm <- rio::import("clipboard", header = TRUE)

# Example using PimaIndiansDiabetes dataset
data(PimaIndiansDiabetes)
predcm <- PimaIndiansDiabetes

# RFE control
control <- rfeControl(
  functions = caretFuncs,
  method = "cv",
  number = 5
)

# Run SVM-RFE for classification
results <- rfe(
  predcm[, 1:7],
  as.factor(predcm[, 8]),  # Outcome must be factor
  sizes = c(1:7),
  rfeControl = control,
  method = "svmRadial"
)

# View results
print(results)
predictors(results)
plot(results, type = c("g", "o"))

# ============================================================
# Part 8: RFE on custom dataset (datExpr.csv)
# ============================================================

# Set working directory and read data
setwd("D:/desktop/FeatureSelection")
rt <- read.csv("datExpr.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Convert status: Dead = 1, Alive = 0
rt$status <- ifelse(rt$status == "Dead", 1, 0)

# ============================================================
# Part 9: Linear Regression RFE on custom data
# ============================================================

# Setup parallel processing
plan("multisession", workers = 8)

# RFE control
set.seed(1)
rfeControl_lm_custom <- rfeControl(
  functions = lmFuncs,
  method = "cv",
  saveDetails = TRUE,
  number = 5,
  allowParallel = TRUE
)

# Run RFE
set.seed(1)
lmProfile_custom <- rfe(
  rt[, 3:82],   # Predictors (columns 3-82)
  rt[, 1],      # Outcome (column 1: status)
  sizes = c(3:82),
  rfeControl = rfeControl_lm_custom
)

# View results
lmProfile_custom
predictors(lmProfile_custom)

# Plot performance
ggplot(data = lmProfile_custom, metric = "RMSE") + theme_bw()
ggplot(data = lmProfile_custom, metric = "MAE") + theme_bw()

# Variable importance plot
varimp_data_lm_custom <- data.frame(
  feature = row.names(varImp(lmProfile_custom))[1:95],
  importance = varImp(lmProfile_custom)[1:95, 1]
)

ggplot(data = varimp_data_lm_custom, aes(x = reorder(feature, -importance), 
                                         y = importance, fill = feature)) +
  geom_bar(stat = "identity") + 
  labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust = 1.6, color = "white", size = 4) + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# Part 10: SVM-RFE on custom data (Method 1)
# ============================================================

set.seed(1)

# RFE control
control_custom <- rfeControl(
  functions = caretFuncs,
  method = "cv",
  number = 5
)

# Run SVM-RFE
results_custom <- rfe(
  rt[, 3:82],
  rt[, 1],
  sizes = c(3:82),
  rfeControl = control_custom,
  method = "svmRadial"
)

# View results
print(results_custom)
predictors(results_custom)
plot(results_custom, type = c("g", "o"))

# Performance plots
ggplot(data = results_custom, metric = "RMSE") + theme_bw()
ggplot(data = results_custom, metric = "MAE") + theme_bw()

# Variable importance
varimp_data_svm1 <- data.frame(
  feature = row.names(varImp(results_custom))[1:95],
  importance = varImp(results_custom)[1:95, 1]
)

ggplot(data = varimp_data_svm1, aes(x = reorder(feature, -importance), 
                                    y = importance, fill = feature)) +
  geom_bar(stat = "identity") + 
  labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust = 1.6, color = "white", size = 4) + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# Part 11: SVM-RFE on custom data (Method 2 - with allowParallel)
# ============================================================

set.seed(1)
rfeControl_svm2 <- rfeControl(
  functions = caretFuncs,
  method = "cv",
  number = 5
)

# Setup parallel processing
plan("multisession", workers = 8)

set.seed(1)
svmProfile_custom <- rfe(
  rt[, 3:82],
  rt[, 1],
  sizes = c(3:82),
  method = "svmRadial",
  rfeControl = rfeControl_svm2,
  allowParallel = TRUE
)

# View results
svmProfile_custom
predictors(svmProfile_custom)

# Performance plots
ggplot(data = svmProfile_custom, metric = "RMSE") + theme_bw()
ggplot(data = svmProfile_custom, metric = "MAE") + theme_bw()

# Variable importance
varimp_data_svm2 <- data.frame(
  feature = row.names(varImp(svmProfile_custom))[1:95],
  importance = varImp(svmProfile_custom)[1:95, 1]
)

ggplot(data = varimp_data_svm2, aes(x = reorder(feature, -importance), 
                                    y = importance, fill = feature)) +
  geom_bar(stat = "identity") + 
  labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust = 1.6, color = "white", size = 4) + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# Part 12: Random Forest RFE on custom data
# ============================================================

set.seed(1)
rfeControl_rf_custom <- rfeControl(
  functions = rfFuncs,
  method = "cv",
  number = 5
)

# Setup parallel processing
plan("multisession", workers = 8)

set.seed(1)
rfProfile_custom <- rfe(
  rt[, 3:82],
  rt[, 1],
  sizes = c(3:82),
  rfeControl = rfeControl_rf_custom,
  allowParallel = TRUE
)

# View results
rfProfile_custom
predictors(rfProfile_custom)

# Performance plots
ggplot(data = rfProfile_custom, metric = "RMSE") + theme_bw()
ggplot(data = rfProfile_custom, metric = "MAE") + theme_bw()

# Variable importance
varimp_data_rf <- data.frame(
  feature = row.names(varImp(rfProfile_custom))[1:95],
  importance = varImp(rfProfile_custom)[1:95, 1]
)

ggplot(data = varimp_data_rf, aes(x = reorder(feature, -importance), 
                                  y = importance, fill = feature)) +
  geom_bar(stat = "identity") + 
  labs(x = "Features", y = "Variable Importance") + 
  geom_text(aes(label = round(importance, 2)), vjust = 1.6, color = "white", size = 4) + 
  theme_bw() + 
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================
# Part 13: Summary report
# ============================================================

cat("\n", rep("=", 60), "\n", sep = "")
cat("RFE FEATURE SELECTION WORKFLOW COMPLETE\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Algorithms executed:\n")
cat("  1. Linear Regression RFE (lmFuncs)\n")
cat("  2. Random Forest RFE (rfFuncs)\n")
cat("  3. SVM-RFE (caretFuncs with svmRadial)\n\n")

cat("Analysis completed at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 60), "\n", sep = "")