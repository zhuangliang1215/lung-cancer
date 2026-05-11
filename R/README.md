# Code for "A Multimodal Ensemble Learning Model for NSCLC Survival Prediction"

This repository contains all R scripts used in this study.

## Scripts and their purposes

| File | Description |
|------|-------------|
| `TCGA_Clinical_Data_Processing.R` | Processing TCGA clinical data |
| `RFE_Feature_Selection_Workflow.R` | SVM-RFE feature selection |
| `LASSO_Cox_Survival_Analysis.R` | Lasso-Cox regression analysis |
| `GRS_Cox_univariate_multivariate.R` | Gene risk score calculation |
| `KM_curve_analysis.R` | Kaplan-Meier survival analysis |
| `ROC_curve_analysis_1_3_5_years.R` | ROC curve analysis |
| `RSF-clinic.R`, `RSF-gene.R`, `RSF-pet.R` | Random Survival Forest models |
| `XGBoost-clinic.R`, `XGBoost-gene.R` | XGBoost models |
| `Cox-clinic.R`, `Cox-gene.R`, `Cox-pet.R` | Cox proportional hazards models |

## Requirements
- R version 4.2.5 or higher
- Required packages: survival, randomForestSRC, xgboost, glmnet, DESeq2, pROC

## Reproducibility
Run the scripts in the order listed above to reproduce the main results reported in the paper.
