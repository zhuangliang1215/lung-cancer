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
| `RSF-clinic.R` | Random Survival Forest (clinical modality) |
| `RSF-gene.R` | Random Survival Forest (genomic modality) |
| `RSF-pet.R` | Random Survival Forest (PET/CT modality) |
| `XGBoost-clinic.R` | XGBoost (clinical modality) |
| `XGBoost-gene.R` | XGBoost (genomic modality) |
| `XGBoost-pet.R` | XGBoost (PET/CT modality) |
| `Cox-clinic.R` | Cox regression (clinical modality) |
| `Cox-gene.R` | Cox regression (genomic modality) |
| `Cox-pet.R` | Cox regression (PET/CT modality) |

## Requirements

- R version 4.2.5 or higher
- Required packages: survival, randomForestSRC, xgboost, glmnet, DESeq2, pROC

## Reproducibility

Run the scripts in the order listed above to reproduce the main results reported in the paper.
