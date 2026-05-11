# ============================================================================
# TCGA_Clinical_Data_Processing_Cox.R
# 
# Description: Complete workflow for TCGA clinical data processing and Cox regression
#              Includes data extraction, preprocessing, Table 1 generation,
#              univariate and multivariate Cox analysis, and forest plot
# 
# Output: Processed data (CSV/Rdata), clinical tables, Cox results, forest plots
# ============================================================================

# ============================================================
# Part 1: Load packages
# ============================================================
library(R.utils)      
library(data.table)   
library(plyr)        
library(survival)    
library(survminer)    
library(tableone)    
library(forestplot)  

# ============================================================
# Part 2: Set working directory and read data
# ============================================================
setwd("D:/desktop/NewFolder")

cl <- fread("TCGA-LUSC.GDC_phenotype.tsv.gz")

# View data structure
colnames(cl)
dim(cl)
tmp <- as.data.frame(colnames(cl))

# ============================================================
# Part 3: Extract clinical variables of interest
# ============================================================
# Locate target columns
which(colnames(cl) == 'submitter_id.samples')
which(colnames(cl) == 'vital_status.demographic')
table(cl[, 'age_at_initial_pathologic_diagnosis'])

# Extract clinical metadata (adjust column indices as needed)
meta <- as.data.frame(cl[, c(1, 80, 98, 42:44, 6, 76, 86, 78, 61, 103, 104, 102, 36, 53)])

# Check extracted data
colnames(meta)
meta[1:4, 1:4]

# Save and reload metadata
write.csv(meta, 'meta.csv', row.names = TRUE)
meta <- read.csv("meta.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Save as Rdata format
save(cl, meta, file = 'GDC_TCGA_LUAD_clinical_df.Rdata')
load(file = 'GDC_TCGA_LUAD_clinical_df.Rdata')

# ============================================================
# Part 4: Data cleaning and survival time calculation
# ============================================================
phe <- meta

# Fill missing values
phe[, 8][is.na(phe[, 8])] <- 0
phe[, 9][is.na(phe[, 9])] <- 0

# Calculate survival time (days and months)
phe$days <- as.numeric(phe[, 9]) + as.numeric(phe[, 8])
phe$time <- phe$days / 30

# Age grouping based on median
phe$age_at_initial_pathologic_diagnosis <- as.numeric(phe$age_at_initial_pathologic_diagnosis)
phe$age_group <- ifelse(phe$age_at_initial_pathologic_diagnosis > median(phe$age_at_initial_pathologic_diagnosis, na.rm = TRUE), 
                        'older', 'younger')
table(phe$age_group)

# ============================================================
# Part 5: Stage and pathologic variable cleaning
# ============================================================
# Tumor stage
table(phe$tumor_stage.diagnoses)
phe <- phe[phe$tumor_stage.diagnoses != 'not reported', ]
phe$tumor_stage.diagnoses <- gsub('[stage]', '', phe$tumor_stage.diagnoses)
phe$tumor_stage.diagnoses <- gsub('[a-c]', '', phe$tumor_stage.diagnoses)
phe$tumor_stage.diagnoses <- toupper(phe$tumor_stage.diagnoses)

# T category
phe <- phe[!phe$pathologic_T == '', ]
phe$pathologic_T <- gsub('[a-b]', '', phe$pathologic_T)

# N category (remove NX)
phe <- phe[!phe$pathologic_N == 'NX', ]
phe$pathologic_N <- as.factor(phe$pathologic_N)

# M category (remove MX)
phe <- phe[!phe$pathologic_M == '', ]
phe <- phe[!phe$pathologic_M == 'MX', ]
phe$pathologic_M <- gsub('[a-b]', '', phe$pathologic_M)

# Save cleaned data
write.csv(phe, 'phe.csv', row.names = TRUE)
save(phe, file = 'GDC_TCGA_LUAD_clinical_ok.Rdata')

# ============================================================
# Part 6: Create Table 1 (Clinical characteristics)
# ============================================================
load(file = 'GDC_TCGA_LUAD_clinical_ok.Rdata')
clin <- phe

# Convert variables to ordered factors
clin$age_at_initial_pathologic_diagnosis <- as.numeric(clin$age_at_initial_pathologic_diagnosis)
clin$AGE <- factor(ifelse(clin$age_at_initial_pathologic_diagnosis > 60, '>60', '<=60'), ordered = TRUE)
clin$gender.demographic <- factor(toupper(clin$gender.demographic), ordered = TRUE)
clin$stage <- factor(clin$tumor_stage.diagnoses, ordered = TRUE)
clin$pathologic_T <- factor(clin$pathologic_T, ordered = TRUE)
clin$pathologic_N <- factor(clin$pathologic_N, ordered = TRUE)
clin$pathologic_M <- factor(clin$pathologic_M, ordered = TRUE)
clin$vital_status.demographic <- factor(toupper(clin$vital_status.demographic), ordered = TRUE)

# Define variables for table
myVars <- colnames(clin)[c(2:8)] 
catVars <- myVars[c(2:8)]

# Generate Table 1
tb_all <- CreateTableOne(vars = myVars, data = clin, factorVars = catVars)
tab_out_all <- print(tb_all, catDigits = 1, contDigits = 2, pDigits = 3,
                     quote = FALSE, missing = TRUE, explain = TRUE, 
                     printToggle = TRUE, test = TRUE, smd = TRUE,
                     showAllLevels = TRUE)
write.csv(tab_out_all, file = "TCGA-LUAD_clinical_table_all.csv")

# ============================================================
# Part 7: Prepare data for Cox regression (convert to numeric)
# ============================================================
load(file = 'GDC_TCGA_LUAD_clinical_ok.Rdata')
clin <- phe

# Convert variables to numeric for Cox regression
clin$vital_status.demographic <- ifelse(clin$vital_status.demographic == 'Alive', 0, 1)
clin$tumor_stage.diagnoses <- as.numeric(as.factor(clin$tumor_stage.diagnoses))
clin$age_group <- as.numeric(as.factor(clin$age_group))
clin$gender.demographic <- as.numeric(as.factor(clin$gender.demographic))
clin$pathologic_T <- as.numeric(as.factor(clin$pathologic_T))
clin$pathologic_N <- as.numeric(as.factor(clin$pathologic_N))
clin$pathologic_M <- as.numeric(as.factor(clin$pathologic_M))

# Save Cox-ready data
write.csv(clin, 'GDC_TCGA_LUAD_clinical_cox.csv', row.names = TRUE)
save(clin, file = 'GDC_TCGA_LUAD_clinical_cox.Rdata')

# ============================================================
# Part 8: Univariate Cox regression analysis
# ============================================================
load(file = 'GDC_TCGA_LUAD_clinical_cox.Rdata')
y <- Surv(clin$time, clin$vital_status.demographic)

# Univariate Cox function
Unicox_model <- function(x) {
  FML <- as.formula(paste0("y ~ ", x))
  cox <- coxph(FML, data = clin)
  cox_sum <- summary(cox)
  HR <- round(cox_sum$coefficients[, 2], 2)
  Pvalue <- round(cox_sum$coefficients[, 5], 3)
  CI5 <- round(cox_sum$conf.int[, 3], 2)
  CI95 <- round(cox_sum$conf.int[, 4], 2)
  CI <- paste0(HR, " (", CI5, "-", CI95, ")")
  Variable <- rownames(cox_sum$coefficients)
  return(data.frame(Variable = Variable, HR_CI = CI, Pvalue = Pvalue))
}

# Define variables for univariate analysis
variables <- c("tumor_stage.diagnoses", "pathologic_M", "pathologic_N", 
               "pathologic_T", "age_at_initial_pathologic_diagnosis", 
               "gender.demographic", "tobacco_smoking_history")

# Run univariate analysis
Unicox_list <- lapply(variables, Unicox_model)
Unicox <- ldply(Unicox_list, data.frame)
colnames(Unicox) <- c("Variable", "HR (95% CI)", "Pvalue")
View(Unicox)

# Save univariate results
write.csv(Unicox, 'Univariate_Cox_results.csv', row.names = FALSE)

# ============================================================
# Part 9: Multivariate Cox regression analysis
# ============================================================
# Run multivariate Cox model
res <- coxph(Surv(time, vital_status.demographic) ~ tumor_stage.diagnoses + 
               pathologic_M + pathologic_N + pathologic_T + 
               age_at_initial_pathologic_diagnosis + gender.demographic, 
             data = clin)
mul_cox <- summary(res)

# Extract results
mul_HR <- round(mul_cox$coefficients[, 2], 2)
mul_Pvalue <- round(mul_cox$coefficients[, 5], 4)
mul_CI5 <- round(mul_cox$conf.int[, 3], 2)
mul_CI95 <- round(mul_cox$conf.int[, 4], 2)
mul_CI <- paste0(mul_HR, ' (', mul_CI5, '-', mul_CI95, ')')
Variable <- rownames(mul_cox$coefficients)

# Create results data frame
mulcox_res <- data.frame(Variable, mul_HR, mul_CI5, mul_CI95, mul_CI, mul_Pvalue)
colnames(mulcox_res) <- c("Variable", "HR", "CI5", "CI95", "HR (95% CI)", "Pvalue")
View(mulcox_res)

# Save multivariate results
write.csv(mulcox_res, 'Multivariate_Cox_results.csv', row.names = FALSE)

# ============================================================
# Part 10: Forest plot for multivariate results
# ============================================================
# Prepare data for forest plot with header
dat_forest <- rbind(c("Variable", NA, NA, NA, "HR (95% CI)", "Pvalue"), 
                    mulcox_res)

# 10.1 TIFF version (300 dpi for submission)
tiff("Forest_plot_multivariate_Cox.tiff", 
     width = 8, height = 6, units = "in", 
     res = 300, compression = "lzw")

forestplot(dat_forest[, c(1, 5, 6)],   
           mean = as.numeric(dat_forest[, 2]), 
           lower = as.numeric(dat_forest[, 3]), 
           upper = as.numeric(dat_forest[, 4]), 
           zero = 1,                
           boxsize = 0.2,           
           graph.pos = 3,         
           xticks = c(0, 1, 2, 3), 
           txt_gp = fpTxtGp(label = gpar(cex = 0.8), ticks = gpar(cex = 0.6)),
           hrzl_lines = list("1" = gpar(lty = 1, lwd = 1.5),
                             "2" = gpar(lty = 1, lwd = 1.5),
                             "7" = gpar(lty = 1, lwd = 1.5)),
           col = fpColors(box = "blue", lines = "black", zero = "grey"),
           lwd.zero = 1,
           lwd.ci = 1.5,
           lty.ci = 2,             
           ci.vertices.height = 0.1)

dev.off()
cat("TIFF saved: Forest_plot_multivariate_Cox.tiff (300 dpi)\n")

# 10.2 PNG version (300 dpi for preview)
png("Forest_plot_multivariate_Cox.png", 
    width = 8, height = 6, units = "in", res = 300)

forestplot(dat_forest[, c(1, 5, 6)],
           mean = as.numeric(dat_forest[, 2]),
           lower = as.numeric(dat_forest[, 3]),
           upper = as.numeric(dat_forest[, 4]),
           zero = 1, boxsize = 0.2, graph.pos = 3,
           xticks = c(0, 1, 2, 3),
           txt_gp = fpTxtGp(label = gpar(cex = 0.8), ticks = gpar(cex = 0.6)),
           hrzl_lines = list("1" = gpar(lty = 1, lwd = 1.5),
                             "2" = gpar(lty = 1, lwd = 1.5),
                             "7" = gpar(lty = 1, lwd = 1.5)),
           col = fpColors(box = "blue", lines = "black", zero = "grey"),
           lwd.zero = 1, lwd.ci = 1.5, lty.ci = 2, ci.vertices.height = 0.1)
dev.off()
cat("PNG saved: Forest_plot_multivariate_Cox.png (300 dpi)\n")

# ============================================================
# Part 11: Summary report
# ============================================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("TCGA CLINICAL DATA PROCESSING AND COX ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Generated files:\n")
cat("1. meta.csv                                    - Extracted metadata\n")
cat("2. GDC_TCGA_LUAD_clinical_df.Rdata             - Raw clinical data\n")
cat("3. phe.csv                                     - Cleaned phenotype data\n")
cat("4. GDC_TCGA_LUAD_clinical_ok.Rdata             - Cleaned clinical data\n")
cat("5. TCGA-LUAD_clinical_table_all.csv            - Table 1\n")
cat("6. GDC_TCGA_LUAD_clinical_cox.csv              - Cox-ready data\n")
cat("7. GDC_TCGA_LUAD_clinical_cox.Rdata            - Cox-ready Rdata\n")
cat("8. Univariate_Cox_results.csv                  - Univariate results\n")
cat("9. Multivariate_Cox_results.csv                - Multivariate results\n")
cat("10. Forest_plot_multivariate_Cox.tiff          - Forest plot (300 dpi)\n")
cat("11. Forest_plot_multivariate_Cox.png           - Forest plot preview\n\n")

cat("Analysis completed at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 60), "\n", sep = "")