# ============================================================================
# Clinical_Cox_univariate_multivariate.R
# 
# Description: Univariate and multivariate Cox regression analysis for 
#              clinical features (Stage, T category, N category, M category,
#              Age, Gender, Smoking history)
# 
# Special feature: TX/NX/MX categories are preserved as independent groups
# 
# Output: Forest plots (TIFF 500 dpi, PNG 300 dpi, PDF vector)
#         Results tables (CSV format)
# ============================================================================

# Load required packages
library(survival)
library(plyr)
library(forestplot)
library(grid)
library(dplyr)

# Set working directory (modify to your path)
setwd("C:/Users/Administrator/Desktop/clinic")

# ============================================================
# Part 1: Data preprocessing (Preserve TX/NX/MX)
# ============================================================
cat("Starting data preprocessing...\n")

# Read clinical data
clin <- read.csv("clinic.csv", header = TRUE, row.names = 1, check.names = FALSE)

# View data structure
cat("\nColumn names:\n")
print(names(clin))
cat("\nData dimensions:", dim(clin), "\n")

# ============================================================
# Variable transformation (Preserve TX/NX/MX)
# ============================================================

# Survival status: Alive -> 0, Dead -> 1
clin$status <- ifelse(clin$vital_status == 'Alive', 0, 1)
cat("\nSurvival status conversion complete\n")
print(table(clin$vital_status, clin$status, dnn = c("Original", "Converted")))

# Gender: female -> 0, male -> 1
clin$gender_num <- ifelse(clin$gender == 'female', 0, 1)
cat("\nGender conversion complete\n")
print(table(clin$gender, clin$gender_num, dnn = c("Original", "Converted")))

# AJCC stage conversion
cat("\nAJCC stage distribution:\n")
print(table(clin$ajcc_stage))
clin$stage_num <- as.numeric(factor(clin$ajcc_stage, 
                                    levels = c("Stage IA", "Stage IB", "Stage IIA", "Stage IIB", 
                                               "Stage IIIA", "Stage IIIB", "Stage IV")))
cat("AJCC stage conversion complete\n")

# ==================== T category (Preserve TX) ====================
cat("\n", paste(rep("-", 30), collapse = ""), "\n")
cat("T category processing (TX preserved as independent category)\n")
cat(paste(rep("-", 30), collapse = ""), "\n")
cat("Original T category distribution:\n")
print(table(clin$T_category, useNA = "ifany"))

# T category numeric coding: T1,T2,T3,T4,TX -> 1,2,3,4,5
clin$T_num <- as.numeric(factor(clin$T_category, 
                                levels = c("T1", "T2", "T3", "T4", "TX")))
cat("\nConverted T category distribution:\n")
print(table(clin$T_num, clin$T_category, dnn = c("Numeric", "Original")))
cat("T category conversion complete\n")

# ==================== N category (Preserve NX) ====================
cat("\n", paste(rep("-", 30), collapse = ""), "\n")
cat("N category processing (NX preserved as independent category)\n")
cat(paste(rep("-", 30), collapse = ""), "\n")
cat("Original N category distribution:\n")
print(table(clin$N_category, useNA = "ifany"))

# N category numeric coding: N0,N1,N2,N3,NX -> 0,1,2,3,4
clin$N_num <- as.numeric(factor(clin$N_category, 
                                levels = c("N0", "N1", "N2", "N3", "NX"))) - 1
cat("\nConverted N category distribution:\n")
print(table(clin$N_num, clin$N_category, dnn = c("Numeric", "Original")))
cat("N category conversion complete\n")

# ==================== M category (Preserve MX) ====================
cat("\n", paste(rep("-", 30), collapse = ""), "\n")
cat("M category processing (MX preserved as independent category)\n")
cat(paste(rep("-", 30), collapse = ""), "\n")
cat("Original M category distribution:\n")
print(table(clin$M_category, useNA = "ifany"))

# M category numeric coding: M0 -> 0, M1 -> 1, MX -> 2
clin$M_num <- case_when(
  clin$M_category == "M0" ~ 0,
  clin$M_category == "M1" ~ 1,
  clin$M_category == "MX" ~ 2,
  TRUE ~ NA_real_
)
cat("\nConverted M category distribution:\n")
print(table(clin$M_num, clin$M_category, dnn = c("Numeric", "Original")))
cat("M category conversion complete\n")

# Smoking history
clin$smoking_num <- clin$smoking
cat("\nSmoking history distribution:\n")
print(table(clin$smoking))

# Age
clin$age_num <- clin$age

# Survival time
clin$time <- clin$survival_days

# ============================================================
# Create analysis dataset
# ============================================================
analysis_data <- data.frame(
  time = clin$time,
  status = clin$status,
  Age = clin$age_num,
  Sex = clin$gender_num,
  Stage = clin$stage_num,
  T = clin$T_num,
  N = clin$N_num,
  M = clin$M_num,
  Smoking = clin$smoking_num
)

# Check for missing values
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Missing value check (should all be 0)\n")
cat(paste(rep("=", 50), collapse = ""), "\n")
missing_counts <- colSums(is.na(analysis_data))
print(missing_counts)

if(sum(missing_counts) == 0) {
  cat("No missing values detected\n")
} else {
  cat("Warning: Missing values detected\n")
}

# Save processed data
write.csv(analysis_data, 'clinical_data_processed.csv', row.names = TRUE)
save(analysis_data, file = 'clinical_data_processed.Rdata')
cat("\nData preprocessing complete\n")

# ============================================================
# Part 2: Univariate Cox regression analysis
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Starting Univariate Cox regression analysis...\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Create survival object
y <- Surv(analysis_data$time, analysis_data$status)

# Univariate Cox analysis function
Unicox_model <- function(x){
  FML <- as.formula(paste0("y ~ ", x))
  cox <- coxph(FML, data = analysis_data)
  cox1 <- summary(cox)
  
  HR <- round(cox1$coefficients[,2], 2)
  Pvalue <- round(cox1$coefficients[,5], 4)
  CI5 <- round(cox1$conf.int[,3], 2)
  CI95 <- round(cox1$conf.int[,4], 2)
  CI <- paste0(HR, " (", CI5, "-", CI95, ")")
  
  Variable <- x
  Unicox_model <- data.frame(Variable, HR, CI5, CI95, CI, Pvalue)
  return(Unicox_model)
} 

# Define analysis variables
variable <- c("Age", "Sex", "Stage", "T", "N", "M", "Smoking")

# Run univariate analysis
Unicox <- lapply(variable, Unicox_model)
Unicox <- ldply(Unicox, data.frame)
colnames(Unicox) <- c("Variable", "HR", "CI5", "CI95", "HR (95% CI)", "Pvalue")

# Save univariate results
write.csv(Unicox, 'Univariate_Cox_results.csv', row.names = FALSE)
cat("Univariate Cox analysis complete. Results saved to Univariate_Cox_results.csv\n")

# Display univariate results
cat("\nUnivariate Cox regression results:\n")
print(Unicox)

# ============================================================
# Part 3: Multivariate Cox regression analysis
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Starting Multivariate Cox regression analysis...\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Run multivariate Cox analysis
res <- coxph(Surv(time, status) ~ Age + Stage + T + N, data = analysis_data)
mul_cox <- summary(res)

# Extract multivariate results
mul_HR <- round(mul_cox$coefficients[,2], 2) 
mul_Pvalue <- round(mul_cox$coefficients[,5], 4)
mul_CI5 <- round(mul_cox$conf.int[,3], 2)
mul_CI95 <- round(mul_cox$conf.int[,4], 2)
mul_CI <- paste0(mul_HR, ' (', mul_CI5, '-', mul_CI95, ')')

# Get variable names
Variable <- row.names(data.frame(mul_cox$coefficients))

mulcox_res <- data.frame(Variable = Variable, 
                         HR = mul_HR, 
                         CI5 = mul_CI5, 
                         CI95 = mul_CI95, 
                         CI = mul_CI, 
                         Pvalue = mul_Pvalue)

colnames(mulcox_res) <- c("Variable", "HR", "CI5", "CI95", "HR (95% CI)", "Pvalue")

# Save multivariate results
write.csv(mulcox_res, 'Multivariate_Cox_results.csv', row.names = FALSE)
cat("Multivariate Cox analysis complete. Results saved to Multivariate_Cox_results.csv\n")

# Display multivariate results
cat("\nMultivariate Cox regression results:\n")
print(mulcox_res)

# ============================================================
# Part 4: Figure settings (SCI standard)
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Setting figure parameters (SCI standard)...\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Fixed dimensions for SCI submission
FIG_WIDTH <- 7.5    # 7.5 inches
FIG_HEIGHT <- 6     # 6 inches

# Unified color scheme
UNIFIED_COLORS <- list(
  box = "#1c61b6",
  lines = "#1c61b6",
  zero = "#666666"
)

# Unified graphic parameters
UNIFIED_SETTINGS <- list(
  boxsize = 0.18,
  graph.pos = 2,
  zero = 1,
  xlog = TRUE,
  lwd.zero = 1.5,
  lwd.ci = 1.5,
  ci.vertices = TRUE,
  ci.vertices.height = 0.12,
  xticks = c(0.5, 1, 2, 4, 8),
  clip = c(0.1, 10)
)

# ============================================================
# Part 5: Univariate forest plot
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Generating Univariate forest plot...\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Pretty variable names for publication
pretty_names <- c(
  "Age" = "Age",
  "Sex" = "Gender",
  "Stage" = "AJCC stage",
  "T" = "T category (T1,T2,T3,T4,TX)",
  "N" = "N category (N0,N1,N2,N3,NX)",
  "M" = "M category (M0,M1,MX)",
  "Smoking" = "Smoking history"
)

# Prepare table text
uni_tabletext <- cbind(
  c("Variable", pretty_names[as.character(Unicox$Variable)]),
  c("HR (95% CI)", as.character(Unicox$`HR (95% CI)`)),
  c("P Value", ifelse(Unicox$Pvalue < 0.001, "<0.001", 
                      format(round(Unicox$Pvalue, 4), nsmall = 4)))
)

# Font settings
txt_gp_uni <- fpTxtGp(
  label = gpar(cex = 0.85, fontfamily = "Helvetica"),
  ticks = gpar(cex = 0.75, fontfamily = "Helvetica"),
  xlab = gpar(cex = 0.9, fontfamily = "Helvetica", fontface = "bold"),
  title = gpar(cex = 1.0, fontfamily = "Helvetica", fontface = "bold")
)

# Horizontal line settings
uni_hrzl_lines <- list(
  "1" = gpar(lty = 1, lwd = 1.5, col = "black"),
  "2" = gpar(lty = 1, lwd = 1, col = "black")
)
uni_last_line <- as.character(nrow(Unicox) + 2)
uni_hrzl_lines[[uni_last_line]] <- gpar(lty = 1, lwd = 1.5, col = "black")

# 5.1 TIFF version (500 dpi for submission)
cat("  Generating TIFF version (500 dpi for submission)...\n")
tiff("Figure_S1_Univariate_Cox.tiff", 
     width = FIG_WIDTH, 
     height = FIG_HEIGHT, 
     units = "in", 
     res = 500,
     compression = "lzw")

forestplot(
  labeltext = uni_tabletext,
  mean = c(NA, Unicox$HR),
  lower = c(NA, Unicox$CI5),
  upper = c(NA, Unicox$CI95),
  is.summary = c(TRUE, rep(FALSE, nrow(Unicox))),
  zero = UNIFIED_SETTINGS$zero,
  xlog = UNIFIED_SETTINGS$xlog,
  boxsize = UNIFIED_SETTINGS$boxsize,
  graph.pos = UNIFIED_SETTINGS$graph.pos,
  xticks = UNIFIED_SETTINGS$xticks,
  clip = UNIFIED_SETTINGS$clip,
  xlab = "Hazard Ratio (HR)",
  txt_gp = txt_gp_uni,
  hrzl_lines = uni_hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(
    box = UNIFIED_COLORS$box,
    lines = UNIFIED_COLORS$lines,
    zero = UNIFIED_COLORS$zero
  ),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Univariate Cox Regression Analysis"
)
dev.off()
cat("  OK TIFF saved: Figure_S1_Univariate_Cox.tiff (500 dpi)\n")

# 5.2 PNG version (300 dpi for preview)
cat("  Generating PNG version (300 dpi for preview)...\n")
png("Figure_S1_Univariate_Cox.png", 
    width = FIG_WIDTH * 300, 
    height = FIG_HEIGHT * 300, 
    res = 300,
    bg = "white")

forestplot(
  labeltext = uni_tabletext,
  mean = c(NA, Unicox$HR),
  lower = c(NA, Unicox$CI5),
  upper = c(NA, Unicox$CI95),
  is.summary = c(TRUE, rep(FALSE, nrow(Unicox))),
  zero = UNIFIED_SETTINGS$zero,
  xlog = UNIFIED_SETTINGS$xlog,
  boxsize = UNIFIED_SETTINGS$boxsize,
  graph.pos = UNIFIED_SETTINGS$graph.pos,
  xticks = UNIFIED_SETTINGS$xticks,
  clip = UNIFIED_SETTINGS$clip,
  xlab = "Hazard Ratio (HR)",
  txt_gp = txt_gp_uni,
  hrzl_lines = uni_hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(
    box = UNIFIED_COLORS$box,
    lines = UNIFIED_COLORS$lines,
    zero = UNIFIED_COLORS$zero
  ),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Univariate Cox Regression Analysis"
)
dev.off()
cat("  OK PNG saved: Figure_S1_Univariate_Cox.png (300 dpi)\n")

# 5.3 PDF version (vector graphic)
cat("  Generating PDF version (vector graphic)...\n")
pdf("Figure_S1_Univariate_Cox.pdf", width = FIG_WIDTH, height = FIG_HEIGHT)

forestplot(
  labeltext = uni_tabletext,
  mean = c(NA, Unicox$HR),
  lower = c(NA, Unicox$CI5),
  upper = c(NA, Unicox$CI95),
  is.summary = c(TRUE, rep(FALSE, nrow(Unicox))),
  zero = UNIFIED_SETTINGS$zero,
  xlog = UNIFIED_SETTINGS$xlog,
  boxsize = UNIFIED_SETTINGS$boxsize,
  graph.pos = UNIFIED_SETTINGS$graph.pos,
  xticks = UNIFIED_SETTINGS$xticks,
  clip = UNIFIED_SETTINGS$clip,
  xlab = "Hazard Ratio (HR)",
  txt_gp = txt_gp_uni,
  hrzl_lines = uni_hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(
    box = UNIFIED_COLORS$box,
    lines = UNIFIED_COLORS$lines,
    zero = UNIFIED_COLORS$zero
  ),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Univariate Cox Regression Analysis"
)
dev.off()
cat("  OK PDF saved: Figure_S1_Univariate_Cox.pdf (vector)\n")

# ============================================================
# Part 6: Multivariate forest plot
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Generating Multivariate forest plot...\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Pretty variable names for multivariate
multi_pretty_names <- c(
  "Age" = "Age",
  "Stage" = "AJCC stage",
  "T" = "T category",
  "N" = "N category"
)

# Prepare table text
multi_tabletext <- cbind(
  c("Variable", multi_pretty_names[as.character(mulcox_res$Variable)]),
  c("HR (95% CI)", as.character(mulcox_res$`HR (95% CI)`)),
  c("P Value", ifelse(mulcox_res$Pvalue < 0.001, "<0.001", 
                      format(round(mulcox_res$Pvalue, 4), nsmall = 4)))
)

# Font settings
txt_gp_multi <- fpTxtGp(
  label = gpar(cex = 0.85, fontfamily = "Helvetica"),
  ticks = gpar(cex = 0.75, fontfamily = "Helvetica"),
  xlab = gpar(cex = 0.9, fontfamily = "Helvetica", fontface = "bold"),
  title = gpar(cex = 1.0, fontfamily = "Helvetica", fontface = "bold")
)

# Horizontal line settings
multi_hrzl_lines <- list(
  "1" = gpar(lty = 1, lwd = 1.5, col = "black"),
  "2" = gpar(lty = 1, lwd = 1, col = "black")
)
multi_last_line <- as.character(nrow(mulcox_res) + 2)
multi_hrzl_lines[[multi_last_line]] <- gpar(lty = 1, lwd = 1.5, col = "black")

# 6.1 TIFF version (500 dpi for submission)
cat("  Generating TIFF version (500 dpi for submission)...\n")
tiff("Figure_1_Multivariate_Cox.tiff", 
     width = FIG_WIDTH, 
     height = FIG_HEIGHT, 
     units = "in", 
     res = 500,
     compression = "lzw")

forestplot(
  labeltext = multi_tabletext,
  mean = c(NA, mulcox_res$HR),
  lower = c(NA, mulcox_res$CI5),
  upper = c(NA, mulcox_res$CI95),
  is.summary = c(TRUE, rep(FALSE, nrow(mulcox_res))),
  zero = UNIFIED_SETTINGS$zero,
  xlog = UNIFIED_SETTINGS$xlog,
  boxsize = UNIFIED_SETTINGS$boxsize,
  graph.pos = UNIFIED_SETTINGS$graph.pos,
  xticks = UNIFIED_SETTINGS$xticks,
  clip = UNIFIED_SETTINGS$clip,
  xlab = "Hazard Ratio (HR)",
  txt_gp = txt_gp_multi,
  hrzl_lines = multi_hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(
    box = UNIFIED_COLORS$box,
    lines = UNIFIED_COLORS$lines,
    zero = UNIFIED_COLORS$zero
  ),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Multivariate Cox Regression Analysis"
)
dev.off()
cat("  OK TIFF saved: Figure_1_Multivariate_Cox.tiff (500 dpi)\n")

# 6.2 PNG version (300 dpi for preview)
cat("  Generating PNG version (300 dpi for preview)...\n")
png("Figure_1_Multivariate_Cox.png", 
    width = FIG_WIDTH * 300, 
    height = FIG_HEIGHT * 300, 
    res = 300,
    bg = "white")

forestplot(
  labeltext = multi_tabletext,
  mean = c(NA, mulcox_res$HR),
  lower = c(NA, mulcox_res$CI5),
  upper = c(NA, mulcox_res$CI95),
  is.summary = c(TRUE, rep(FALSE, nrow(mulcox_res))),
  zero = UNIFIED_SETTINGS$zero,
  xlog = UNIFIED_SETTINGS$xlog,
  boxsize = UNIFIED_SETTINGS$boxsize,
  graph.pos = UNIFIED_SETTINGS$graph.pos,
  xticks = UNIFIED_SETTINGS$xticks,
  clip = UNIFIED_SETTINGS$clip,
  xlab = "Hazard Ratio (HR)",
  txt_gp = txt_gp_multi,
  hrzl_lines = multi_hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(
    box = UNIFIED_COLORS$box,
    lines = UNIFIED_COLORS$lines,
    zero = UNIFIED_COLORS$zero
  ),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Multivariate Cox Regression Analysis"
)
dev.off()
cat("  OK PNG saved: Figure_1_Multivariate_Cox.png (300 dpi)\n")

# 6.3 PDF version (vector graphic)
cat("  Generating PDF version (vector graphic)...\n")
pdf("Figure_1_Multivariate_Cox.pdf", width = FIG_WIDTH, height = FIG_HEIGHT)

forestplot(
  labeltext = multi_tabletext,
  mean = c(NA, mulcox_res$HR),
  lower = c(NA, mulcox_res$CI5),
  upper = c(NA, mulcox_res$CI95),
  is.summary = c(TRUE, rep(FALSE, nrow(mulcox_res))),
  zero = UNIFIED_SETTINGS$zero,
  xlog = UNIFIED_SETTINGS$xlog,
  boxsize = UNIFIED_SETTINGS$boxsize,
  graph.pos = UNIFIED_SETTINGS$graph.pos,
  xticks = UNIFIED_SETTINGS$xticks,
  clip = UNIFIED_SETTINGS$clip,
  xlab = "Hazard Ratio (HR)",
  txt_gp = txt_gp_multi,
  hrzl_lines = multi_hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(
    box = UNIFIED_COLORS$box,
    lines = UNIFIED_COLORS$lines,
    zero = UNIFIED_COLORS$zero
  ),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Multivariate Cox Regression Analysis"
)
dev.off()
cat("  OK PDF saved: Figure_1_Multivariate_Cox.pdf (vector)\n")

# ============================================================
# Part 7: Create SCI format tables
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("Creating SCI format tables...\n")
cat(paste(rep("=", 50), collapse = ""), "\n")

# Table S1: Univariate Cox analysis
table_uni <- data.frame(
  Variable = pretty_names[as.character(Unicox$Variable)],
  `HR (95% CI)` = Unicox$`HR (95% CI)`,
  `P Value` = ifelse(Unicox$Pvalue < 0.001, "<0.001", 
                     format(round(Unicox$Pvalue, 4), nsmall = 4))
)
write.csv(table_uni, 'Supplementary_Table_S1_Univariate_Cox.csv', row.names = FALSE)
cat("Supplementary Table S1 saved: Supplementary_Table_S1_Univariate_Cox.csv\n")

# Table 1: Multivariate Cox analysis
table_multi <- data.frame(
  Variable = multi_pretty_names[as.character(mulcox_res$Variable)],
  `HR (95% CI)` = mulcox_res$`HR (95% CI)`,
  `P Value` = ifelse(mulcox_res$Pvalue < 0.001, "<0.001", 
                     format(round(mulcox_res$Pvalue, 4), nsmall = 4))
)
write.csv(table_multi, 'Table_1_Multivariate_Cox.csv', row.names = FALSE)
cat("Table 1 saved: Table_1_Multivariate_Cox.csv\n")

# ============================================================
# Part 8: Summary report
# ============================================================
cat("\n", paste(rep("=", 50), collapse = ""), "\n")
cat("COX REGRESSION ANALYSIS COMPLETE\n")
cat(paste(rep("=", 50), collapse = ""), "\n\n")

cat("Data preprocessing complete (TX/NX/MX preserved as independent categories)\n")
cat("Univariate Cox analysis complete (7 variables)\n")
cat("Multivariate Cox analysis complete (4 variables)\n")
cat("SCI format figures and tables generated\n\n")

cat("Generated files:\n")
cat("1. clinical_data_processed.csv              - Processed clinical data\n")
cat("2. clinical_data_processed.Rdata            - R format data\n")
cat("3. Univariate_Cox_results.csv               - Univariate Cox results\n")
cat("4. Multivariate_Cox_results.csv             - Multivariate Cox results\n")
cat("5. Figure_S1_Univariate_Cox.tiff            - Univariate forest plot (500 dpi)\n")
cat("6. Figure_S1_Univariate_Cox.png             - Univariate forest plot (preview)\n")
cat("7. Figure_S1_Univariate_Cox.pdf             - Univariate forest plot (vector)\n")
cat("8. Figure_1_Multivariate_Cox.tiff           - Multivariate forest plot (500 dpi)\n")
cat("9. Figure_1_Multivariate_Cox.png            - Multivariate forest plot (preview)\n")
cat("10. Figure_1_Multivariate_Cox.pdf           - Multivariate forest plot (vector)\n")
cat("11. Supplementary_Table_S1_Univariate_Cox.csv - Table S1\n")
cat("12. Table_1_Multivariate_Cox.csv            - Table 1\n\n")

cat("Figure parameters:\n")
cat("TIFF resolution: 500 dpi (Elsevier compliant)\n")
cat("PNG resolution: 300 dpi (preview)\n")
cat("PDF format: Vector (scale-free)\n")
cat("Compression: LZW lossless for TIFF\n\n")

cat("Variable names (Lung Cancer journal compliant):\n")
cat("  - AJCC stage\n")
cat("  - T category (T1,T2,T3,T4,TX) - TX preserved\n")
cat("  - N category (N0,N1,N2,N3,NX) - NX preserved\n")
cat("  - M category (M0,M1,MX) - MX preserved\n")
cat("  - Age\n")
cat("  - Gender\n")
cat("  - Smoking history\n\n")

cat("Analysis completed at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 50), collapse = ""), "\n")