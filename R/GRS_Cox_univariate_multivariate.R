# ============================================================================
# GRS_Cox_univariate_multivariate.R
# 
# Description: Univariate and multivariate Cox regression analysis with GRS
#              (Genetic Risk Score)
# 
# Functions:
#   1. Univariate Cox analysis for each clinical variable + GRS
#   2. Multivariate Cox analysis with selected variables
#   3. Forest plot for visualization (SCI quality, 500 DPI)
# 
# Output: 
#   - Forest plots (TIFF/PNG/PDF)
#   - Results tables (CSV)
# ============================================================================

# ============================================================
# 1. Load packages
# ============================================================
library(survival)
library(plyr)
library(forestplot)
library(grid)
library(dplyr)

# ============================================================
# 2. Set working directory (modify to your path)
# ============================================================
setwd("C:/Users/Administrator/Desktop/GRS")

# ============================================================
# 3. Data preprocessing
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Starting data preprocessing...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Read data
data <- read.csv("clinic-GRS.csv", header = TRUE, row.names = 1, check.names = FALSE)

# View data structure
cat("Data columns:\n")
print(names(data))
cat("\nData dimensions:", dim(data), "\n")
cat("\nFirst few rows:\n")
print(head(data))

# ============================================================
# 4. Variable transformation
# ============================================================

# Survival status: Alive -> 0, Dead -> 1
data$status <- ifelse(data$status == 'Alive', 0, 1)
cat("\nSurvival status conversion complete\n")
print(table(data$status, dnn = "Status (0=Alive, 1=Dead)"))

# Gender: female -> 0, male -> 1
data$gender_num <- ifelse(data$sex == 'female', 0, 1)
cat("\nGender conversion complete\n")
print(table(data$sex, data$gender_num, dnn = c("Original", "Converted")))

# AJCC stage conversion
data$stage_num <- as.numeric(factor(data$stage, 
                                    levels = c("Stage IA", "Stage IB", "Stage IIA", "Stage IIB", 
                                               "Stage IIIA", "Stage IIIB", "Stage IV")))
cat("\nAJCC stage conversion complete\n")
print(table(data$stage, data$stage_num, dnn = c("Original stage", "Converted")))

# T stage conversion
data$T_num <- as.numeric(factor(data$pathologic_T, 
                                levels = c("T1", "T2", "T3", "T4", "TX")))
cat("\nT stage conversion complete\n")
print(table(data$pathologic_T, data$T_num, dnn = c("Original T", "Converted")))

# N stage conversion
data$N_num <- as.numeric(factor(data$pathologic_N, 
                                levels = c("N0", "N1", "N2", "N3", "NX"))) - 1
cat("\nN stage conversion complete\n")
print(table(data$pathologic_N, data$N_num, dnn = c("Original N", "Converted")))

# M stage conversion
data$M_num <- case_when(
  data$pathologic_M == "M0" ~ 0,
  data$pathologic_M == "M1" ~ 1,
  data$pathologic_M == "MX" ~ 2,
  TRUE ~ NA_real_
)
cat("\nM stage conversion complete\n")
print(table(data$pathologic_M, data$M_num, dnn = c("Original M", "Converted")))

# Smoking history
data$smoking_num <- data$smoking_history

# Age
data$age_num <- data$age

# Survival time
data$time <- data$time

# GRS risk score (Genetic Risk Score)
data$GRS <- data$riskscore
cat("\nGRS (Genetic Risk Score) extracted\n")
summary(data$GRS)

# ============================================================
# 5. Create analysis dataset
# ============================================================
analysis_data <- data.frame(
  time = data$time,
  status = data$status,
  Age = data$age_num,
  Sex = data$gender_num,
  Stage = data$stage_num,
  T = data$T_num,
  N = data$N_num,
  M = data$M_num,
  Smoking = data$smoking_num,
  GRS = data$GRS
)

# Check missing values
cat("\nMissing value summary:\n")
missing_counts <- colSums(is.na(analysis_data))
print(missing_counts)

# Save processed data
write.csv(analysis_data, 'GDC_TCGA_LUAD_clinical_GRS.csv', row.names = TRUE)
save(analysis_data, file = 'GDC_TCGA_LUAD_clinical_GRS.Rdata')
cat("\nData preprocessing complete and saved\n")

# ============================================================
# 6. Univariate Cox regression analysis (including GRS)
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Starting Univariate Cox regression analysis (including GRS)...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

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

# Define analysis variables (including GRS)
variable <- c("Age", "Sex", "Stage", "T", "N", "M", "Smoking", "GRS")

# Run univariate analysis
Unicox <- lapply(variable, Unicox_model)
Unicox <- ldply(Unicox, data.frame)
colnames(Unicox) <- c("Variable", "HR", "CI5", "CI95", "HR (95% CI)", "Pvalue")

# Save univariate results
write.csv(Unicox, 'Unicox.csv', row.names = FALSE)
cat("Univariate Cox analysis complete. Results saved to Unicox.csv\n")

# Display univariate results
cat("\nUnivariate Cox regression results (including GRS):\n")
print(Unicox)

# ============================================================
# 7. Multivariate Cox regression analysis (including GRS)
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Starting Multivariate Cox regression analysis (including GRS)...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Run multivariate Cox analysis (including GRS)
res <- coxph(Surv(time, status) ~ Age + Stage + T + N + GRS, data = analysis_data)
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
write.csv(mulcox_res, 'mulcox_res.csv', row.names = FALSE)
cat("Multivariate Cox analysis complete. Results saved to mulcox_res.csv\n")

# Display multivariate results
cat("\nMultivariate Cox regression results (including GRS):\n")
print(mulcox_res)

# ============================================================
# 8. Figure settings (SCI standard)
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Setting figure parameters (SCI standard)...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Fixed dimensions for SCI submission
FIG_WIDTH <- 7.5    # 7.5 inches
FIG_HEIGHT <- 6     # 6 inches

# Unified color scheme
UNIFIED_COLORS <- list(
  box = "#1c61b6",      # Professional blue
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
# 9. Univariate forest plot (including GRS)
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Generating Univariate forest plot (including GRS)...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Pretty variable names for publication
pretty_names <- c(
  "Age" = "Age",
  "Sex" = "Gender",
  "Stage" = "AJCC stage",
  "T" = "T category",
  "N" = "N category",
  "M" = "M category",
  "Smoking" = "Smoking history",
  "GRS" = "Genetic risk score"
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

# 9.1 TIFF version (500 dpi for submission)
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

# 9.2 PNG version (300 dpi for preview)
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

# 9.3 PDF version (vector graphic)
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
# 10. Multivariate forest plot (including GRS)
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Generating Multivariate forest plot (including GRS)...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Pretty variable names for multivariate
multi_pretty_names <- c(
  "Age" = "Age",
  "Stage" = "AJCC stage",
  "T" = "T category",
  "N" = "N category",
  "GRS" = "Genetic risk score"
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

# 10.1 TIFF version (500 dpi for submission)
cat("  Generating TIFF version (500 dpi for submission)...\n")
tiff("Figure_2_Multivariate_Cox.tiff", 
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
cat("  OK TIFF saved: Figure_2_Multivariate_Cox.tiff (500 dpi)\n")

# 10.2 PNG version (300 dpi for preview)
cat("  Generating PNG version (300 dpi for preview)...\n")
png("Figure_2_Multivariate_Cox.png", 
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
cat("  OK PNG saved: Figure_2_Multivariate_Cox.png (300 dpi)\n")

# 10.3 PDF version (vector graphic)
cat("  Generating PDF version (vector graphic)...\n")
pdf("Figure_2_Multivariate_Cox.pdf", width = FIG_WIDTH, height = FIG_HEIGHT)

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
cat("  OK PDF saved: Figure_2_Multivariate_Cox.pdf (vector)\n")

# ============================================================
# 11. Create SCI format tables
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Creating SCI format tables...\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

# Table 1: Univariate Cox analysis
table_uni <- data.frame(
  Variable = pretty_names[as.character(Unicox$Variable)],
  `HR (95% CI)` = Unicox$`HR (95% CI)`,
  `P Value` = ifelse(Unicox$Pvalue < 0.001, "<0.001", 
                     format(round(Unicox$Pvalue, 4), nsmall = 4))
)
write.csv(table_uni, 'Supplementary_Table1_Univariate_Cox.csv', row.names = FALSE)
cat("Supplementary Table 1 saved: Supplementary_Table1_Univariate_Cox.csv\n")

# Table 2: Multivariate Cox analysis
table_multi <- data.frame(
  Variable = multi_pretty_names[as.character(mulcox_res$Variable)],
  `HR (95% CI)` = mulcox_res$`HR (95% CI)`,
  `P Value` = ifelse(mulcox_res$Pvalue < 0.001, "<0.001", 
                     format(round(mulcox_res$Pvalue, 4), nsmall = 4))
)
write.csv(table_multi, 'Table2_Multivariate_Cox.csv', row.names = FALSE)
cat("Table 2 saved: Table2_Multivariate_Cox.csv\n")

# ============================================================
# 12. Summary report
# ============================================================
cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("COX REGRESSION ANALYSIS COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("OK Data preprocessing complete\n")
cat("OK Univariate Cox analysis complete (8 variables, including GRS)\n")
cat("OK Multivariate Cox analysis complete (5 variables, including GRS)\n")
cat("OK SCI format figures and tables generated\n\n")

cat("Generated files:\n")
cat("1. GDC_TCGA_LUAD_clinical_GRS.csv      - Processed data\n")
cat("2. GDC_TCGA_LUAD_clinical_GRS.Rdata    - R format data\n")
cat("3. Unicox.csv                          - Univariate Cox results\n")
cat("4. mulcox_res.csv                      - Multivariate Cox results\n")
cat("5. Figure_S1_Univariate_Cox.tiff       - Univariate forest plot (500 dpi)\n")
cat("6. Figure_S1_Univariate_Cox.png        - Univariate forest plot (preview)\n")
cat("7. Figure_S1_Univariate_Cox.pdf        - Univariate forest plot (vector)\n")
cat("8. Figure_2_Multivariate_Cox.tiff      - Multivariate forest plot (500 dpi)\n")
cat("9. Figure_2_Multivariate_Cox.png       - Multivariate forest plot (preview)\n")
cat("10. Figure_2_Multivariate_Cox.pdf      - Multivariate forest plot (vector)\n")
cat("11. Supplementary_Table1_Univariate_Cox.csv - Supplementary Table 1\n")
cat("12. Table2_Multivariate_Cox.csv        - Table 2 (main results)\n\n")

cat("Figure parameters:\n")
cat("OK Figure type: Forest plots\n")
cat("OK TIFF resolution: 500 dpi (Elsevier compliant)\n")
cat("OK PNG resolution: 300 dpi (preview)\n")
cat("OK PDF format: Vector (scale-free)\n")
cat("OK Compression: LZW lossless for TIFF\n\n")

cat("Variable names (Lung Cancer journal compliant):\n")
cat("  - Age\n")
cat("  - Gender\n")
cat("  - AJCC stage\n")
cat("  - T category\n")
cat("  - N category\n")
cat("  - M category\n")
cat("  - Smoking history\n")
cat("  - Genetic risk score (GRS)\n\n")

cat("Analysis completed at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 60), collapse = ""), "\n")