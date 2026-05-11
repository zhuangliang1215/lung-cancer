# ============================================================================
# PET_metabolic_Cox_univariate_multivariate.R
# 
# Description: Univariate and multivariate Cox regression analysis for 
#              PET/CT metabolic parameters (SUVmean, SUVmax, MTV, TLG)
# 
# Output: Forest plots (TIFF 500 dpi for submission, PNG 300 dpi for preview)
#         Results tables (CSV format)
# ============================================================================

# Load required packages
library("R.utils")
library(data.table)
library(survival)
library(plyr)
library(forestplot)
library(grid)

# Set working directory
setwd("C:/Users/Administrator/Desktop/Pet1")

# ============================================================
# Part 1: Data preprocessing
# ============================================================
cat("Starting data preprocessing...\n")

# Read clinical data
clin <- read.csv("pet.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Convert survival status: 'Alive' -> 0, others (Dead) -> 1
clin$status <- ifelse(clin$status == 'Alive', 0, 1)

# Save processed data
write.csv(clin, 'PET_clinical_data_processed.csv', row.names = TRUE)
save(clin, file = 'PET_clinical_data_processed.Rdata')

cat("Data preprocessing complete\n")

# ============================================================
# Part 2: Univariate Cox regression analysis
# ============================================================
cat("\nStarting Univariate Cox regression analysis...\n")

# Load data
load(file = 'PET_clinical_data_processed.Rdata')

# Create survival object
y <- Surv(clin$time, clin$status)

# Univariate Cox analysis function
Unicox_model <- function(x){
  # Handle variable names with spaces
  if(grepl(" ", x)) {
    x_clean <- paste0("`", x, "`")
  } else {
    x_clean <- x
  }
  
  FML <- as.formula(paste0("y ~ ", x_clean))
  cox <- coxph(FML, data = clin)
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

# Define PET metabolic parameters for analysis
variable <- c("SUVmean", "SUVmax", "MTV", "TLG")

# Run univariate analysis
Unicox <- lapply(variable, Unicox_model)
Unicox <- ldply(Unicox, data.frame)
colnames(Unicox) <- c("Variable", "HR", "CI5", "CI95", "HR (95% CI)", "Pvalue")

# Save univariate results
write.csv(Unicox, 'Univariate_Cox_PET_results.csv', row.names = FALSE)
cat("Univariate Cox analysis complete. Results saved to Univariate_Cox_PET_results.csv\n")

# Display univariate results
cat("\nUnivariate Cox regression results (PET parameters):\n")
print(Unicox)

# ============================================================
# Part 3: Multivariate Cox regression analysis
# ============================================================
cat("\nStarting Multivariate Cox regression analysis...\n")

# Run multivariate Cox analysis (selecting meaningful variables)
res <- coxph(Surv(time, status) ~ SUVmean + SUVmax, data = clin)
mul_cox <- summary(res)

# Extract multivariate results
mul_HR <- round(mul_cox$coefficients[,2], 2)
mul_Pvalue <- round(mul_cox$coefficients[,5], 4)
mul_CI5 <- round(mul_cox$conf.int[,3], 2)
mul_CI95 <- round(mul_cox$conf.int[,4], 2)
mul_CI <- paste0(mul_HR, ' (', mul_CI5, '-', mul_CI95, ')')

# Get and clean variable names
Variable <- row.names(data.frame(mul_cox$coefficients))
Variable_clean <- gsub("`", "", Variable)

mulcox_res <- data.frame(Variable = Variable_clean,
                         HR = mul_HR,
                         CI5 = mul_CI5,
                         CI95 = mul_CI95,
                         CI = mul_CI,
                         Pvalue = mul_Pvalue)

colnames(mulcox_res) <- c("Variable", "HR", "CI5", "CI95", "HR (95% CI)", "Pvalue")

# Save multivariate results
write.csv(mulcox_res, 'Multivariate_Cox_PET_results.csv', row.names = FALSE)
cat("Multivariate Cox analysis complete. Results saved to Multivariate_Cox_PET_results.csv\n")

# Display multivariate results
cat("\nMultivariate Cox regression results (PET parameters):\n")
print(mulcox_res)

# ============================================================
# Part 4: Unified figure settings (Elsevier compliant)
# ============================================================

# Fixed dimensions
FIG_WIDTH <- 5.76    # Width (inches)
FIG_HEIGHT <- 4.61   # Height (inches)

# Unified color scheme
UNIFIED_COLORS <- list(
  box = "#1c61b6",    # Professional blue
  lines = "#1c61b6",  # Line color
  zero = "#666666"    # Reference line color
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
# Part 5: Univariate forest plot (Supplementary material)
# ============================================================
cat("\nGenerating Univariate forest plot (Supplementary)...\n")

# Prepare table text
uni_tabletext <- cbind(
  c("Variable", as.character(Unicox$Variable)),
  c("HR (95% CI)", as.character(Unicox$`HR (95% CI)`)),
  c("P Value", ifelse(Unicox$Pvalue < 0.001, "<0.001", 
                      format(round(Unicox$Pvalue, 4), nsmall = 4)))
)

# Font settings
txt_gp <- fpTxtGp(
  label = gpar(cex = 0.9, fontfamily = "Helvetica", fontface = "plain"),
  ticks = gpar(cex = 0.8, fontfamily = "Helvetica"),
  xlab = gpar(cex = 1.0, fontfamily = "Helvetica", fontface = "bold"),
  title = gpar(cex = 1.1, fontfamily = "Helvetica", fontface = "bold")
)

# Horizontal line settings
hrzl_lines <- list(
  "1" = gpar(lty = 1, lwd = 1.5, col = "black"),
  "2" = gpar(lty = 1, lwd = 1, col = "black")
)
last_line <- as.character(nrow(Unicox) + 2)
hrzl_lines[[last_line]] <- gpar(lty = 1, lwd = 1.5, col = "black")

# 5.1 TIFF version (500 dpi for submission)
cat("  Generating TIFF version (500 dpi for submission)...\n")
tiff("Figure_S1_Univariate_Cox_PET.tiff", 
     width = FIG_WIDTH, 
     height = FIG_HEIGHT, 
     units = "in", 
     res = 500,
     compression = "lzw",
     pointsize = 10)

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
  txt_gp = txt_gp,
  hrzl_lines = hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(box = UNIFIED_COLORS$box, lines = UNIFIED_COLORS$lines, zero = UNIFIED_COLORS$zero),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Univariate Cox Regression Analysis (PET Parameters)"
)
dev.off()
cat("  OK Saved: Figure_S1_Univariate_Cox_PET.tiff (500 dpi)\n")

# 5.2 PNG version (300 dpi for preview)
cat("  Generating PNG version (300 dpi for preview)...\n")
png("Figure_S1_Univariate_Cox_PET.png", 
    width = FIG_WIDTH * 300, 
    height = FIG_HEIGHT * 300, 
    res = 300)

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
  txt_gp = txt_gp,
  hrzl_lines = hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(box = UNIFIED_COLORS$box, lines = UNIFIED_COLORS$lines, zero = UNIFIED_COLORS$zero),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Univariate Cox Regression Analysis (PET Parameters)"
)
dev.off()
cat("  OK Saved: Figure_S1_Univariate_Cox_PET.png (300 dpi)\n")

# ============================================================
# Part 6: Multivariate forest plot (Main figure)
# ============================================================
cat("\nGenerating Multivariate forest plot (Main figure)...\n")

# Prepare table text
multi_tabletext <- cbind(
  c("Variable", as.character(mulcox_res$Variable)),
  c("HR (95% CI)", as.character(mulcox_res$`HR (95% CI)`)),
  c("P Value", ifelse(mulcox_res$Pvalue < 0.001, "<0.001", 
                      format(round(mulcox_res$Pvalue, 4), nsmall = 4)))
)

# Multivariate horizontal line settings
multi_hrzl_lines <- list(
  "1" = gpar(lty = 1, lwd = 1.5, col = "black"),
  "2" = gpar(lty = 1, lwd = 1, col = "black")
)
multi_last_line <- as.character(nrow(mulcox_res) + 2)
multi_hrzl_lines[[multi_last_line]] <- gpar(lty = 1, lwd = 1.5, col = "black")

# 6.1 TIFF version (500 dpi for submission)
cat("  Generating TIFF version (500 dpi for submission)...\n")
tiff("Figure_1_Multivariate_Cox_PET.tiff", 
     width = FIG_WIDTH, 
     height = FIG_HEIGHT, 
     units = "in", 
     res = 500,
     compression = "lzw",
     pointsize = 10)

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
  txt_gp = txt_gp,
  hrzl_lines = multi_hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(box = UNIFIED_COLORS$box, lines = UNIFIED_COLORS$lines, zero = UNIFIED_COLORS$zero),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Multivariate Cox Regression Analysis (PET Parameters)"
)
dev.off()
cat("  OK Saved: Figure_1_Multivariate_Cox_PET.tiff (500 dpi)\n")

# 6.2 PNG version (300 dpi for preview)
cat("  Generating PNG version (300 dpi for preview)...\n")
png("Figure_1_Multivariate_Cox_PET.png", 
    width = FIG_WIDTH * 300, 
    height = FIG_HEIGHT * 300, 
    res = 300)

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
  txt_gp = txt_gp,
  hrzl_lines = multi_hrzl_lines,
  grid = structure(c(1), gp = gpar(lty = 2, lwd = 0.5, col = "#EFEFEF")),
  col = fpColors(box = UNIFIED_COLORS$box, lines = UNIFIED_COLORS$lines, zero = UNIFIED_COLORS$zero),
  lwd.zero = UNIFIED_SETTINGS$lwd.zero,
  lwd.ci = UNIFIED_SETTINGS$lwd.ci,
  ci.vertices = UNIFIED_SETTINGS$ci.vertices,
  ci.vertices.height = UNIFIED_SETTINGS$ci.vertices.height,
  title = "Multivariate Cox Regression Analysis (PET Parameters)"
)
dev.off()
cat("  OK Saved: Figure_1_Multivariate_Cox_PET.png (300 dpi)\n")

# ============================================================
# Part 7: Create SCI format tables
# ============================================================
cat("\nCreating SCI format tables...\n")

# Table S1: Univariate Cox analysis (Supplementary)
table_uni <- data.frame(
  Variable = Unicox$Variable,
  `HR (95% CI)` = Unicox$`HR (95% CI)`,
  `P Value` = ifelse(Unicox$Pvalue < 0.001, "<0.001", 
                     format(round(Unicox$Pvalue, 4), nsmall = 4))
)
write.csv(table_uni, 'Supplementary_Table_S1_Univariate_Cox_PET.csv', row.names = FALSE)
cat("  OK Saved: Supplementary_Table_S1_Univariate_Cox_PET.csv\n")

# Table 1: Multivariate Cox analysis (Main text)
table_multi <- data.frame(
  Variable = mulcox_res$Variable,
  `HR (95% CI)` = mulcox_res$`HR (95% CI)`,
  `P Value` = ifelse(mulcox_res$Pvalue < 0.001, "<0.001", 
                     format(round(mulcox_res$Pvalue, 4), nsmall = 4))
)
write.csv(table_multi, 'Table_1_Multivariate_Cox_PET.csv', row.names = FALSE)
cat("  OK Saved: Table_1_Multivariate_Cox_PET.csv\n")

# ============================================================
# Part 8: Verification and summary
# ============================================================
cat("\n")
cat(rep("=", 60), "\n", sep = "")
cat("FILE GENERATION COMPLETE\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Generated files:\n\n")

cat("[TIFF files for submission - 500 dpi]\n")
cat("  Figure_S1_Univariate_Cox_PET.tiff  (Supplementary)\n")
cat("  Figure_1_Multivariate_Cox_PET.tiff  (Main figure)\n\n")

cat("[PNG files for preview - 300 dpi]\n")
cat("  Figure_S1_Univariate_Cox_PET.png\n")
cat("  Figure_1_Multivariate_Cox_PET.png\n\n")

cat("[CSV tables]\n")
cat("  Supplementary_Table_S1_Univariate_Cox_PET.csv\n")
cat("  Table_1_Multivariate_Cox_PET.csv\n\n")

cat("[Analysis results]\n")
cat("  Univariate_Cox_PET_results.csv\n")
cat("  Multivariate_Cox_PET_results.csv\n")
cat("  PET_clinical_data_processed.csv\n")
cat("  PET_clinical_data_processed.Rdata\n\n")

cat(rep("=", 60), "\n", sep = "")
cat("Dimension verification:\n")
cat(rep("=", 60), "\n", sep = "")
cat(sprintf("Physical dimensions: %.2f in x %.2f in\n", FIG_WIDTH, FIG_HEIGHT))
cat(sprintf("TIFF (500 dpi) pixels: %.0f x %.0f\n", FIG_WIDTH * 500, FIG_HEIGHT * 500))
cat(sprintf("PNG (300 dpi) pixels: %.0f x %.0f\n", FIG_WIDTH * 300, FIG_HEIGHT * 300))
cat(rep("=", 60), "\n", sep = "")

cat("\nAll files generated successfully!\n")
cat("TIFF files (500 dpi) are ready for journal submission\n")
cat("PNG files (300 dpi) are ready for preview\n")
cat("CSV tables can be directly copied to Word\n\n")

cat("Analysis completed at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")