# ============================================================================
# KM_curve_analysis.R
# 
# Description: Kaplan-Meier survival curve analysis with log-rank test
#              High-resolution output for SCI publication (500 dpi)
# 
# Output: Survival curve with risk table, censoring plot, TIFF/PNG files
# ============================================================================

# ================================
# Set working directory and load packages
# ================================
setwd("C:/Users/Administrator/Desktop/KM")

library(survival)   # Survival analysis
library(survminer)  # Kaplan-Meier plotting

# ================================
# Read data
# ================================
dat <- read.csv("dat.csv", header = TRUE, row.names = 1, check.names = FALSE)

# Convert time to months (assuming original unit is days)
dat$time <- dat$time / 30

# Optional: Convert risk factor to numeric
# dat$risk <- ifelse(dat$risk == "high", 1, 0)

# ================================
# KM survival fit
# ================================
fit <- survfit(Surv(time, status) ~ risk, data = dat)

print(fit)
summary(fit)

# ================================
# Extract survival table data
# ================================
d <- data.frame(
  time    = fit$time,
  n.risk  = fit$n.risk,
  n.event = fit$n.event,
  n.censor = fit$n.censor,
  surv    = fit$surv,
  upper   = fit$upper,
  lower   = fit$lower
)
head(d)

# ================================
# Log-rank test
# ================================
sur_diff <- survdiff(Surv(time, status) ~ risk, data = dat)
print(sur_diff)

# ================================
# High-quality KM curve (SCI standard, high DPI)
# ================================

# Define colors (customizable)
my_colors <- c("#E7B800", "#2E9FDF")

# Create KM plot
km_plot <- ggsurvplot(
  fit,
  pval = TRUE,               # Show log-rank P value
  conf.int = TRUE,           # Show confidence interval
  pval.method = TRUE,        # Show test method
  pval.size = 3.3,           # P value font size
  risk.table = TRUE,         # Show risk table
  risk.table.col = "strata", # Color risk table by group
  linetype = "strata",       # Different line types for groups
  surv.median.line = "hv",   # Show median survival time line
  ggtheme = theme_bw(),      # Base theme
  xlab = "Time (months)",    # X-axis label
  legend.labs = c("High risk", "Low risk"),  # Legend labels
  legend.title = "NSCLC",
  ncensor.plot = TRUE,       # Show censoring plot
  tables.height = 0.25,      # Risk table height ratio
  fontsize = 3,              # Risk table font size
  palette = my_colors        # Custom colors
)

# Display plot
print(km_plot)

# ================================
# Save as high-resolution TIFF (SCI recommended)
# ================================
tiff("KM_curve_SCI.tiff", 
     width = 8, height = 7, units = "in", 
     res = 500, compression = "lzw")

print(km_plot)
dev.off()

# Save as high-resolution PNG
png("KM_curve_SCI.png", width = 8, height = 7, units = "in", res = 500)
print(km_plot)
dev.off()

# ================================
# Optional: Example using lung dataset
# ================================
if (FALSE) {  # Not run by default, example only
  fit_lung <- survfit(Surv(time, status) ~ sex, data = lung)
  sur_diff_lung <- survdiff(Surv(time, status) ~ sex, data = lung)
  print(sur_diff_lung)
  
  col_lung <- c("#E7B800", "#2E9FDF")
  
  km_lung <- ggsurvplot(
    fit_lung,
    xlim = c(0, 1100),
    pval = TRUE,
    conf.int = TRUE,
    conf.int.style = "ribbon",
    linetype = "strata",
    surv.median.line = "hv",
    ggtheme = theme_classic(),
    legend.labs = c("F", "M"),
    legend.title = "Sex",
    legend = c(0.80, 0.88),
    palette = col_lung
  )
  
  tiff("KM_curve_lung_SCI.tiff", width = 8, height = 7, units = "in", res = 500)
  print(km_lung)
  dev.off()
}

# ================================
# Summary report
# ================================
cat("\n", rep("=", 60), "\n", sep = "")
cat("KAPLAN-MEIER SURVIVAL ANALYSIS COMPLETE\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Generated files:\n")
cat("  KM_curve_SCI.tiff - 500 dpi, TIFF format (for submission)\n")
cat("  KM_curve_SCI.png   - 500 dpi, PNG format (for preview)\n\n")

cat("Figure parameters:\n")
cat("  Dimensions: 8 x 7 inches\n")
cat("  Resolution: 500 dpi\n")
cat("  Compression: LZW (lossless)\n\n")

cat("Statistical results:\n")
print(sur_diff)

cat("\nAnalysis completed at: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(rep("=", 60), "\n", sep = "")