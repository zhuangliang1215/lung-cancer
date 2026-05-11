# ============================================================================
# LASSO_Cox_Survival_Analysis_SCI.R
# 
# Description: LASSO-Cox regression for gene selection and risk score calculation
#              SCI-level visualization (600 DPI, 8.5x7 cm figures)
# 
# Output: LASSO path plot, CV plot, risk score boxplot, ROC curve, jitter plot
# ============================================================================

# Load required packages
library(glmnet)
library(survival)
library(ggplot2)
library(ggpubr)
library(ROCR)
library(caret)

# ============================================================
# Part 1: Data preparation
# ============================================================
cat("[1/9] Loading and preparing data...\n")
cat("----------------------------------------\n")

rt <- read.csv("datExpr.csv", header = TRUE, row.names = 1, check.names = FALSE)

cat(sprintf("  Data dimensions: %d rows x %d columns\n", nrow(rt), ncol(rt)))

rt$status <- ifelse(rt$status == "Dead", 1, 0)
cat(sprintf("  Survival status: %d dead, %d alive\n", 
            sum(rt$status == 1), sum(rt$status == 0)))

# ============================================================
# Part 2: LASSO-Cox regression model
# ============================================================
cat("\n[2/9] Building LASSO-Cox regression model...\n")
cat("----------------------------------------\n")

set.seed(2)

x <- as.matrix(rt[, c(3:ncol(rt))])
y <- data.matrix(Surv(rt$time, rt$status))

cat(sprintf("  Feature matrix: %d samples x %d genes\n", nrow(x), ncol(x)))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
cvfit <- cv.glmnet(x, y, family = "cox", maxit = 1000, nfolds = 10)

coef_min <- coef(fit, s = cvfit$lambda.min)
selected_idx <- which(coef_min != 0)
selected_genes <- rownames(coef_min)[selected_idx]
actCoef <- coef_min[selected_idx]
n_genes <- length(selected_genes)

cat(sprintf("  Selected genes: %d\n", n_genes))
cat(sprintf("  lambda.min: %.4f\n", cvfit$lambda.min))
cat(sprintf("  lambda.1se: %.4f\n", cvfit$lambda.1se))

# ============================================================
# Part 3: SCI-level LASSO path plot
# ============================================================
cat("\n[3/9] Generating SCI-level LASSO path plot...\n")
cat("----------------------------------------\n")

lambda_values <- log(fit$lambda)
nonzero_counts <- apply(coef(fit), 2, function(x) sum(x != 0))
coef_matrix <- as.matrix(coef(fit))

# Smart selection of annotation positions
if (length(lambda_values) > 5) {
  change_points <- which(diff(nonzero_counts) != 0)
  if (length(change_points) >= 4) {
    key_positions <- unique(sort(c(1, change_points, length(lambda_values))))
    if (length(key_positions) > 6) {
      lambda_diff <- diff(lambda_values[key_positions])
      large_gaps <- which(lambda_diff > quantile(lambda_diff, 0.6))
      key_positions <- key_positions[unique(c(1, large_gaps + 1, length(key_positions)))]
    }
  } else {
    key_positions <- round(seq(1, length(lambda_values), length.out = 5))
  }
  key_positions <- head(key_positions, 5)
} else {
  key_positions <- 1:length(lambda_values)
}

lambda_key <- lambda_values[key_positions]
nonzero_key <- nonzero_counts[key_positions]

# SCI color palette
path_colors <- c(
  "#1E88E5", "#D32F2F", "#388E3C", "#7B1FA2", "#F57C00", "#0097A7",
  "#C2185B", "#00796B", "#512DA8", "#FBC02D", "#5D4037"
)

# Calculate Y-axis range
coef_range <- range(coef_matrix, na.rm = TRUE)
y_min <- floor(coef_range[1] * 1.1 / 0.05) * 0.05
y_max <- ceiling(coef_range[2] * 1.1 / 0.05) * 0.05

if (y_min > -0.05) y_min <- -0.1
if (y_max < 0.05) y_max <- 0.1
if (abs(y_min) > abs(y_max)) y_max <- abs(y_min)
if (abs(y_max) > abs(y_min)) y_min <- -abs(y_max)

# PDF version
pdf("Fig1_LASSO_path_SCI.pdf", 
    width = 8.5/2.54,
    height = 7/2.54,
    pointsize = 8,
    useDingbats = TRUE)

par(mar = c(3.8, 8.0, 3.5, 4.2),
    mgp = c(4.2, 0.8, 0),
    tcl = -0.25,
    cex.axis = 0.75,
    cex.lab = 0.85,
    las = 1,
    lwd = 0.8,
    family = "sans",
    bty = "l",
    xaxs = "i",
    yaxs = "i")

plot(NA, 
     xlim = range(lambda_values, na.rm = TRUE),
     ylim = c(y_min, y_max),
     xlab = "", ylab = "", main = "",
     bty = "n", xaxt = "n", yaxt = "n")

grid(col = "gray90", lty = "solid", lwd = 0.5)

num_paths <- min(nrow(coef_matrix), length(path_colors))
for(i in 1:num_paths) {
  lines(lambda_values, coef_matrix[i,], 
        col = path_colors[i], lwd = 0.8)
}

abline(h = 0, lty = "dashed", col = "gray50", lwd = 0.8)

lambda_min_x <- log(cvfit$lambda.min)
lambda_1se_x <- log(cvfit$lambda.1se)

abline(v = lambda_min_x, lty = "solid", lwd = 1.2, col = "#D32F2F")
if (abs(lambda_1se_x - lambda_min_x) > 0.08) {
  abline(v = lambda_1se_x, lty = "dashed", lwd = 1.0, col = "#1976D2")
}

# X-axis
x_ticks <- pretty(lambda_values, n = 6)
axis(1, at = x_ticks, 
     labels = format(round(x_ticks, 2), nsmall = 2),
     cex.axis = 0.75, col.axis = "#333333", col = "#666666", 
     lwd = 0.5, tcl = -0.25, lwd.ticks = 0.5, padj = -0.8, line = -0.5)

mtext(expression(paste("log(", lambda, ")")), 
      side = 1, line = 1.8, cex = 0.85, col = "#333333")

# Y-axis
y_ticks <- seq(y_min, y_max, by = 0.05)
y_labels <- format(y_ticks, digits = 2, nsmall = 2)
axis(2, at = y_ticks, labels = y_labels,
     cex.axis = 0.75, col.axis = "#333333", col = "#666666", 
     lwd = 0.5, tcl = -0.25, lwd.ticks = 0.5, hadj = 0.8, las = 1)

mtext("Coefficients", side = 2, line = 2.5, cex = 0.85, col = "#333333")

# Top annotation
mtext("Number of non-zero coefficients", 
      side = 3, line = 2.2, cex = 0.75, font = 2, col = "#333333")

top_margin <- par("usr")[4]
text_y_pos <- top_margin + (top_margin - par("usr")[3]) * 0.07

if (length(lambda_key) > 0) {
  adjusted_positions <- lambda_key
  if (length(lambda_key) > 1) {
    for (i in 2:length(lambda_key)) {
      if (lambda_key[i] - lambda_key[i-1] < 0.15) {
        adjusted_positions[i] <- lambda_key[i-1] + 0.15
      }
    }
  }
  
  for (i in 1:length(adjusted_positions)) {
    if (adjusted_positions[i] > par("usr")[1] && 
        adjusted_positions[i] < par("usr")[2]) {
      text(adjusted_positions[i], text_y_pos, nonzero_key[i],
           cex = 0.65, col = "#1E88E5", font = 2, xpd = TRUE)
    }
  }
}

# Lambda annotations
lambda_y_pos <- top_margin * 0.92
if (lambda_min_x > par("usr")[1] && lambda_min_x < par("usr")[2]) {
  text(lambda_min_x, lambda_y_pos, 
       expression(bold(paste(lambda[min]))),
       col = "#D32F2F", cex = 0.7, pos = 3, font = 2, xpd = TRUE)
}

if (abs(lambda_1se_x - lambda_min_x) > 0.08 && 
    lambda_1se_x > par("usr")[1] && 
    lambda_1se_x < par("usr")[2]) {
  text(lambda_1se_x, lambda_y_pos, 
       expression(bold(paste(lambda["1se"]))),
       col = "#1976D2", cex = 0.7, pos = 3, font = 2, xpd = TRUE)
}

# Selected genes count
text_x_pos <- par("usr")[2] - diff(par("usr")[1:2]) * 0.02
text_y_pos <- par("usr")[3] + diff(par("usr")[3:4]) * 0.08
text(text_x_pos, text_y_pos,
     sprintf("Selected: %d", n_genes),
     col = "#388E3C", cex = 0.7, pos = 2, font = 2, xpd = TRUE)

box(lwd = 0.8, col = "#666666")
dev.off()

cat("  Fig1_LASSO_path_SCI.pdf saved\n")

# PNG version
png("Fig1_LASSO_path_SCI.png",
    width = 8.5/2.54 * 600,
    height = 7/2.54 * 600,
    res = 600,
    pointsize = 8)

par(mar = c(3.8, 8.0, 3.5, 4.2),
    mgp = c(4.2, 0.8, 0),
    tcl = -0.25,
    cex.axis = 0.75,
    cex.lab = 0.85,
    las = 1,
    lwd = 0.8,
    family = "sans",
    bty = "l",
    xaxs = "i",
    yaxs = "i")

plot(NA, xlim = range(lambda_values), ylim = c(y_min, y_max),
     xlab = "", ylab = "", main = "", bty = "n", xaxt = "n", yaxt = "n")

grid(col = "gray90", lty = "solid", lwd = 0.5)

for(i in 1:num_paths) {
  lines(lambda_values, coef_matrix[i,], col = path_colors[i], lwd = 0.8)
}

abline(h = 0, lty = "dashed", col = "gray50", lwd = 0.8)
abline(v = lambda_min_x, lty = "solid", lwd = 1.2, col = "#D32F2F")
if (abs(lambda_1se_x - lambda_min_x) > 0.08) {
  abline(v = lambda_1se_x, lty = "dashed", lwd = 1.0, col = "#1976D2")
}

axis(1, at = x_ticks, labels = format(round(x_ticks, 2), nsmall = 2),
     cex.axis = 0.75, col.axis = "#333333", col = "#666666", 
     lwd = 0.5, tcl = -0.25, lwd.ticks = 0.5, padj = -0.8, line = -0.5)

axis(2, at = y_ticks, labels = y_labels,
     cex.axis = 0.75, col.axis = "#333333", col = "#666666", 
     lwd = 0.5, tcl = -0.25, lwd.ticks = 0.5, hadj = 0.8, las = 1)

mtext(expression(paste("log(", lambda, ")")), side = 1, line = 1.8, cex = 0.85, col = "#333333")
mtext("Coefficients", side = 2, line = 2.5, cex = 0.85, col = "#333333")
mtext("Number of non-zero coefficients", side = 3, line = 2.2, cex = 0.75, font = 2, col = "#333333")

top_margin <- par("usr")[4]
text_y_pos <- top_margin + (top_margin - par("usr")[3]) * 0.07

if (length(lambda_key) > 0) {
  for (i in 1:length(lambda_key)) {
    if (lambda_key[i] > par("usr")[1] && lambda_key[i] < par("usr")[2]) {
      text(lambda_key[i], text_y_pos, nonzero_key[i],
           cex = 0.65, col = "#1E88E5", font = 2, xpd = TRUE)
    }
  }
}

lambda_y_pos <- top_margin * 0.92
if (lambda_min_x > par("usr")[1] && lambda_min_x < par("usr")[2]) {
  text(lambda_min_x, lambda_y_pos, expression(bold(paste(lambda[min]))),
       col = "#D32F2F", cex = 0.7, pos = 3, font = 2, xpd = TRUE)
}

text_x_pos <- par("usr")[2] - diff(par("usr")[1:2]) * 0.02
text_y_pos <- par("usr")[3] + diff(par("usr")[3:4]) * 0.08
text(text_x_pos, text_y_pos, sprintf("Selected: %d", n_genes),
     col = "#388E3C", cex = 0.7, pos = 2, font = 2, xpd = TRUE)

box(lwd = 0.8, col = "#666666")
dev.off()

cat("  Fig1_LASSO_path_SCI.png saved\n")

# ============================================================
# Part 4: SCI-level cross-validation plot
# ============================================================
cat("\n[4/9] Generating SCI-level cross-validation plot...\n")
cat("----------------------------------------\n")

# Define color scheme
cv_line_color <- "#1E88E5"
cv_point_color <- "#D32F2F"
cv_band_color <- "#90CAF9"
cv_lambda_min_color <- "#D32F2F"
cv_lambda_1se_color <- "#1976D2"
cv_grid_color <- "#F5F5F5"

# PDF version
pdf("Fig2_CV_plot_SCI.pdf", 
    width = 8.5/2.54,
    height = 7/2.54,
    pointsize = 8,
    useDingbats = TRUE)

par(mar = c(3.8, 5.5, 3.5, 2.0),
    mgp = c(3.2, 0.7, 0),
    tcl = -0.25,
    cex.axis = 0.75,
    cex.lab = 0.85,
    las = 0,
    lwd = 0.8,
    family = "sans",
    bty = "l",
    xaxs = "i",
    yaxs = "i")

plot(cvfit, lwd = 1.2, col = cv_line_color, pch = 16, cex = 0.6,
     xlab = "", ylab = "", main = "", bty = "n", axes = FALSE)

grid(col = cv_grid_color, lty = "solid", lwd = 0.5)

lambda_min_x_cv <- log(cvfit$lambda.min)
lambda_1se_x_cv <- log(cvfit$lambda.1se)

abline(v = lambda_min_x_cv, lty = "solid", lwd = 1.5, col = cv_lambda_min_color)
if (abs(lambda_1se_x_cv - lambda_min_x_cv) > 0.08) {
  abline(v = lambda_1se_x_cv, lty = "dashed", lwd = 1.2, col = cv_lambda_1se_color)
}

if (!is.null(cvfit$cvup) && !is.null(cvfit$cvlo)) {
  polygon(c(rev(log(cvfit$lambda)), log(cvfit$lambda)), 
          c(rev(cvfit$cvup), cvfit$cvlo), 
          col = adjustcolor(cv_band_color, alpha.f = 0.3),
          border = NA)
}

lines(log(cvfit$lambda), cvfit$cvm, lwd = 1.2, col = cv_line_color)
points(log(cvfit$lambda), cvfit$cvm, pch = 16, col = cv_line_color, cex = 0.6)
points(lambda_min_x_cv, min(cvfit$cvm), 
       pch = 21, bg = cv_point_color, col = cv_point_color, cex = 0.9, lwd = 1.0)

x_ticks_cv <- pretty(log(cvfit$lambda), n = 6)
axis(1, at = x_ticks_cv, 
     labels = format(round(x_ticks_cv, 2), nsmall = 2),
     cex.axis = 0.75, col.axis = "#333333", col = "#666666", 
     lwd = 0.5, tcl = -0.25, lwd.ticks = 0.5, padj = -0.8, line = -0.5)

mtext(expression(paste("log(", lambda, ")")), side = 1, line = 1.8, cex = 0.85, col = "#333333")

y_range_cv <- range(cvfit$cvm)
y_pretty_cv <- pretty(y_range_cv, n = 6)
axis(2, at = y_pretty_cv, 
     labels = format(round(y_pretty_cv, 3), nsmall = 3),
     cex.axis = 0.75, col.axis = "#333333", col = "#666666", 
     lwd = 0.5, tcl = -0.25, lwd.ticks = 0.5, hadj = 0.8, las = 1)

mtext("Partial Likelihood Deviance", side = 2, line = 2.8, cex = 0.85, col = "#333333")

if (length(cvfit$lambda) > 0 && length(cvfit$nzero) > 0) {
  n_indices <- min(8, length(cvfit$lambda))
  step <- max(1, floor(length(cvfit$lambda) / n_indices))
  nzero_key_indices <- seq(1, length(cvfit$lambda), by = step)
  
  mtext("Number of non-zero coefficients", side = 3, line = 2.2, cex = 0.75, font = 2, col = "#333333")
  
  top_margin_cv <- par("usr")[4]
  text_y_pos_cv <- top_margin_cv + (top_margin_cv - par("usr")[3]) * 0.07
  
  for (idx in nzero_key_indices) {
    if (idx <= length(cvfit$lambda)) {
      lambda_val <- log(cvfit$lambda[idx])
      nzero_val <- cvfit$nzero[idx]
      if (lambda_val > par("usr")[1] && lambda_val < par("usr")[2]) {
        text(lambda_val, text_y_pos_cv, nzero_val,
             cex = 0.6, col = "#1E88E5", font = 2, xpd = TRUE)
      }
    }
  }
}

lambda_y_pos_cv <- par("usr")[4] * 0.92
text(lambda_min_x_cv, lambda_y_pos_cv, 
     expression(bold(paste(lambda[min]))),
     col = cv_lambda_min_color, cex = 0.7, pos = 3, font = 2, xpd = TRUE)

if (abs(lambda_1se_x_cv - lambda_min_x_cv) > 0.08) {
  text(lambda_1se_x_cv, lambda_y_pos_cv, 
       expression(bold(paste(lambda["1se"]))),
       col = cv_lambda_1se_color, cex = 0.7, pos = 3, font = 2, xpd = TRUE)
}

box(lwd = 0.8, col = "#666666")
dev.off()

cat("  Fig2_CV_plot_SCI.pdf saved\n")

# PNG version
png("Fig2_CV_plot_SCI.png",
    width = 8.5/2.54 * 600,
    height = 7/2.54 * 600,
    res = 600,
    pointsize = 8)

par(mar = c(3.8, 5.5, 3.5, 2.0),
    mgp = c(3.2, 0.7, 0),
    tcl = -0.25,
    cex.axis = 0.75,
    cex.lab = 0.85,
    las = 0,
    lwd = 0.8,
    family = "sans",
    bty = "l",
    xaxs = "i",
    yaxs = "i")

plot(cvfit, lwd = 1.2, col = cv_line_color, pch = 16, cex = 0.6,
     xlab = "", ylab = "", main = "", bty = "n", axes = FALSE)

grid(col = cv_grid_color, lty = "solid", lwd = 0.5)
abline(v = lambda_min_x_cv, lty = "solid", lwd = 1.5, col = cv_lambda_min_color)
if (abs(lambda_1se_x_cv - lambda_min_x_cv) > 0.08) {
  abline(v = lambda_1se_x_cv, lty = "dashed", lwd = 1.2, col = cv_lambda_1se_color)
}

if (!is.null(cvfit$cvup) && !is.null(cvfit$cvlo)) {
  polygon(c(rev(log(cvfit$lambda)), log(cvfit$lambda)), 
          c(rev(cvfit$cvup), cvfit$cvlo), 
          col = adjustcolor(cv_band_color, alpha.f = 0.3), border = NA)
}

lines(log(cvfit$lambda), cvfit$cvm, lwd = 1.2, col = cv_line_color)
points(log(cvfit$lambda), cvfit$cvm, pch = 16, col = cv_line_color, cex = 0.6)
points(lambda_min_x_cv, min(cvfit$cvm), 
       pch = 21, bg = cv_point_color, col = cv_point_color, cex = 0.9, lwd = 1.0)

axis(1, at = x_ticks_cv, labels = format(round(x_ticks_cv, 2), nsmall = 2),
     cex.axis = 0.75, col.axis = "#333333", col = "#666666", 
     lwd = 0.5, tcl = -0.25, lwd.ticks = 0.5, padj = -0.8, line = -0.5)

axis(2, at = y_pretty_cv, labels = format(round(y_pretty_cv, 3), nsmall = 3),
     cex.axis = 0.75, col.axis = "#333333", col = "#666666", 
     lwd = 0.5, tcl = -0.25, lwd.ticks = 0.5, hadj = 0.8, las = 1)

mtext(expression(paste("log(", lambda, ")")), side = 1, line = 1.8, cex = 0.85, col = "#333333")
mtext("Partial Likelihood Deviance", side = 2, line = 2.8, cex = 0.85, col = "#333333")

if (length(cvfit$lambda) > 0 && length(cvfit$nzero) > 0) {
  mtext("Number of non-zero coefficients", side = 3, line = 2.2, cex = 0.75, font = 2, col = "#333333")
  top_margin_cv <- par("usr")[4]
  text_y_pos_cv <- top_margin_cv + (top_margin_cv - par("usr")[3]) * 0.07
  
  for (idx in nzero_key_indices) {
    if (idx <= length(cvfit$lambda)) {
      lambda_val <- log(cvfit$lambda[idx])
      nzero_val <- cvfit$nzero[idx]
      if (lambda_val > par("usr")[1] && lambda_val < par("usr")[2]) {
        text(lambda_val, text_y_pos_cv, nzero_val,
             cex = 0.6, col = "#1E88E5", font = 2, xpd = TRUE)
      }
    }
  }
}

text(lambda_min_x_cv, lambda_y_pos_cv, expression(bold(paste(lambda[min]))),
     col = cv_lambda_min_color, cex = 0.7, pos = 3, font = 2, xpd = TRUE)

box(lwd = 0.8, col = "#666666")
dev.off()

cat("  Fig2_CV_plot_SCI.png saved\n")

# ============================================================
# Part 5: Calculate risk score
# ============================================================
cat("\n[5/9] Calculating risk score...\n")
cat("----------------------------------------\n")

if (n_genes > 0) {
  geneCoef <- cbind(Gene = selected_genes, Coef = actCoef)
  write.csv(geneCoef, 'geneCoef.csv', row.names = TRUE)
  
  FinalGeneExp <- rt[, selected_genes, drop = FALSE]
  myFun <- function(x){crossprod(as.numeric(x), actCoef)}
  riskScore <- apply(FinalGeneExp, 1, myFun)
  risk <- as.vector(ifelse(riskScore > median(riskScore), "high", "low"))
  
  dat <- cbind(rt[, c("time", "status", selected_genes)], 
               riskScore = as.vector(riskScore), 
               risk = risk,
               status_label = ifelse(rt$status == 1, "Dead", "Alive"))
  
  write.csv(dat, 'dat.csv', row.names = TRUE)
  
  cat(sprintf("  Risk score range: %.3f - %.3f\n", min(riskScore), max(riskScore)))
  cat(sprintf("  Median risk score: %.3f\n", median(riskScore)))
  cat(sprintf("  High-risk group: %d\n", sum(risk == "high")))
  cat(sprintf("  Low-risk group: %d\n", sum(risk == "low")))
}

# ============================================================
# Part 6: SCI-level risk score boxplot
# ============================================================
cat("\n[6/9] Generating SCI-level risk score boxplot...\n")
cat("----------------------------------------\n")

box_alive_color <- "#1E88E5"
box_dead_color <- "#D32F2F"

p_val <- wilcox.test(riskScore ~ status, data = dat)$p.value
p_label <- ifelse(p_val < 0.001, "p < 0.001", 
                  ifelse(p_val < 0.01, "p < 0.01",
                         ifelse(p_val < 0.05, "p < 0.05",
                                sprintf("p = %.3f", p_val))))

p_boxplot <- ggplot(dat, aes(x = factor(status_label, levels = c("Alive", "Dead")), 
                             y = riskScore, fill = status_label)) +
  geom_boxplot(width = 0.5, alpha = 0.85, outlier.shape = NA, lwd = 0.25) +
  geom_jitter(width = 0.1, height = 0, alpha = 0.6, size = 1.0, 
              aes(color = status_label)) +
  scale_fill_manual(values = c(box_alive_color, box_dead_color)) +
  scale_color_manual(values = c(box_alive_color, box_dead_color)) +
  labs(x = "Survival Status", y = "Risk Score", 
       subtitle = paste0("Wilcoxon test: ", p_label)) +
  theme_classic(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        axis.text = element_text(color = "black", size = 7),
        axis.title = element_text(color = "black", size = 8, face = "bold"),
        legend.position = "none",
        plot.subtitle = element_text(size = 7, color = "black", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5)) +
  annotate("segment", x = 1, xend = 2, y = max(dat$riskScore) * 1.05, 
           yend = max(dat$riskScore) * 1.05, color = "black", size = 0.25) +
  annotate("text", x = 1.5, y = max(dat$riskScore) * 1.08, 
           label = p_label, size = 2.5, color = "black")

ggsave("Fig3_RiskScore_boxplot_SCI.pdf", p_boxplot, 
       width = 8.5/2.54, height = 7/2.54, dpi = 600)
ggsave("Fig3_RiskScore_boxplot_SCI.png", p_boxplot, 
       width = 8.5/2.54, height = 7/2.54, dpi = 600)

cat(sprintf("  Fig3_RiskScore_boxplot_SCI saved (p = %.4f)\n", p_val))

# ============================================================
# Part 7: SCI-level ROC curve
# ============================================================
cat("\n[7/9] Generating SCI-level ROC curve...\n")
cat("----------------------------------------\n")

pred <- prediction(dat$riskScore, dat$status)
perf <- performance(pred, "tpr", "fpr")
auc_value <- performance(pred, "auc")@y.values[[1]]

roc_data <- data.frame(FPR = perf@x.values[[1]], TPR = perf@y.values[[1]])
roc_line_color <- "#7B1FA2"

p_roc <- ggplot(roc_data, aes(x = FPR, y = TPR)) +
  geom_line(color = roc_line_color, size = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "gray60", size = 0.25, alpha = 0.7) +
  labs(x = "False Positive Rate (1 - Specificity)",
       y = "True Positive Rate (Sensitivity)",
       subtitle = sprintf("AUC = %.3f", auc_value)) +
  theme_classic(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        axis.text = element_text(color = "black", size = 7),
        axis.title = element_text(color = "black", size = 8, face = "bold"),
        legend.position = "none",
        plot.subtitle = element_text(size = 7, color = "black", hjust = 0.5),
        plot.margin = margin(5, 5, 5, 5)) +
  coord_equal() +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text", x = 0.6, y = 0.2, 
           label = sprintf("AUC = %.3f", auc_value),
           size = 3.0, fontface = "bold", color = roc_line_color)

ggsave("Fig4_ROC_curve_SCI.pdf", p_roc, width = 8.5/2.54, height = 7/2.54, dpi = 600)
ggsave("Fig4_ROC_curve_SCI.png", p_roc, width = 8.5/2.54, height = 7/2.54, dpi = 600)

cat(sprintf("  Fig4_ROC_curve_SCI saved (AUC = %.3f)\n", auc_value))

# ============================================================
# Part 8: SCI-level jitter plot
# ============================================================
cat("\n[8/9] Generating SCI-level jitter plot...\n")
cat("----------------------------------------\n")

p_jitter <- ggpubr::ggboxplot(dat, x = "status_label", y = "riskScore", 
                              color = "status_label", 
                              palette = c(box_alive_color, box_dead_color),
                              add = "jitter", size = 0.25) +
  ggpubr::stat_compare_means(size = 2.5, label.y = max(dat$riskScore) * 1.1) +
  labs(x = "Survival Status", y = "Risk Score") +
  theme_classic(base_size = 8) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black", linewidth = 0.25),
        axis.ticks = element_line(color = "black", linewidth = 0.25),
        axis.text = element_text(color = "black", size = 7),
        axis.title = element_text(color = "black", size = 8, face = "bold"),
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5))

ggsave("Fig5_Jitter_plot_SCI.pdf", p_jitter, width = 8.5/2.54, height = 7/2.54, dpi = 600)
ggsave("Fig5_Jitter_plot_SCI.png", p_jitter, width = 8.5/2.54, height = 7/2.54, dpi = 600)

cat("  Fig5_Jitter_plot_SCI saved\n")

# ============================================================
# Part 9: Summary report
# ============================================================
cat("\n[9/9] SCI-level results summary...\n")
cat("========================================\n")
cat("\nLASSO-Cox Survival Analysis Results\n")
cat("----------------------------------------\n")
cat(sprintf("  Total samples: %d\n", nrow(rt)))
cat(sprintf("  Total features: %d\n", ncol(x)))
cat(sprintf("  Selected genes: %d\n", n_genes))
cat(sprintf("  lambda.min: %.4f (log: %.4f)\n", cvfit$lambda.min, log(cvfit$lambda.min)))
cat(sprintf("  lambda.1se: %.4f (log: %.4f)\n", cvfit$lambda.1se, log(cvfit$lambda.1se)))
cat(sprintf("  ROC AUC: %.3f\n", auc_value))
cat(sprintf("  Risk score p-value: %.4f\n", p_val))

cat("\nGenerated files:\n")
cat("----------------------------------------\n")
files_list <- c(
  "Fig1_LASSO_path_SCI.pdf",
  "Fig1_LASSO_path_SCI.png",
  "Fig2_CV_plot_SCI.pdf",
  "Fig2_CV_plot_SCI.png",
  "Fig3_RiskScore_boxplot_SCI.pdf",
  "Fig3_RiskScore_boxplot_SCI.png",
  "Fig4_ROC_curve_SCI.pdf",
  "Fig4_ROC_curve_SCI.png",
  "Fig5_Jitter_plot_SCI.pdf",
  "Fig5_Jitter_plot_SCI.png",
  "geneCoef.csv",
  "dat.csv"
)

for (file in files_list) {
  if (file.exists(file)) {
    file_size <- file.info(file)$size / 1024
    cat(sprintf("  OK %-35s %6.1f KB\n", file, file_size))
  }
}

cat("\n========================================\n")
cat("LASSO-Cox Survival Analysis Complete\n")
cat("All figures are optimized for SCI submission (600 DPI, 8.5x7 cm)\n")
cat("========================================\n")

sink("session_info.txt")
sessionInfo()
sink()
cat("\n  Session info saved: session_info.txt\n")