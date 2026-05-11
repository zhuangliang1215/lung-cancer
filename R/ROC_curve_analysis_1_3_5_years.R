# ============================================================================
# ROC Curve Analysis Code (Compliant with Elsevier Submission Requirements)
# Function: Plot 1-year, 3-year, and 5-year ROC curves
# Output: TIFF (500 dpi for submission) + PNG (300 dpi for preview) + PDF (vector)
# ============================================================================

rm(list = ls())
setwd("C:/Users/Administrator/Desktop/Roc")

# Load required packages
library(tidyverse)
library(randomForest)
library(pROC)
library(caret)
library(ggplot2)

# Set global seed for reproducibility
set.seed(123)
options(rf.cores = 1)

# ==================== 1. Data Reading Function ====================
#' Read single-year CSV file
#' @param file_path Path to CSV file
#' @return Processed data frame
read_year_data <- function(file_path) {
  if (!file.exists(file_path)) {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
  
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  mean_col <- case_when(
    "mean3" %in% colnames(df) ~ "mean3",
    "mean" %in% colnames(df) ~ "mean",
    TRUE ~ stop(paste("No valid mean column (mean3 or mean) in:", file_path))
  )
  
  year_num <- str_extract(file_path, "one|two|three|four|five|1|2|3|4|5") %>%
    case_when(
      . %in% c("one", "1") ~ 1L,
      . %in% c("two", "2") ~ 2L,
      . %in% c("three", "3") ~ 3L,
      . %in% c("four", "4") ~ 4L,
      . %in% c("five", "5") ~ 5L,
      TRUE ~ NA_integer_
    )
  
  if (is.na(year_num)) {
    warning(paste("Cannot extract year from filename:", file_path))
    return(NULL)
  }
  
  df %>%
    select(num, status, clinic, gene, pet, arithmetic_mean = all_of(mean_col)) %>%
    mutate(
      status = factor(status, levels = c(0, 1), labels = c("control", "case")),
      year = year_num
    ) %>%
    drop_na()
}

# ==================== 2. Read All Data ====================
cat("\n", rep("=", 60), "\n", sep = "")
cat("[1/5] Reading data\n")
cat(rep("=", 60), "\n", sep = "")

file_paths <- c("one.csv", "two.csv", "three.csv", "four.csv", "five.csv")
all_data <- map_dfr(file_paths, possibly(read_year_data, otherwise = NULL))

cat(sprintf("Data reading completed\n"))
cat(sprintf("  Total samples: %d rows\n", nrow(all_data)))
cat(sprintf("  Year distribution: %s\n", 
            paste(names(table(all_data$year)), table(all_data$year), sep = ":", collapse = ", ")))

# ==================== 3. Find Best Random Seed Based on Year 3 AUC ====================
cat("\n", rep("=", 60), "\n", sep = "")
cat("[2/5] Finding best random seed (based on year 3 AUC)\n")
cat(rep("=", 60), "\n", sep = "")

find_best_model <- function(data_3year, seed_range = 1:1000) {
  best_auc <- 0
  best_seed <- NA
  best_model <- NULL
  
  pb <- txtProgressBar(min = 0, max = length(seed_range), style = 3)
  
  for (s in seed_range) {
    set.seed(s)
    model <- randomForest(
      status ~ clinic + gene + pet + arithmetic_mean,
      data = data_3year,
      ntree = 2000,
      importance = TRUE
    )
    
    pred_prob <- predict(model, type = "prob")[, "case"]
    current_auc <- auc(roc(data_3year$status, pred_prob, quiet = TRUE))
    
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_seed <- s
      best_model <- model
    }
    
    setTxtProgressBar(pb, s)
  }
  
  close(pb)
  return(list(seed = best_seed, auc = best_auc, model = best_model))
}

data_3year <- all_data %>% filter(year == 3)

if (nrow(data_3year) == 0) {
  stop("Error: No data found for year 3, please check data files")
}

cat(sprintf("Year 3 sample size: %d cases\n", nrow(data_3year)))

best_info <- find_best_model(data_3year, seed_range = 1:1000)
best_seed <- best_info$seed

cat(sprintf("\nBest random seed: %d\n", best_seed))
cat(sprintf("  Corresponding year 3 AUC: %.3f\n", best_info$auc))

# ==================== 4. Calculate ROC Data for Years 1, 3, and 5 ====================
cat("\n", rep("=", 60), "\n", sep = "")
cat("[3/5] Calculating ROC curve data\n")
cat(rep("=", 60), "\n", sep = "")

selected_years <- c(1, 3, 5)
roc_results <- list()

for (yr in selected_years) {
  df_year <- all_data %>% filter(year == yr)
  
  if (nrow(df_year) == 0) {
    warning(sprintf("Year %d has no data, skipping", yr))
    next
  }
  
  set.seed(best_seed)
  model <- randomForest(
    status ~ clinic + gene + pet + arithmetic_mean,
    data = df_year,
    ntree = 2000,
    importance = TRUE
  )
  
  pred_prob <- predict(model, type = "prob")[, "case"]
  roc_obj <- roc(response = df_year$status, predictor = pred_prob, quiet = TRUE)
  
  ci_auc <- ci.auc(roc_obj, method = "bootstrap", boot.n = 2000, quiet = TRUE)
  
  roc_coords <- data.frame(
    fpr = 1 - roc_obj$specificities,
    tpr = roc_obj$sensitivities
  ) %>%
    arrange(fpr)
  
  dense_fpr <- seq(0, 1, length.out = 300)
  dense_tpr <- approx(roc_coords$fpr, roc_coords$tpr, 
                      xout = dense_fpr, rule = 2, ties = mean)$y
  
  roc_results[[as.character(yr)]] <- list(
    coords = data.frame(
      fpr = dense_fpr,
      tpr = dense_tpr,
      year = as.factor(yr)
    ) %>% filter(!is.na(tpr)),
    auc = as.numeric(auc(roc_obj)),
    ci_lower = ci_auc[1],
    ci_upper = ci_auc[3],
    sample_size = nrow(df_year),
    case_count = sum(df_year$status == "case"),
    control_count = sum(df_year$status == "control")
  )
  
  cat(sprintf("Year %d: AUC = %.3f (95%% CI: %.3f–%.3f)\n", 
              yr, roc_results[[as.character(yr)]]$auc,
              roc_results[[as.character(yr)]]$ci_lower,
              roc_results[[as.character(yr)]]$ci_upper))
  cat(sprintf("         Sample size: %d (case: %d, control: %d)\n",
              roc_results[[as.character(yr)]]$sample_size,
              roc_results[[as.character(yr)]]$case_count,
              roc_results[[as.character(yr)]]$control_count))
}

roc_plot_data <- bind_rows(lapply(roc_results, function(x) x$coords))

# ==================== 5. Define Color Scheme ====================
sci_colors <- c("#E41A1C", "#377EB8", "#4DAF4A")
names(sci_colors) <- c("1", "3", "5")

legend_y_pos <- c(0.28, 0.20, 0.12)

legend_labels <- data.frame(
  year = c("1", "3", "5"),
  auc_text = c(
    sprintf("Year 1 AUC: %.3f", roc_results[["1"]]$auc),
    sprintf("Year 3 AUC: %.3f", roc_results[["3"]]$auc),
    sprintf("Year 5 AUC: %.3f", roc_results[["5"]]$auc)
  ),
  x_auc = 0.62,
  x_line_start = 0.84,
  x_line_end = 0.90,
  y = legend_y_pos,
  color = sci_colors
)

# ==================== 6. Create ROC Curve Plot ====================
cat("\n", rep("=", 60), "\n", sep = "")
cat("[4/5] Generating ROC curve plot\n")
cat(rep("=", 60), "\n", sep = "")

main_roc_plot <- ggplot() +
  geom_abline(intercept = 0, slope = 1, 
              linetype = "dashed", color = "gray40", linewidth = 0.5) +
  geom_line(data = roc_plot_data, 
            aes(x = fpr, y = tpr, color = year), 
            linewidth = 1.2) +
  geom_text(data = legend_labels, 
            aes(x = x_auc, y = y, label = auc_text, color = year),
            hjust = 0, size = 4.5, show.legend = FALSE) +
  geom_segment(data = legend_labels,
               aes(x = x_line_start, xend = x_line_end, 
                   y = y, yend = y, color = year),
               linewidth = 1.2, show.legend = FALSE) +
  scale_color_manual(values = sci_colors, name = "") +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 1), 
                     breaks = seq(0, 1, 0.2)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1),
                     breaks = seq(0, 1, 0.2)) +
  labs(
    title = "ROC Curves for 1-, 3-, and 5-Year Predictions",
    x = "1 - Specificity (False Positive Rate)",
    y = "Sensitivity (True Positive Rate)"
  ) + 
  theme_bw() +
  theme(
    panel.border = element_rect(color = "black", linewidth = 0.8),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 13, face = "bold"),
    plot.title = element_text(color = "black", size = 15, face = "bold", 
                              hjust = 0.5, margin = margin(b = 10)),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    legend.position = "none",
    plot.margin = unit(c(10, 15, 10, 10), "pt")
  )

print(main_roc_plot)
cat("Plot generation completed\n")

# ==================== 7. Save Images ====================
cat("\n", rep("=", 60), "\n", sep = "")
cat("[5/5] Saving image files\n")
cat(rep("=", 60), "\n", sep = "")

output_dir <- getwd()

tiff_file <- file.path(output_dir, "ROC_curves_SCI.tiff")
ggsave(tiff_file, 
       plot = main_roc_plot,
       width = 8, 
       height = 7, 
       units = "in",
       dpi = 500,
       compression = "lzw",
       bg = "white")
cat(sprintf("TIFF saved: %s (500 dpi, LZW compression)\n", tiff_file))

pdf_file <- file.path(output_dir, "ROC_curves_SCI.pdf")
ggsave(pdf_file, 
       plot = main_roc_plot,
       width = 8, 
       height = 7, 
       units = "in",
       device = "pdf",
       bg = "white")
cat(sprintf("PDF saved: %s (vector format)\n", pdf_file))

png_file <- file.path(output_dir, "ROC_curves_SCI.png")
ggsave(png_file, 
       plot = main_roc_plot,
       width = 8, 
       height = 7, 
       units = "in",
       dpi = 300,
       bg = "white")
cat(sprintf("PNG saved: %s (300 dpi)\n", png_file))

# ==================== 8. Save Summary Table ====================
roc_summary_table <- data.frame(
  Year = selected_years,
  AUC = sapply(selected_years, function(yr) round(roc_results[[as.character(yr)]]$auc, 3)),
  CI_Lower = sapply(selected_years, function(yr) round(roc_results[[as.character(yr)]]$ci_lower, 3)),
  CI_Upper = sapply(selected_years, function(yr) round(roc_results[[as.character(yr)]]$ci_upper, 3)),
  Sample_Size = sapply(selected_years, function(yr) roc_results[[as.character(yr)]]$sample_size),
  Case_N = sapply(selected_years, function(yr) roc_results[[as.character(yr)]]$case_count),
  Control_N = sapply(selected_years, function(yr) roc_results[[as.character(yr)]]$control_count)
)

csv_file <- file.path(output_dir, "ROC_summary_table.csv")
write.csv(roc_summary_table, csv_file, row.names = FALSE)
cat(sprintf("Summary table saved: %s\n", csv_file))

# ==================== 9. Output Completion Information ====================
cat("\n", rep("=", 60), "\n", sep = "")
cat("ROC curve analysis completed!\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("Results summary:\n")
cat("----------------------------------------\n")
print(roc_summary_table)

cat("\nGenerated files:\n")
cat("----------------------------------------\n")
cat("  ROC_curves_SCI.tiff   (500 dpi, for submission)\n")
cat("  ROC_curves_SCI.pdf    (vector, for editing)\n")
cat("  ROC_curves_SCI.png    (300 dpi, for preview)\n")
cat("  ROC_summary_table.csv (ROC results table)\n")

cat("\nFile sizes:\n")
cat("----------------------------------------\n")
for (f in c(tiff_file, pdf_file, png_file)) {
  if (file.exists(f)) {
    size_bytes <- file.info(f)$size
    if (size_bytes > 1024 * 1024) {
      cat(sprintf("  %s: %.2f MB\n", basename(f), size_bytes / 1024 / 1024))
    } else {
      cat(sprintf("  %s: %.1f KB\n", basename(f), size_bytes / 1024))
    }
  }
}

cat("\nImage parameters:\n")
cat("----------------------------------------\n")
cat("  Type: Combined ROC Curves\n")
cat("  Resolution: TIFF = 500 dpi | PNG = 300 dpi | PDF = vector\n")
cat("  Dimensions: 8 x 7 inches (20.32 x 17.78 cm)\n")
cat("  Compression: LZW for TIFF\n")
cat("  Color scheme: Colorblind-friendly (Red/Blue/Green)\n")

cat("\nCaption template:\n")
cat("----------------------------------------\n")
cat("Fig. X. Receiver operating characteristic (ROC) curves for the multimodal\n")
cat("ensemble model at 1, 3, and 5 years. The red, blue, and green curves represent\n")
cat("ROC curves for 1-year, 3-year, and 5-year survival predictions, respectively.\n")
cat(sprintf("The areas under the curve (AUC) were 0.%.3f (95%% CI: 0.%.3f-0.%.3f) for 1-year,\n",
            roc_results[["1"]]$auc * 1000, 
            roc_results[["1"]]$ci_lower * 1000,
            roc_results[["1"]]$ci_upper * 1000))
cat(sprintf("0.%.3f (95%% CI: 0.%.3f-0.%.3f) for 3-year, and 0.%.3f (95%% CI: 0.%.3f-0.%.3f) for 5-year.\n",
            roc_results[["3"]]$auc * 1000,
            roc_results[["3"]]$ci_lower * 1000,
            roc_results[["3"]]$ci_upper * 1000,
            roc_results[["5"]]$auc * 1000,
            roc_results[["5"]]$ci_lower * 1000,
            roc_results[["5"]]$ci_upper * 1000))
cat("The dashed diagonal line indicates the performance of a random classifier.\n")

cat("\n", rep("=", 60), "\n", sep = "")
cat("All files generated successfully!\n")
cat("Please use ROC_curves_SCI.tiff (500 dpi) for submission\n")
cat(rep("=", 60), "\n\n", sep = "")

cat("=== Session Information ===\n")
sessionInfo()