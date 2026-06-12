# SCRIPT FOR COMPARING MODELFITS


# FUNCTIONS | GENERAL

extract_thinning_value <- function(filename) {
    parts <- strsplit(filename, "_")[[1]] 
    thinning_value <- as.numeric(parts[length(parts) - 1])
    print(thinning_value)
    return(thinning_value)
}


# FUNCTIONS | MODEL QUALITY

create_modelfit_dataframe <- function(base_df, category, modelfit_results) {
    modelfit_dataframe <- base_df
    modelfit_dataframe$value <- c(modelfit_results$Natura$explanatory_power[,category],
                                  modelfit_results$Natura$predictive_power_transect[,category],
                                  modelfit_results$Natura$predictive_power_year[,category],
                                  modelfit_results$Corine$explanatory_power[,category],
                                  modelfit_results$Corine$predictive_power_transect[,category],
                                  modelfit_results$Corine$predictive_power_year[,category])
    return(modelfit_dataframe)
}


create_modelfit_comparison_plot <- function(modelfit_results, category_title, average_performance = NULL) {
    plot <- ggplot(modelfit_results,
                   aes(x = metric, y = value, fill = model)) +
        geom_boxplot(position = position_dodge()) +
        labs(x = "Metric", y = category_title) +
        scale_fill_manual(values = c("Natura" = "lightblue",
                                     "Corine" = "dodgerblue")) +
        theme_minimal() +
        geom_hline(yintercept = average_performance, linetype = "dashed", color = "red") + 
        coord_flip()
    return(plot)
}
    


# SCRIPT STARTS





# MODEL QUALITY EVALUATION

# Load in modelfits
modelfit_files <- list.files(dir_modelfits, pattern="*.RData", full.names=TRUE)
modelfit_results <- list()
model_names <- c("Corine", "Natura")

for (modelfit_file_number in 1:length(modelfit_files)) {
    
    # Get model information
    load(modelfit_files[modelfit_file_number])
    model_name <- model_names[modelfit_file_number]
    thinning_value <- extract_thinning_value(model_name)
    modelfit_file <- file.path(dir_modelfits, sprintf("modelfit_%s.RData", model_name))
    
    
    # Load fit results for model
    modelfit_results[[model_name]] <- list(explanatory_power = explanatory_power,
                                           predictive_power_transect = predictive_power_transect,
                                           predictive_power_year = predictive_power_year,
                                           waic = waic,
                                           waic_by_column = waic_by_column)

    
}


# Format fit results as dataframe
categories <- c("TjurR2", "AUC", "RMSE")
number_of_rows <- nrow(modelfit_results$Natura$explanatory_power)
number_of_categories <- length(categories)
number_of_models <- length(modelfit_files)
species_list <- rownames(modelfit_results$Natura$explanatory_power)

results_dataframe_base <- data.frame(model = c(rep("Natura", number_of_rows * number_of_categories),
                                               rep("Corine", number_of_rows * number_of_categories)),
                                     metric = rep(c(rep("Explanatory_power", number_of_rows),
                                                    rep("Predictive_power_transect", number_of_rows),
                                                    rep("Predictive_power_year", number_of_rows)),
                                                  number_of_models),
                                     species = rep(species_list, number_of_models * number_of_categories))


tjurr2_results <- create_modelfit_dataframe(results_dataframe_base, categories[1], modelfit_results)
auc_results <- create_modelfit_dataframe(results_dataframe_base, categories[2], modelfit_results)
rmse_results <- create_modelfit_dataframe(results_dataframe_base, categories[3], modelfit_results)
waic_results <- data.frame(model = c(rep("Natura", number_of_rows),
                                     rep("Corine", number_of_rows)),
                           value = c(modelfit_results$Natura$waic_by_column,
                                     modelfit_results$Corine$waic_by_column),
                           metric = "WAIC")






# Plot modelfit comparisons
overall_tjurr2_comparison_plot <- create_modelfit_comparison_plot(tjurr2_results, "Tjur R²", 0)
overall_auc_comparison_plot <- create_modelfit_comparison_plot(auc_results, "AUC", 0.5)
overall_rmse_comparison_plot <- create_modelfit_comparison_plot(rmse_results, "RMSE")
overall_waic_comparison_plot <- create_modelfit_comparison_plot(waic_results, "WAIC")



#grid.arrange(overall_tjurr2_comparison_plot,
#             overall_auc_comparison_plot,
#             overall_rmse_comparison_plot,
#             overall_waic_comparison_plot,
#             heights = c(3, 3, 3, 1.7),
#             ncol = 1)


grid.arrange(overall_auc_comparison_plot, 
             overall_rmse_comparison_plot)

print(overall_waic_comparison_plot)


hist(waic_results[waic_results$model == "Natura",]$value)
hist(waic_results[waic_results$model == "Corine",]$value)
qqnorm(waic_results[waic_results$model == "Natura",]$value)
qqline(waic_results[waic_results$model == "Natura",]$value)
qqnorm(waic_results[waic_results$model == "Corine",]$value)
qqline(waic_results[waic_results$model == "Corine",]$value)

mean(waic_results[waic_results$model == "Natura",]$value)
median(waic_results[waic_results$model == "Natura",]$value)
mean(waic_results[waic_results$model == "Corine",]$value)
median(waic_results[waic_results$model == "Corine",]$value)

sd(waic_results[waic_results$model == "Natura",]$value)
sd(waic_results[waic_results$model == "Corine",]$value)

quantile(waic_results[waic_results$model == "Natura",]$value, 
         c(0.25, 0.75), 
         na.rm = TRUE)
quantile(waic_results[waic_results$model == "Corine",]$value, 
         c(0.25, 0.75), 
         na.rm = TRUE)


sd(waic_results[waic_results$model == "Natura",]$value) / mean(waic_results[waic_results$model == "Natura",]$value)
sd(waic_results[waic_results$model == "Corine",]$value) / mean(waic_results[waic_results$model == "Corine",]$value)

sum(waic_results[waic_results$model == "Natura",]$value < waic_results[waic_results$model == "Corine",]$value)
length(waic_results[waic_results$model == "Corine",]$value)






# UNSCALED PAIRED T-TEST
metric_order <- c("average", 
                  "average_rmse",
                  "average_tjurr2",
                  "average_auc",
                  "average_explanatory_power",
                  "average_predictive_power_transect",
                  "average_predictive_power_year",
                  "exp_pw_TjurR2",
                  "pred_pw_transect_TjurR2",
                  "pred_pw_year_TjurR2",
                  "exp_pw_AUC",
                  "pred_pw_transect_AUC",
                  "pred_pw_year_AUC",
                  "scaled_exp_pw_RMSE",
                  "scaled_pred_pw_transect_RMSE",
                  "scaled_pred_pw_year_RMSE",
                  "scaled_waic")
all_true <- rep(TRUE, times = length(metric_order))
all_false <- rep(FALSE, times = length(metric_order))
only_averages <- all_false
only_averages[1:7] <- TRUE
no_averages <- all_true
no_averages[1:7] <- FALSE



# COMPARE MODEL PERFORMANCE USING LINEAR MODELS

# Response variable: fit metrics
# Fixed effect: model
# Random effect: species

# First transform data to wide format for t-test
metric_dataframe <- data.frame()
for (metric_name in metric_order[no_averages]) {
    # Transform data to wide format
    metric_values_natura <- long_df[long_df$metric == metric_name & long_df$type == "Natura", c("value")]
    metric_values_corine <- long_df[long_df$metric == metric_name & long_df$type == "Corine", c("value")]
    species <- long_df[long_df$metric == metric_name & long_df$type == "Natura", c("species")]
    metric_dataframe <- rbind(metric_dataframe, data.frame(values_natura = metric_values_natura,
                                                           values_corine = metric_values_corine,
                                                           species = species,
                                                           metric = rep(metric_name, length(species))))
}

# Check for skew in the differences
difference_between_corine_and_natura <- metric_dataframe$values_natura - metric_dataframe$values_corine
hist(difference_between_corine_and_natura)


# Perform t-test
result_total <- t.test(metric_dataframe$values_natura, metric_dataframe$values_corine, paired = TRUE)
print(result_total)

# Test waic
result_waic <- t.test(metric_dataframe[metric_dataframe$metric == "scaled_waic", ]$values_natura, 
                       metric_dataframe[metric_dataframe$metric == "scaled_waic", ]$values_corine, 
                       paired = TRUE)
print(result_waic)







# MODEL COMPARISON


par(mfrow = c(1, 1))
plot(ecdf(waic_results[waic_results$model == "Natura",]$value - waic_results[waic_results$model == "Corine",]$value),
     do.points = FALSE,
     verticals = TRUE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of specieswise WAIC difference",
     xlab = "difference",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))






dev.off()


