# SCRIPT FOR COMPARING MODELFITS


# FUNCTIONS | GENERAL
extract_thinning_value <- function(filename) {
    parts <- strsplit(filename, "_")[[1]] 
    thinning_value <- as.numeric(parts[length(parts) - 1])
    print(thinning_value)
    return(thinning_value)
}


# FUNCTIONS | MODEL QUALITY


create_modelfit_plot <- function(explanatory_power,
                                 predictive_power, 
                                 predictive_power_type, 
                                 type, 
                                 model_name, 
                                 thinning_value) {
    
    if (!type %in% colnames(explanatory_power)) {
        return("Type not yet calculated or does not apply for this model")
    }
    if (!is.null(explanatory_power[,type])) {
        plot(explanatory_power[,type],
             predictive_power[,type],
             xlim = c(-1,1),
             ylim = c(-1,1),
             xlab = "explanatory power",
             ylab = sprintf("predictive power %s", predictive_power_type),
             cex.main = 0.8,
             main = sprintf("%s\n thin = %s: %s. \nexp pw: [%.3f, %.3f], mean = %.3f, sd = %.3f \npred pw: [%.3f, %.3f], mean = %.3f, sd = %.3f\n",
                            model_name,
                            as.character(thinning_value),
                            type,
                            min(explanatory_power[,type]),
                            max(explanatory_power[,type]),
                            mean(explanatory_power[,type], na.rm = TRUE),
                            sd(explanatory_power[,type]),
                            min(predictive_power[,type]),
                            max(predictive_power[,type]),
                            mean(predictive_power[,type], na.rm = TRUE),
                            sd(predictive_power[,type])))
        abline(h = 0)
        abline(v = 0)
    }
    
}


plot_fits_per_species <- function(type, data) {
    values <- data[,type]
    names(values) <- data$Species
    ordered_values <- sort(values, decreasing = FALSE)
    barplot(ordered_values,
            horiz = TRUE,
            las = 1,                    
            main = type)
}







# FUNCTIONS | MODEL COMPARISON

# FUNCTION FOR PLOTTING MODEL PERFORMANCE METRICS
plot_metric <- function(data) {
    metric_name <- data$Metric[1]
    ggplot(data, 
           aes(x = Value, y = reorder(Model, Value))) +
        geom_bar(stat = "identity", fill = "dodgerblue") +
        geom_text(aes(label = round(Value, 3)), hjust = -0.2, size = 4) +
        labs(x = "Value", y = "Model", title = paste("Model Comparison -", metric_name)) +
        theme_minimal()
}



# FUNCTION FOR PLOTTING SCALED PERFORMANCE VALUES 
plot_performance <- function(data, 
                             aggregation_column, 
                             orders, 
                             which_types, 
                             which_natura, 
                             which_corine, 
                             which_difference,
                             colors,
                             plot_type) {
    
    aggregated_natura <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                   data = subset(data, type == "Natura"),
                                   FUN = mean)

    aggregated_corine <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                   data = subset(data, type == "Corine"),
                                   FUN = mean)
    
    aggregated_difference <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                       data = subset(data, type == "Difference"),
                                       FUN = mean)
    
    
    aggregated_difference$winners <- ifelse(aggregated_corine$value > aggregated_natura$value, "corine", "natura")
    aggregated_difference$value <- ifelse(aggregated_difference$winners == "corine", 
                                          -1 * aggregated_difference$value, 
                                          aggregated_difference$value)
    
    aggregated_natura <- subset(aggregated_natura, metric %in% which_natura)
    aggregated_corine <- subset(aggregated_corine, metric %in% which_corine)
    aggregated_difference <- subset(aggregated_difference, metric %in% which_difference)
    
    if (plot_type == "line") {
        natura_plots <- create_metric_plot_line(aggregated_natura, aggregation_column, orders$natura, colors, "Natura metric averages")
        corine_plots <- create_metric_plot_line(aggregated_corine, aggregation_column, orders$corine, colors, "Corine metric averages")
        difference_plots <- create_difference_plot_line(aggregated_difference, aggregation_column, orders$difference, colors, "Metric differences")
    } else {
        natura_plots <- create_metric_plots(aggregated_natura, aggregation_column, orders$natura)
        corine_plots <- create_metric_plots(aggregated_corine, aggregation_column, orders$corine)
        difference_plots <- create_difference_plots(aggregated_difference, aggregation_column, orders$difference)
    }

    
    plots <- list(natura_plots, corine_plots, difference_plots)
    
    
    grid.arrange(grobs = plots[which_types],
                 ncol = 1)
    
}

create_metric_plots <- function(data, aggregation_column, order) {
    
    if (!is.null(order)) {
        data[,aggregation_column] <- factor(data[,aggregation_column],
                                            levels = order)
    }
    
    plots <- ggplot(data, aes(x = .data[[aggregation_column]], y = value, fill = metric)) +
        geom_col(width = 0.6) +
        coord_flip() +
        facet_wrap(~ metric, nrow = 1) +
        labs(x = aggregation_column, y = "Value") +
        theme_minimal()
    
    return (plots)
} 


create_metric_plot_line <- function(data, aggregation_column, order, colors, title) {
    
    
    if (!is.null(order)) {
        data[,aggregation_column] <- factor(data[,aggregation_column],
                                            levels = order)
    }
    
    plot <- ggplot(data, aes(x = .data[[aggregation_column]], y = value, group = metric, color = metric)) +
        geom_line() +
        scale_colour_manual(values = colors) +
        theme_minimal() +
        labs(
            title = sprintf("%s for each %s value", title, aggregation_column),
            x = aggregation_column,
            y = "Average metric value"
        )
    
    return (plot)
    
}
    
    
create_difference_plots <- function(data, aggregation_column, order) {
    
    if (!is.null(order)) {
        data[,aggregation_column] <- factor(data[,aggregation_column],
                                            levels = order)
    }
    
    plots <- ggplot(data, aes(x = .data[[aggregation_column]], y = value, fill = winners)) +
        geom_col(width = 0.6) +
        coord_flip() +
        facet_wrap(~ metric, nrow = 1) +
        labs(x = aggregation_column, y = "Difference") +
        scale_fill_manual(values = c("natura" = "red", "corine" = "blue")) +
        theme_minimal() 
        
    
    
    return (plots)
}


create_difference_plot_line <- function(data, aggregation_column, order, colors, title) {
    
    
    if (!is.null(order)) {
        data[, aggregation_column] <- factor(data[, aggregation_column],
                                             levels = order)
    }
    
    # Heatmap with aggregation_column and metric as axes
    plot <- ggplot(data, aes(x = .data[[aggregation_column]], 
                             y = metric, 
                             fill = value)) +
        geom_tile(color = "white") +  # optional: white borders between tiles
        scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
        theme_minimal() +
        labs(
            title = sprintf("%s heatmap", title),
            x = aggregation_column,
            y = "Metric",
            fill = "Value"
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    return(plot)
    
    
    
}


    






# SCRIPT STARTS

# LOAD IN MODELFITS AND DATA FOR COMPARISONS
modelfit_files <- list.files(dir_modelfits, pattern="*.RData", full.names=TRUE)
load(file.path(dir_data, "species_prevalences.RData"))
load(file.path(dir_data, "trait_data.RData"))
trait_data <- as.data.frame(trait_data)
rownames(trait_data) <- trait_data$Species
load(file = file.path(dir_data, "phylogeny_data.RData"))



# QUALITY OF MODELS





modelfit_results <- list()
model_names <- c("Corine", "Natura")

for (modelfit_file_number in 1:length(modelfit_files)) {
    
    # GET MODEL INFORMATION
    load(modelfit_files[modelfit_file_number])
    model_name <- model_names[modelfit_file_number]
    thinning_value <- extract_thinning_value(model_name)
    modelfit_file <- file.path(dir_modelfits, sprintf("modelfit_%s.RData", model_name))
    
    
    # LOAD MODELFIT RESULTS
    modelfit_results[[model_name]] <- list(explanatory_power = explanatory_power,
                                           predictive_power_transect = predictive_power_transect,
                                           predictive_power_year = predictive_power_year,
                                           waic = waic,
                                           waic_by_column = waic_by_column)

    
}

number_of_rows <- nrow(modelfit_results$Natura$explanatory_power)
number_of_categories <- 3
number_of_models <- 2
species_list <- rownames(modelfit_results$Natura$explanatory_power)

results_dataframe_base <- data.frame(model = c(rep("Natura", number_of_rows * number_of_categories),
                                               rep("Corine", number_of_rows * number_of_categories)),
                                     metric = rep(c(rep("Explanatory_power", number_of_rows),
                                                    rep("Predictive_power_transect", number_of_rows),
                                                    rep("Predictive_power_year", number_of_rows)),
                                                  number_of_models))
tjurr2_results <- results_dataframe_base
tjurr2_results$value <- c(modelfit_results$Natura$explanatory_power$TjurR2,
                          modelfit_results$Natura$predictive_power_transect$TjurR2,
                          modelfit_results$Natura$predictive_power_year$TjurR2,
                          modelfit_results$Corine$explanatory_power$TjurR2,
                          modelfit_results$Corine$predictive_power_transect$TjurR2,
                          modelfit_results$Corine$predictive_power_year$TjurR2)
tjurr2_results$species <- rep(species_list, number_of_models * number_of_categories)

auc_results <- results_dataframe_base
auc_results$value <- c(modelfit_results$Natura$explanatory_power$AUC,
                          modelfit_results$Natura$predictive_power_transect$AUC,
                          modelfit_results$Natura$predictive_power_year$AUC,
                          modelfit_results$Corine$explanatory_power$AUC,
                          modelfit_results$Corine$predictive_power_transect$AUC,
                          modelfit_results$Corine$predictive_power_year$AUC)
auc_results$species <- rep(species_list, number_of_models * number_of_categories)


rmse_results <- results_dataframe_base
rmse_results$value <- c(modelfit_results$Natura$explanatory_power$RMSE,
                          modelfit_results$Natura$predictive_power_transect$RMSE,
                          modelfit_results$Natura$predictive_power_year$RMSE,
                          modelfit_results$Corine$explanatory_power$RMSE,
                          modelfit_results$Corine$predictive_power_transect$RMSE,
                          modelfit_results$Corine$predictive_power_year$RMSE)
rmse_results$species <- rep(species_list, number_of_models * number_of_categories)

waic_results <- data.frame(model = c(rep("Natura", number_of_rows),
                                     rep("Corine", number_of_rows)),
                           value = c(modelfit_results$Natura$waic_by_column,
                                     modelfit_results$Corine$waic_by_column))





# PLOT AVERAGE COMPARISONS

overall_tjurr2_comparison_plot <- ggplot(tjurr2_results, 
                                         aes(x = metric, y = value, fill = model)) +
                                    geom_boxplot(position = position_dodge()) +
                                    labs(x = "Metric", y = "Tjur R²") +
                                    scale_fill_manual(values = c("Natura" = "lightblue",
                                                                 "Corine" = "dodgerblue")) +
                                    theme_minimal() +
                                    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
                                    coord_flip()

overall_auc_comparison_plot <- ggplot(auc_results, 
                                         aes(x = metric, y = value, fill = model)) +
                                    geom_boxplot(position = position_dodge()) +
                                    labs(x = "Metric", y = "AUC") +
                                    scale_fill_manual(values = c("Natura" = "lightblue",
                                                                 "Corine" = "dodgerblue")) +
                                    theme_minimal() +
                                    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
                                    coord_flip() 


overall_rmse_comparison_plot <- ggplot(rmse_results, 
                                         aes(x = metric, y = value, fill = model)) +
                                    geom_boxplot(position = position_dodge()) +
                                    labs(x = "Metric", y = "RMSE") +
                                    scale_fill_manual(values = c("Natura" = "lightblue",
                                                                 "Corine" = "dodgerblue")) +
                                    theme_minimal() +
                                    coord_flip()

overall_waic_comparison_plot <- ggplot(waic_results, 
                                       aes(x = model, y s= value, fill = model)) +
                                    geom_boxplot() +
                                    labs(x = "Model", y = "WAIC") +
                                    scale_fill_manual(values = c("Natura" = "lightblue",
                                                                 "Corine" = "dodgerblue")) +
                                    theme_minimal() +
                                    coord_flip()





grid.arrange(overall_tjurr2_comparison_plot,
             overall_auc_comparison_plot,
             overall_rmse_comparison_plot,
             overall_waic_comparison_plot,
             heights = c(3, 3, 3, 1.7),
             ncol = 1)





# SPECIES COMPARISON


# FIRST SCALE THE METRICS
# For metrics where smaller is better = multiply by -1 to change the direction
# Scale all metrics to 0-1 for comparison

max_waic <- max(max(modelfit_results$Natura$waic_by_column), max(modelfit_results$Corine$waic_by_column))
max_exp_rmse <- max(max(modelfit_results$Natura$explanatory_power$RMSE), max(modelfit_results$Corine$explanatory_power$RMSE))
max_pred_transect_rmse <- max(max(modelfit_results$Natura$predictive_power_transect$RMSE), max(modelfit_results$Corine$predictive_power_transect$RMSE))
max_pred_year_rmse <- max(max(modelfit_results$Natura$predictive_power_year$RMSE), max(modelfit_results$Corine$predictive_power_year$RMSE))

natura_species_comparison_df <- data.frame(scaled_waic = (-1 * modelfit_results$Natura$waic_by_column) + max_waic,
                                           exp_pw_TjurR2 = modelfit_results$Natura$explanatory_power$TjurR2,
                                           exp_pw_AUC = modelfit_results$Natura$explanatory_power$AUC,
                                           scaled_exp_pw_RMSE = (-1 * modelfit_results$Natura$explanatory_power$RMSE) + max_exp_rmse,
                                           pred_pw_transect_TjurR2 = modelfit_results$Natura$predictive_power_transect$TjurR2,
                                           pred_pw_transect_AUC = modelfit_results$Natura$predictive_power_transect$AUC,
                                           scaled_pred_pw_transect_RMSE = (-1 * modelfit_results$Natura$predictive_power_transect$RMSE) + max_pred_transect_rmse,
                                           pred_pw_year_TjurR2 = modelfit_results$Natura$predictive_power_year$TjurR2,
                                           pred_pw_year_AUC = modelfit_results$Natura$predictive_power_year$AUC,
                                           scaled_pred_pw_year_RMSE = (-1 * modelfit_results$Natura$predictive_power_year$RMSE) + max_pred_year_rmse)

corine_species_comparison_df <- data.frame(scaled_waic = (-1 * modelfit_results$Corine$waic_by_column) + max_waic,
                                           exp_pw_TjurR2 = modelfit_results$Corine$explanatory_power$TjurR2,
                                           exp_pw_AUC = modelfit_results$Corine$explanatory_power$AUC,
                                           scaled_exp_pw_RMSE = (-1 * modelfit_results$Corine$explanatory_power$RMSE) + max_exp_rmse,
                                           pred_pw_transect_TjurR2 = modelfit_results$Corine$predictive_power_transect$TjurR2,
                                           pred_pw_transect_AUC = modelfit_results$Corine$predictive_power_transect$AUC,
                                           scaled_pred_pw_transect_RMSE = (-1 * modelfit_results$Corine$predictive_power_transect$RMSE) + max_pred_transect_rmse,
                                           pred_pw_year_TjurR2 = modelfit_results$Corine$predictive_power_year$TjurR2,
                                           pred_pw_year_AUC = modelfit_results$Corine$predictive_power_year$AUC,
                                           scaled_pred_pw_year_RMSE = (-1 * modelfit_results$Corine$predictive_power_year$RMSE) + max_pred_year_rmse)

natura_corine_absolute_difference <- abs(abs(natura_species_comparison_df) - abs(corine_species_comparison_df))

explanatory_power_columns <- c("exp_pw_TjurR2", "exp_pw_AUC", "scaled_exp_pw_RMSE")
predictive_power_transect_columns <- c("pred_pw_transect_TjurR2", "pred_pw_transect_AUC", "scaled_pred_pw_transect_RMSE")
predictive_power_year_columns <- c("pred_pw_year_TjurR2", "pred_pw_year_AUC", "scaled_pred_pw_year_RMSE")
rmse_columns <- c("scaled_exp_pw_RMSE", "scaled_pred_pw_transect_RMSE", "scaled_pred_pw_year_RMSE")
auc_columns <- c("exp_pw_AUC", "pred_pw_transect_AUC", "pred_pw_year_AUC")
tjurr2_columns <- c("exp_pw_TjurR2", "pred_pw_transect_TjurR2", "pred_pw_year_TjurR2")

natura_species_comparison_df$average <- apply(natura_species_comparison_df, 1, mean)
natura_species_comparison_df$average_explanatory_power <- apply(natura_species_comparison_df[explanatory_power_columns], 1, mean)
natura_species_comparison_df$average_predictive_power_transect <- apply(natura_species_comparison_df[predictive_power_transect_columns], 1, mean)
natura_species_comparison_df$average_predictive_power_year <- apply(natura_species_comparison_df[predictive_power_year_columns], 1, mean)
natura_species_comparison_df$average_rmse <- apply(natura_species_comparison_df[rmse_columns], 1, mean)
natura_species_comparison_df$average_auc <- apply(natura_species_comparison_df[auc_columns], 1, mean)
natura_species_comparison_df$average_tjurr2 <- apply(natura_species_comparison_df[tjurr2_columns], 1, mean)

corine_species_comparison_df$average <- apply(corine_species_comparison_df, 1, mean)
corine_species_comparison_df$average_explanatory_power <- apply(corine_species_comparison_df[explanatory_power_columns], 1, mean)
corine_species_comparison_df$average_predictive_power_transect <- apply(corine_species_comparison_df[predictive_power_transect_columns], 1, mean)
corine_species_comparison_df$average_predictive_power_year <- apply(corine_species_comparison_df[predictive_power_year_columns], 1, mean)
corine_species_comparison_df$average_rmse <- apply(corine_species_comparison_df[rmse_columns], 1, mean)
corine_species_comparison_df$average_auc <- apply(corine_species_comparison_df[auc_columns], 1, mean)
corine_species_comparison_df$average_tjurr2 <- apply(corine_species_comparison_df[tjurr2_columns], 1, mean)


natura_corine_absolute_difference$average <- apply(natura_corine_absolute_difference, 1, mean)
natura_corine_absolute_difference$average_explanatory_power <- apply(natura_corine_absolute_difference[explanatory_power_columns], 1, mean)
natura_corine_absolute_difference$average_predictive_power_transect <- apply(natura_corine_absolute_difference[predictive_power_transect_columns], 1, mean)
natura_corine_absolute_difference$average_predictive_power_year <- apply(natura_corine_absolute_difference[predictive_power_year_columns], 1, mean)
natura_corine_absolute_difference$average_rmse <- apply(natura_corine_absolute_difference[rmse_columns], 1, mean)
natura_corine_absolute_difference$average_auc <- apply(natura_corine_absolute_difference[auc_columns], 1, mean)
natura_corine_absolute_difference$average_tjurr2 <- apply(natura_corine_absolute_difference[tjurr2_columns], 1, mean)

rownames(natura_species_comparison_df) <- species_list
rownames(corine_species_comparison_df) <- species_list
rownames(natura_corine_absolute_difference) <- species_list

natura_wins <- natura_species_comparison_df$average > corine_species_comparison_df$average


natura_bigger <- 1 * (natura_species_comparison_df > corine_species_comparison_df)
wins_natura <- apply(natura_bigger, 2, sum) 
signed_difference <- natura_corine_absolute_difference
signed_difference[!natura_bigger] <- -1 * natura_corine_absolute_difference[!natura_bigger]




# TRANSFORM TO LONG FORMAT
types <- c("Natura", "Corine", "Difference")
metrics <- colnames(natura_species_comparison_df)
species <- rownames(natura_species_comparison_df)
number_of_types <- length(types)
number_of_metrics <- length(metrics)
number_of_species <- length(species)
long_df <- data.frame(type = rep(types, each = number_of_metrics * number_of_species),
                      metric = rep(metrics, each = number_of_species, times = number_of_types),
                      species = rep(species, times = number_of_types * number_of_metrics),
                      value = c(as.vector(as.matrix(natura_species_comparison_df)),
                                as.vector(as.matrix(corine_species_comparison_df)),
                                as.vector(as.matrix(natura_corine_absolute_difference))))


# ADD TRAITS TO DF
long_df$sample_prevalence <- species_prevalences[long_df$species, c("samples")]
long_df$transect_prevalence <- species_prevalences[long_df$species, c("transects")]
long_df$year_prevalence <- species_prevalences[long_df$species, c("years")]
long_df$feeding <- trait_data[long_df$species, c("Feeding")]
long_df$mass <- trait_data[long_df$species, c("Mass")]
long_df$logmass <- trait_data[long_df$species, c("LogMass")]
long_df$migration <- trait_data[long_df$species, c("Mig")]

# Add transect prevalence group with equal prevalence range in each group
min_val <- min(long_df$transect_prevalence, na.rm = TRUE)
max_val <- max(long_df$transect_prevalence, na.rm = TRUE)
long_df$prevalence_group_equal_interval <- cut(
    long_df$transect_prevalence,
    breaks = seq(min_val, max_val, length.out = 6), # 5 groups
    include.lowest = TRUE,
    labels = FALSE
)

# Add transect prevalence group with equal number of species in each group
long_df$prevalence_group_equal_species <- cut(
    long_df$transect_prevalence,
    breaks = quantile(long_df$transect_prevalence, probs = seq(0, 1, 0.2), na.rm = TRUE),
    include.lowest = TRUE,
    labels = FALSE
)





# CREATE FACTORS FOR PLOT ORDERING

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

long_df$metric <- factor(long_df$metric, levels = metric_order)

species_orders <- list()
species_orders$natura <- species[order(natura_species_comparison_df$average, decreasing = TRUE)]
species_orders$corine <- species[order(corine_species_comparison_df$average, decreasing = TRUE)]
species_orders$difference <- species[order(natura_wins, 
                                           natura_corine_absolute_difference$average,
                                           decreasing = c(TRUE, TRUE))]






# CREATE SPECIES VS FIT PLOTS
all_true <- rep(TRUE, times = length(metrics))
all_false <- rep(FALSE, times = length(metrics))
only_averages <- all_false
only_averages[1:7] <- TRUE
no_averages <- all_true
no_averages[1:7] <- FALSE



plot_performance(long_df, "species", species_orders, c(FALSE, FALSE, TRUE), NULL, NULL, metric_order[only_averages])




max_absolute_difference <- max(natura_corine_absolute_difference)
color_palette_breaks <- seq(-max_absolute_difference, 
                            max_absolute_difference,
                            length.out = 11)
pheatmap(natura_bigger[,metric_order[8:length(metric_order)]],
         cluster_cols = FALSE)
pheatmap(signed_difference[,metric_order[8:length(metric_order)]],
         color = colorRampPalette(c("blue", "white", "red"))(10),
         breaks = color_palette_breaks,
         cluster_cols = FALSE)










# CREATE FIT VS PREVALENCE PLOTS
    



colors <- c("steelblue1", "steelblue2", "steelblue3", 
            "yellow1", "yellow2", "yellow3", 
            "green1", "darkgreen", "forestgreen", 
            "red")



plot_performance(long_df, 
                 "transect_prevalence", 
                 NULL, 
                 c(TRUE, TRUE, TRUE), 
                 metric_order[no_averages], 
                 metric_order[no_averages],
                 metric_order[no_averages], 
                 colors,
                 "line")
plot_performance(long_df, 
                 "feeding", 
                 NULL, 
                 c(TRUE, TRUE, TRUE), 
                 metric_order[no_averages], 
                 metric_order[no_averages],
                 metric_order[no_averages], 
                 colors,
                 "line")


plot_performance(long_df, "sample_prevalence", NULL, c(TRUE, FALSE, TRUE), c("scaled_waic", "exp_pw_AUC"), NULL,  c("average", "scaled_waic", "exp_pw_AUC"))

plot_performance(long_df, "year_prevalence", NULL, c(TRUE, FALSE, TRUE), c("scaled_waic", "exp_pw_TjurR2"), NULL,  c("average", "scaled_waic", "exp_pw_TjurR2"))

plot_performance(long_df, "prevalence_group_equal_interval", NULL, c(TRUE, TRUE, TRUE), metric_order[only_averages], metric_order[only_averages], metric_order[only_averages])
plot_performance(long_df, "prevalence_group_equal_interval", NULL, c(TRUE, TRUE, TRUE), metric_order[no_averages], metric_order[no_averages], metric_order[no_averages])
plot_performance(long_df, "prevalence_group_equal_species", NULL, c(TRUE, TRUE, TRUE), metric_order[only_averages], metric_order[only_averages], metric_order[only_averages])
plot_performance(long_df, "prevalence_group_equal_species", NULL, c(TRUE, TRUE, TRUE), metric_order[no_averages], metric_order[no_averages], metric_order[no_averages])



# How many species have each prevalence value?
barplot(table(species_prevalences[unique(long_df$species), c("samples")]))
barplot(table(species_prevalences[unique(long_df$species), c("transects")]))
barplot(table(species_prevalences[unique(long_df$species), c("years")]))






# CREATE FIT VS TRAIT PLOTS

plot_performance(long_df, "feeding", NULL)
plot_performance(long_df, "migration", NULL)






# COMPARE MODEL PERFORMANCE USING LINEAR MODELS

# Response variable: fit metrics
# Fixed effect: model
# Random effect: species

# First transform data to wide format for t-test
metrics_dataframes <- list()
for (metric_name in metrics) {
    # Transform data to wide format
    metric_values_natura <- long_df[long_df$metric == metric_name & long_df$type == "Natura", c("value")]
    metric_values_corine <- long_df[long_df$metric == metric_name & long_df$type == "Corine", c("value")]
    species <- long_df[long_df$metric == metric_name & long_df$type == "Natura", c("species")]
    metric_dataframe <- data.frame(values_natura = metric_values_natura,
                                   values_corine = metric_values_corine)
    rownames(metric_dataframe) <- species
    metrics_dataframes[[metric_name]] <- metric_dataframe
}

# Check for skew in the differences
for (metric_name in names(metrics_dataframes)) {
    data <- metrics_dataframes[[metric_name]]
    difference_between_corine_and_natura <- data$values_natura - data$values_corine
    hist(difference_between_corine_and_natura,
         main = metric_name)
}


# Perform t-test
test_results <- list()
for (metric_name in names(metrics_dataframes)) {
    data <- metrics_dataframes[[metric_name]]
    result <- t.test(data$values_natura, data$values_corine, paired = TRUE)
    test_results[[metric_name]] <- result
}

test_results




# CORRELATIONS FOR METRICS


calculate_between_model_correlations <- function(data, aggregation_column) {
    between_model_correlation <- c()
    for (metric_name in metrics) {
        metric_data_natura <- long_df[long_df$type == "Natura" & long_df$metric == metric_name,]
        metric_data_corine <- long_df[long_df$type == "Corine" & long_df$metric == metric_name,]
        if (!is.null(aggregation_column)) {
            metric_data_natura <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                            data = metric_data_natura,
                                            FUN = mean)
            metric_data_corine <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                            data = metric_data_corine,
                                            FUN = mean)
        }
        correlation <- cor(metric_data_natura$value, metric_data_corine$value)
        between_model_correlation <- c(between_model_correlation, correlation)
    }
    names(between_model_correlation) <- metrics
    between_model_correlation_df <- data.frame(metric = names(between_model_correlation),
                                               value = between_model_correlation)
    print(ggplot(between_model_correlation_df, aes(x = reorder(metric, value), y = value)) +
        geom_col(fill = "dodgerblue") +
        coord_flip() +  
        labs(x = "Metric", y = "Value", title = "Metric correlations between models") +
        theme_minimal())
    
}

calculate_between_model_correlations(long_df, NULL)
calculate_between_model_correlations(long_df, "sample_prevalence")
calculate_between_model_correlations(long_df, "transect_prevalence")
calculate_between_model_correlations(long_df, "year_prevalence")




calculate_within_model_correlations <- function(data, aggregation_column, included_columns) {
 
    columns_for_reshaping <- c(aggregation_column, "metric", "value")
  
    aggregated_natura <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                 data = data[data$type == "Natura",],
                                 FUN = mean)
    
    aggregated_corine <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                   data = data[data$type == "Corine",],
                                   FUN = mean)

    natura_wide_df <- reshape(aggregated_natura,
                              idvar = aggregation_column,
                              timevar = "metric",
                              direction = "wide")
    colnames(natura_wide_df) <- sub("^value\\.", "", colnames(natura_wide_df))
    rownames(natura_wide_df) <- natura_wide_df[,aggregation_column]
    natura_wide_df[,aggregation_column] <- NULL
    natura_wide_df <- natura_wide_df[,included_columns]
    
    corine_wide_df <- reshape(aggregated_corine,
                              idvar = aggregation_column,
                              timevar = "metric",
                              direction = "wide")
    colnames(corine_wide_df) <- sub("^value\\.", "", colnames(corine_wide_df))
    rownames(corine_wide_df) <- corine_wide_df[,aggregation_column]
    corine_wide_df[,aggregation_column] <- NULL
    corine_wide_df <- corine_wide_df[,included_columns]
    
    
    print(pheatmap(cor(natura_wide_df),
             display_numbers = TRUE,
             main = "Metric correlations within Natura"))
    print(pheatmap(cor(corine_wide_df),
             display_numbers = TRUE,
             main = "Metric correlations within Corine"))
       
}

calculate_within_model_correlations(long_df, "species", metric_order[no_averages])
calculate_within_model_correlations(long_df, "sample_prevalence", metric_order[no_averages])
calculate_within_model_correlations(long_df, "transect_prevalence", metric_order[no_averages])
calculate_within_model_correlations(long_df, "year_prevalence", metric_order[no_averages])





# COMPARE PREDICTIONS BETWEEN CORINE AND NATURA















dev.off()


