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
plot_performance <- function(data, aggregation_column, orders) {
    
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
    
    natura_plots <- create_metric_plots(aggregated_natura, aggregation_column, orders$natura)
    corine_plots <- create_metric_plots(aggregated_corine, aggregation_column, orders$corine)
    difference_plots <- create_difference_plots(aggregated_difference, aggregation_column, orders$difference)
    
    grid.arrange(natura_plots, 
                 corine_plots,
                 difference_plots,
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


    






# SCRIPT STARTS

# LOAD IN MODELFITS AND DATA FOR COMPARISONS
modelfit_files <- list.files(dir_modelfits, pattern="*.RData", full.names=TRUE)
load(file.path(dir_data, "species_prevalences.RData"))
load(file.path(dir_data, "trait_data.RData"))
trait_data <- as.data.frame(trait_data)
rownames(trait_data) <- trait_data$Species
load(file = file.path(dir_data, "phylogeny_data.RData"))



# QUALITY OF MODELS


# CREATE PDF FOR PLOTTING THE FIT VALUES
pdf(file = file.path(dir_results, sprintf("model_comparisons.pdf")))


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
                                       aes(x = model, y = value, fill = model)) +
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

natura_species_comparison_df$average <- apply(natura_species_comparison_df, 1, mean)
rownames(natura_species_comparison_df) <- species_list
corine_species_comparison_df$average <- apply(corine_species_comparison_df, 1, mean)
rownames(corine_species_comparison_df) <- species_list
natura_corine_absolute_difference$average <- apply(natura_corine_absolute_difference, 1, mean)
rownames(natura_corine_absolute_difference) <- species_list



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


# CREATE FACTORS FOR PLOT ORDERING

metric_order <- c("average", colnames(natura_species_comparison_df[1:length(natura_species_comparison_df) - 1]))

long_df$metric <- factor(long_df$metric, levels = metric_order)

species_orders <- list()
species_orders$natura <- species[order(natura_species_comparison_df$average, decreasing = TRUE)]
species_orders$corine <- species[order(corine_species_comparison_df$average, decreasing = TRUE)]
natura_wins <- natura_species_comparison_df$average > corine_species_comparison_df$average
species_orders$difference <- species[order(natura_wins, 
                                           natura_corine_absolute_difference$average,
                                           decreasing = c(TRUE, TRUE))]






# CREATE PLOTS
plot_performance(long_df, "species", species_orders)
plot_performance(long_df, "sample_prevalence", NULL)
plot_performance(long_df, "transect_prevalence", NULL)
plot_performance(long_df, "year_prevalence", NULL)
plot_performance(long_df, "feeding", NULL)
plot_performance(long_df, "migration", NULL)






# COMPARE MODEL PERFORMANCE USING LINEAR MODELS








dev.off()


