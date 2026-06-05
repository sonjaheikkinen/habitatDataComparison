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

plot_difference_line <- function(data,
                                 aggregation_columns) {
    
    par(mfrow = c(length(aggregation_columns),2))
    
    for (column in aggregation_columns) {
        natura <- subset(data, type == "Natura")
        corine <- subset(data, type == "Corine")
        difference <- natura 
        difference$value <- natura$value - corine$value
        aggregated_difference <- aggregate(reformulate(c(column, "metric"), response = "value"),
                                           data = difference,
                                           FUN = function(x) {
                                               quantile(x, 0.5, na.rm = TRUE)})
        aggregated_difference_25 <- aggregate(reformulate(c(column, "metric"), response = "value"),
                                           data = difference,
                                           FUN = function(x) {
                                               quantile(x, 0.25, na.rm = TRUE)})
        aggregated_difference_75 <- aggregate(reformulate(c(column, "metric"), response = "value"),
                                           data = difference,
                                           FUN = function(x) {
                                               quantile(x, 0.75, na.rm = TRUE)})
        plot(natura[,column],
             natura$value,
             col = "red",
             xlab = "Prevalence",
             ylab = "WAIC",
             main = column)
        points(corine[,column],
               corine$value,
               col = "blue")
        plot(aggregated_difference[,column],
             aggregated_difference$value,
             type = "l",
             xlab = "Prevalence",
             ylab = "WAIC difference",
             main = column)
        lines(aggregated_difference_25[,column],
              aggregated_difference_25$value,
              col = "grey")
        lines(aggregated_difference_75[,column],
              aggregated_difference_75$value,
              col = "grey")
        abline(0, 0)
    }
    

    
}


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
        natura_plots <- create_metric_plot_line(aggregated_natura, aggregation_column, orders$natura, colors, "Natura")
        corine_plots <- create_metric_plot_line(aggregated_corine, aggregation_column, orders$corine, colors, "Corine")
        difference_plots <- create_difference_plot_heatmap(aggregated_difference, aggregation_column, orders$difference, colors)
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
            title = sprintf("Scaled %s metric averages", title),
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


analyze_difference <- function(data,
                               aggregation_column,
                               which_difference) {
    
        

        aggregated_natura <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                       data = subset(data, type == "Natura"),
                                       FUN = mean)
        aggregated_corine <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                       data = subset(data, type == "Corine"),
                                       FUN = mean)
        aggregated_difference <- aggregate(reformulate(c(aggregation_column, "metric"), response = "value"),
                                           data = subset(data, type == "Difference"),
                                           FUN = mean)
    
        aggregated_difference$winners <- ifelse(aggregated_corine$value > aggregated_natura$value, 
                                                "corine", 
                                                "natura")
        aggregated_difference$value <- ifelse(aggregated_difference$winners == "corine", 
                                              -1 * aggregated_difference$value, 
                                              aggregated_difference$value)
        
        aggregated_difference <- subset(aggregated_difference, metric %in% which_difference)
        
        difference_plots <- create_difference_plot_heatmap(aggregated_difference, aggregation_column)
        
        plots <- list(difference_plots)
        
        grid.arrange(grobs = plots,
                     ncol = 1)
        
        
        difference <- subset(data, type == "Difference")
        natura <- subset(data, type == "Natura")
        corine <- subset(data, type == "Corine")
        difference$winner <- ifelse(natura$value > corine$value,
                                           "natura",
                                           "corine")
        difference$value <- ifelse(difference$winner == "corine",
                                   -1 * difference$value,
                                   difference$value)
        difference$natura_winner <- as.numeric(difference$winner == "natura")
        difference <- subset(difference, metric %in% which_difference)
        
        natura_wins <- table(difference[[aggregation_column]], difference$value > 0)
        natura_wins <- data.frame(wins = natura_wins[,"TRUE"],
                                  losses = natura_wins[,"FALSE"],
                                  fraction = natura_wins[,"TRUE"] / rowSums(natura_wins))
        print(natura_wins)
        
        fit <- glm(natura_winner ~ get(aggregation_column), data = difference, family = binomial)
        summary(fit)
        
}



create_difference_plot_heatmap <- function(data, aggregation_column, order = NULL) {
    
    
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
            title = "Scaled metric differences heatmap",
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
load(file.path(dir_data, "env_data_natura.RData"))
load(file.path(dir_data, "env_data_corine.RData"))
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
                                       aes(x = model, y= value, fill = model)) +
                                    geom_boxplot() +
                                    labs(x = "Model", y = "WAIC") +
                                    scale_fill_manual(values = c("Natura" = "lightblue",
                                                                 "Corine" = "dodgerblue")) +
                                    theme_minimal() +
                                    coord_flip()





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






# First transform data to wide format for t-test
metric_dataframe <- data.frame(values_natura = waic_results[waic_results$model == "Natura",]$value,
                               values_corine = waic_results[waic_results$model == "Corine",]$value,
                               species = species_list,
                               metric = rep("WAIC", length(species_list)))
add_metric <- function(dataframe, metric_data, species_list,metric_name) {
    types <- c("Explanatory_power", "Predictive_power_transect", "Predictive_power_year")
    for (metric_type in types) {
        dataframe <- rbind(dataframe, data.frame(values_natura = metric_data[metric_data$model == "Natura" & metric_data$metric == metric_type,]$value,
                                                 values_corine = metric_data[metric_data$model == "Corine" & metric_data$metric == metric_type,]$value,
                                                 species = species_list,
                                                 metric = rep(sprintf("%s_%s", metric_name, metric_type), length(species_list))))
    }
    return (dataframe)
}
metric_dataframe <- add_metric(metric_dataframe, rmse_results, species_list, "RMSE")
metric_dataframe <- add_metric(metric_dataframe, tjurr2_results, species_list, "TjurR2")
metric_dataframe <- add_metric(metric_dataframe, auc_results, species_list, "AUC")


# Check for skew in the differences
difference_between_corine_and_natura <- metric_dataframe$values_natura - metric_dataframe$values_corine
hist(difference_between_corine_and_natura)


# Perform t-test
result_total <- t.test(metric_dataframe$values_natura, metric_dataframe$values_corine, paired = TRUE)
print(result_total)

for (metric in unique(metric_dataframe$metric)) {
    print(metric)
    print(t.test(metric_dataframe[metric_dataframe$metric == metric, ]$values_natura, 
                          metric_dataframe[metric_dataframe$metric == metric, ]$values_corine, 
                          paired = TRUE))
}



# UNSCALED CORRELATION
metric_dataframe_wide_natura <- data.frame(WAIC = metric_dataframe[metric_dataframe$metric == "WAIC",]$values_natura)
metric_dataframe_wide_corine <- data.frame(WAIC = metric_dataframe[metric_dataframe$metric == "WAIC",]$values_corine)
for (metric in unique(metric_dataframe$metric)) {
    metric_dataframe_wide_natura[,metric] <- metric_dataframe[metric_dataframe$metric == metric,]$values_natura
    metric_dataframe_wide_corine[,metric] <- metric_dataframe[metric_dataframe$metric == metric,]$values_corine
}

correlations <- data.frame(metric = unique(metric_dataframe$metric),
                           value = rep(0, length(unique(metric_dataframe$metric))))

for (metric in unique(metric_dataframe$metric)) {
    correlations[correlations$metric == metric,]$value <- cor(metric_dataframe_wide_natura[,metric],
                                                        metric_dataframe_wide_corine[,metric])
}
print(ggplot(correlations, aes(x = reorder(metric, value), y = value)) +
          geom_col(fill = "dodgerblue") +
          coord_flip() +  
          labs(x = "Metric", y = "Value", title = "Metric correlations between models") +
          theme_minimal())


pheatmap(cor(metric_dataframe_wide_corine),
         display_numbers = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)

pheatmap(cor(metric_dataframe_wide_natura),
         display_numbers = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)



# UNSCALED TRAIT COMPARISON
metric_dataframe_wide_difference <- metric_dataframe_wide_natura - metric_dataframe_wide_corine
rownames(metric_dataframe_wide_difference) <- species_list

for (metric in colnames(metric_dataframe_wide_difference)) {
    df_for_metric <- metric_dataframe_wide_difference[,metric,drop=FALSE]
    df_for_metric <- df_for_metric[order(df_for_metric[,metric]),,drop=FALSE]  
    mass <- trait_data[trait_data$Species == df_for_metric$species,]$Mass
    df_for_metric_long <- data.frame(species = rep(rownames(df_for_metric), 2),
                                     attribute = c(rep(metric, nrow(df_for_metric)),
                                                   rep("mass", nrow(df_for_metric))),
                                     value = c(df_for_metric[,metric],
                                               mass / max(mass)))
    df_for_metric_long$species <- factor(rownames(df_for_metric), 
                                       levels = rownames(df_for_metric))
    colors <- rainbow(length(unique(df_for_metric_long$attribute)))
    ggplot(df_for_metric_long, aes(x = species, y = value, group = attribute, color = attribute)) +
        geom_line() +
        scale_colour_manual(values = colors) +
        theme_minimal() +
        labs(title = sprintf("%s", metric),
            x = "Species",
            y = "Value"
        ) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
}

colors <- c("steelblue1", "steelblue2", "steelblue3", 
            "yellow1", "yellow2", "yellow3", 
            "green1", "darkgreen", "forestgreen", 
            "red")

metric_dataframe_long <- data.frame(type = c(rep("Natura", nrow(metric_dataframe)),
                                             rep("Corine", nrow(metric_dataframe))),
                                    value = c(metric_dataframe$values_natura, 
                                              metric_dataframe$values_corine),
                                    species = rep(metric_dataframe$species, 2),
                                    metric = rep(metric_dataframe$metric, 2))
metric_dataframe_long$mass <- trait_data[metric_dataframe_long$species,]$Mass
metric_dataframe_long$logmass <- trait_data[metric_dataframe_long$species,]$LogMass
metric_dataframe_long$feeding <- trait_data[metric_dataframe_long$species,]$Feeding
metric_dataframe_long$mig <- trait_data[metric_dataframe_long$species,]$Mig
metric_dataframe_long$transect_prevalence <- species_prevalences[metric_dataframe_long$species,]$transects
metric_dataframe_long$year_prevalence <- species_prevalences[metric_dataframe_long$species,]$years
metric_dataframe_long$sample_prevalence <- species_prevalences[metric_dataframe_long$species,]$samples

plot_data <- metric_dataframe_long[metric_dataframe_long$metric == "WAIC",]



plot(plot_data$transect_prevalence, plot_data$value,
     xlab = "Prevalence", 
     ylab = "WAIC",
     main = "Transect prevalence")
plot(plot_data$year_prevalence, plot_data$value,
     xlab = "Prevalence", 
     ylab = "WAIC",
     main = "Year prevalence")
plot(plot_data$sample_prevalence, plot_data$value,
     xlab = "Prevalence", 
     ylab = "WAIC",
     main = "Sample prevalence")


plot_difference_line(plot_data, c("transect_prevalence",
                                  "year_prevalence",
                                  "sample_prevalence"))

difference_data <- plot_data[plot_data$type == "Natura",]
difference_data$value <- difference_data$value - plot_data[plot_data$type == "Corine",]$value


par(mfrow = c(1, 3))
boxplot(value ~ feeding,
        data = plot_data[plot_data$type == "Natura",],
        main = "Natura WAIC")
boxplot(value ~ feeding,
        data = plot_data[plot_data$type == "Corine",],
        main = "Corine WAIC")
boxplot(value ~ feeding,
        data = difference_data,
        main = "Difference Natura - Corine")
abline(0, 0)



pheatmap(cor(metric_dataframe_wide_difference, cbind(transect_prevalence = species_prevalences[rownames(metric_dataframe_wide_difference), c("transects")], 
                                                     mass = trait_data[rownames(metric_dataframe_wide_difference),]$Mass)),
         display_numbers = TRUE,
         cluster_rows = FALSE)

                                

sapply(metric_dataframe_wide_natura, min, na.rm = TRUE)
sapply(metric_dataframe_wide_natura, max, na.rm = TRUE)



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

long_df$prevalence_group_equal_species <- cut(long_df$transect_prevalence,
                                              breaks = quantile(long_df$transect_prevalence, 
                                                                probs = seq(0, 1, 0.2), 
                                                                na.rm = TRUE),
                                              include.lowest = TRUE,
                                              labels = FALSE)


# Add mass groups with equal number of species in each group
min_val <- min(long_df$logmass, na.rm = TRUE)
max_val <- max(long_df$logmass, na.rm = TRUE)
long_df$logmass_group_equal_interval <- cut(long_df$logmass,
                                            breaks = seq(min_val, 
                                                         max_val, 
                                                         length.out = 6),
                                            include.lowest = TRUE,
                                            labels = FALSE)




# CREATE FACTORS FOR PLOT ORDERING


long_df$metric <- factor(long_df$metric, levels = metric_order)

species_orders <- list()
species_orders$natura <- species[order(natura_species_comparison_df$average, decreasing = TRUE)]
species_orders$corine <- species[order(corine_species_comparison_df$average, decreasing = TRUE)]
species_orders$difference <- species[order(natura_wins, 
                                           natura_corine_absolute_difference$average,
                                           decreasing = c(TRUE, TRUE))]






# CREATE SPECIES VS FIT PLOTS




max_absolute_difference <- max(natura_corine_absolute_difference)
color_palette_breaks <- seq(-max_absolute_difference, 
                            max_absolute_difference,
                            length.out = 11)
species_division <- pheatmap(natura_bigger[,metric_order[8:length(metric_order)]],
                            cluster_cols = FALSE)
species_clusters <- cluster_assignments <- cutree(species_division$tree_row, k = 2)
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
                 c(FALSE, FALSE, TRUE), 
                 metric_order[no_averages], 
                 metric_order[no_averages],
                 metric_order[no_averages], 
                 colors,
                 "line")
plot_performance(long_df, 
                 "feeding", 
                 NULL, 
                 c(FALSE, FALSE, TRUE), 
                 metric_order[no_averages], 
                 metric_order[no_averages],
                 metric_order[no_averages], 
                 colors,
                 "line")
plot_performance(long_df, 
                 "migration", 
                 NULL, 
                 c(TRUE, TRUE, TRUE), 
                 metric_order[no_averages], 
                 metric_order[no_averages],
                 metric_order[no_averages], 
                 colors,
                 "line")
plot_performance(long_df, 
                 "logmass", 
                 NULL, 
                 c(TRUE, TRUE, TRUE), 
                 metric_order[no_averages], 
                 metric_order[no_averages],
                 metric_order[no_averages], 
                 colors,
                 "line")
plot_performance(long_df, 
                 "logmass_group_equal_interval", 
                 NULL, 
                 c(TRUE, TRUE, TRUE), 
                 metric_order[no_averages], 
                 metric_order[no_averages],
                 metric_order[no_averages], 
                 colors,
                 "line")


long_df$feeding_category <- ifelse(long_df$feeding %in% c("M", "O"),
                                   "generalist",
                                   "specialist")


analyze_difference(long_df, 
                 "feeding_category", 
                 metric_order[no_averages])
analyze_difference(long_df, 
                   "mass", 
                   metric_order[no_averages])







# How many species have each prevalence value?
barplot(table(species_prevalences[unique(long_df$species), c("samples")]))
barplot(table(species_prevalences[unique(long_df$species), c("transects")]))
barplot(table(species_prevalences[unique(long_df$species), c("years")]))




# FIT VS HABITATS PLOTS

fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

corine_habitat_variables <- c("Havumetsät.kivennäismaalla",
                              "Sekametsät.kivennäismaalla",
                              "Sekametsät.turvemaalla",
                              "Lehtimetsät.kivennäismaalla",
                              "Havumetsät.kalliomaalla",
                              "CorinePatchDensity")
env_data_corine$CorinePatchDensity <- env_data_corine$PatchDensity
natura_habitat_variables <- c("Luonnonmetsät",
                              "Tunturikoivikot",
                              "Lehdot",
                              "Tulvametsät",
                              "NaturaPatchDensity")
env_data_natura$NaturaPatchDensity <- env_data_natura$PatchDensity
other_variables <- c("Temperature",
                     "Rainfall")

for (variable in c(corine_habitat_variables, natura_habitat_variables, other_variables)) {
    long_df[,variable] <- rep(0, nrow(long_df))
}

for (species in unique(long_df$species)) {
    load(fitted_models[[1]])
    rows_for_species <- fitted_model$Y[,species] == 1
    env_data_natura_species <- env_data_natura[rows_for_species,]
    env_data_corine_species <- env_data_corine[rows_for_species,]
    for (variable in natura_habitat_variables) {
        long_df[long_df$species == species, variable] <- mean(env_data_natura_species[,variable])
        
    }
    for (variable in corine_habitat_variables) {
        long_df[long_df$species == species, variable] <- mean(env_data_corine_species[,variable])
    }
    for (variable in other_variables) {
        long_df[long_df$species == species, variable] <- mean(env_data_natura_species[,variable])
    }
}


for (variable in c(corine_habitat_variables, natura_habitat_variables, other_variables)) {
    grouped_variable_name <- sprintf("grouped_%s", variable)
    long_df[,grouped_variable_name] <- cut(long_df[,variable],
                                           breaks = seq(min(long_df[,variable]), 
                                                        max(long_df[,variable]), 
                                                        length.out = 6),
                                           include.lowest = TRUE,
                                           labels = FALSE)
    plot_performance(long_df, 
                     grouped_variable_name, 
                     NULL, 
                     c(TRUE, TRUE, TRUE), 
                     metric_order[no_averages], 
                     metric_order[no_averages],
                     metric_order[no_averages], 
                     colors,
                     "line")
}




# CHECK SPECIES TRAIT CORRELATIONS

species_traits_and_ecological_traits <- c("sample_prevalence",
                                          "transect_prevalence",
                                          "year_prevalence",
                                          "feeding",
                                          "mass",
                                          "logmass",
                                          "migration",
                                          "Havumetsät.kivennäismaalla",
                                          "Sekametsät.kivennäismaalla",
                                          "Sekametsät.turvemaalla",
                                          "Lehtimetsät.kivennäismaalla",
                                          "Havumetsät.kalliomaalla",
                                          "NaturaPatchDensity",
                                          "CorinePatchDensity",
                                          "Luonnonmetsät",
                                          "Tunturikoivikot",
                                          "Lehdot",
                                          "Tulvametsät",
                                          "Temperature",
                                          "Rainfall")
species_expanded_traits <- unique(long_df[long_df$type == "Natura", 
                                          c(species_traits_and_ecological_traits, "species")])
rownames(species_expanded_traits) <- species_expanded_traits$species
species_expanded_traits$species <- NULL
                                 
     
plot(species_expanded_traits)
boxplot(mass ~ feeding, 
        data = species_expanded_traits,
        main = "Mass by Feeding Category",
        xlab = "Feeding Category",
        ylab = "Mass",
        col = "lightblue",
        border = "darkblue")



species_expanded_traits$cluster <- rep(0, nrow(species_expanded_traits))
for (species in rownames(species_expanded_traits)) {
    species_expanded_traits[species,]$cluster <- species_clusters[species]
}
species_expanded_traits$cluster <- ifelse(species_expanded_traits$cluster == 1,
                                  "corine",
                                  "natura")

numeric_cols <- sapply(species_expanded_traits, is.numeric)
numeric_names <- names(numeric_cols)[numeric_cols]

par(mfrow = c(3, 3))  # adjust grid depending on how many traits you have
for (trait in numeric_names) {
    print(boxplot(species_expanded_traits[,trait] ~ species_expanded_traits$cluster,
            main = paste("Trait:", trait),
            xlab = "Cluster",
            ylab = trait,
            col = c("lightblue", "lightgreen")))
}
par(mfrow = c(1, 1))  # reset layout


cat_cols <- sapply(species_expanded_traits, function(x) is.character(x) || is.factor(x))
cat_names <- names(cat_cols)[cat_cols]

par(mfrow = c(2, 2))  
for (trait in cat_names) {
    tab <- table(species_expanded_traits[,trait], species_expanded_traits$cluster)
    barplot(tab,
            beside = TRUE,
            main = paste("Trait:", trait),
            xlab = trait,
            ylab = "Count",
            legend.text = TRUE,
            args.legend = list(title = "Cluster", x = "topright"))
}
par(mfrow = c(1, 1))







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

# Test predictive power transect
transect_metrics <- c("scaled_pred_pw_transect_RMSE", "pred_pw_transect_AUC", "pred_pw_transect_TjurR2")
result_transect <- t.test(metric_dataframe[metric_dataframe$metric %in% transect_metrics, ]$values_natura, 
                          metric_dataframe[metric_dataframe$metric %in% transect_metrics, ]$values_corine, 
                          paired = TRUE)
print(result_transect)


# Test predictive power year
year_metrics <- c("scaled_pred_pw_year_RMSE", "pred_pw_year_AUC", "pred_pw_year_TjurR2")
result_year <- t.test(metric_dataframe[metric_dataframe$metric %in% year_metrics, ]$values_natura, 
                          metric_dataframe[metric_dataframe$metric %in% year_metrics, ]$values_corine, 
                          paired = TRUE)
print(result_year)





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





# COMPARE PREDICTIONS TO MODEL FIT

# 















dev.off()


