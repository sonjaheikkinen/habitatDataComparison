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



# FUNCTIONS | MODEL COMPARISON

plot_model_evaluation_ecdf <- function(data, metrics, title, limit, difference = FALSE) {
    
    plot(ecdf(data[data$metric == metrics[1],]$value),
         verticals = TRUE,
         do.points = FALSE,
         col = "blue",
         lwd = 3,
         yaxt = "n",
         main = title,
         xlab = "value",
         ylab = "cumulative probability",
         xlim = limit)
    axis(2, at = c(0, 0.25, 0.5, 0.75, 1), 
         labels = c(0, 0.25, 0.5, 0.75, 1),
         las = 2)
    if (length(metrics) > 1) {
        lines(ecdf(data[data$metric == metrics[2],]$value),
              verticals = TRUE,
              do.points = FALSE,
              col = "lightblue",
              lwd = 3)
    }
    abline(0.5, 0, col = alpha("black", 0.25))
    abline(0.25, 0, col = alpha("black", 0.25))
    abline(0.75, 0, col = alpha("black", 0.25))
    if (length(metrics) > 1) {
        legend("bottomright",
               legend = c("E (hab)", "E (lc)"),
               col = c("lightblue", "blue"),
               lty = 1)
    }
    if (difference) {
        abline(v = 0, col = alpha("black", 0.25))
    }
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










# MODEL COMPARISON


# Load overall_long_df
load(file.path(dir_results, "prediction_dataframe.RData"))
# Load unit_data
load(file.path(dir_results, "unit_data_combined.RData"))



# OVERALL MODEL COMPARISON




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



# Value correlations

par(mfrow = c(1, 3))


plot(overall_long_df[overall_long_df$metric == "natura_occurrence_prob",]$value,
     overall_long_df[overall_long_df$metric == "corine_occurrence_prob",]$value,
     xlab = "E (hab)",
     ylab = "E (lc)",
     main = "Expected occurrence probability E (hab. vs. lc)")

plot(overall_long_df[overall_long_df$metric == "natura_error",]$value,
     overall_long_df[overall_long_df$metric == "corine_error",]$value,
     xlab = "Error (hab)",
     ylab = "Error (lc)",
     main = "Mean absolute error MAE (hab vs. lc)")

plot(overall_long_df[overall_long_df$metric == "natura_uncertainty",]$value,
     overall_long_df[overall_long_df$metric == "corine_uncertainty",]$value,
     xlab = "Uncertainty (hab)",
     ylab = "Uncertainty (lc)",
     main = "Uncertainty UC (hab vs. lc)")





cor(overall_long_df[overall_long_df$metric == "uncertainty_difference",]$value,
    overall_long_df[overall_long_df$metric == "error_difference",]$value)





# Value distributions





par(mfrow = c(2, 3))


plot_model_evaluation_ecdf(overall_long_df, 
                           c("corine_occurrence_prob", "natura_occurrence_prob"), 
                           "CDF of expected occurrence probability E",
                           c(0, 1))

plot_model_evaluation_ecdf(overall_long_df,
                           c("corine_error", "natura_error"), 
                           "CDF of mean absolute error MAE",
                           c(0, 1))

plot_model_evaluation_ecdf(overall_long_df, 
                           c("corine_uncertainty", "natura_uncertainty"),
                           "CDF of uncertainty UC",
                           c(0, 1))


values <- overall_long_df[overall_long_df$metric == "occurrence_prob_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(overall_long_df,
                           "occurrence_prob_difference",
                           "CDF of E difference",
                           c(-limit, limit),
                           difference = TRUE)


values <- overall_long_df[overall_long_df$metric == "error_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(overall_long_df,
                           "error_difference",
                           "CDF of MAE difference",
                           c(-limit, limit),
                           difference = TRUE)

values <- overall_long_df[overall_long_df$metric == "uncertainty_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(overall_long_df,
                           "uncertainty_difference",
                           "CDF of UC difference",
                           c(-limit, limit),
                           difference = TRUE)



# TRANSECTWISE ANALYSIS



transect_long_df <- aggregate(overall_long_df,
                              cbind(value, x, y) ~ transect + metric,
                              FUN = median)



# Names for environmental variables
natura_habitat_variables <- c("Luonnonmetsät",
                              "Tunturikoivikot",
                              "Lehdot",
                              "Tulvametsät",
                              "NaturaPatchDensity")
corine_habitat_variables <- c("Havumetsät.kivennäismaalla",
                              "Sekametsät.kivennäismaalla",
                              "Sekametsät.turvemaalla",
                              "Lehtimetsät.kivennäismaalla",
                              "Havumetsät.kalliomaalla",
                              "CorinePatchDensity")
other_variables <- c("Temperature", "Rainfall")



par(mfrow = c(1, 3))


plot_model_evaluation_ecdf(transect_long_df, 
                           c("corine_occurrence_prob", "natura_occurrence_prob"), 
                           "CDF of expected occurrence probability E",
                           c(0, 1))

plot_model_evaluation_ecdf(transect_long_df,
                           c("corine_error", "natura_error"), 
                           "CDF of mean absolute error MAE",
                           c(0, 1))

plot_model_evaluation_ecdf(transect_long_df, 
                           c("corine_uncertainty", "natura_uncertainty"),
                           "CDF of uncertainty UC",
                           c(0, 1))


values <- transect_long_df[transect_long_df$metric == "occurrence_prob_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(transect_long_df,
                           "occurrence_prob_difference",
                           "CDF of expected occurrence probablity difference \nfor transect medians",
                           c(-limit, limit),
                           difference = TRUE)


values <- transect_long_df[transect_long_df$metric == "error_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(transect_long_df,
                           "error_difference",
                           "CDF of mean absolute error difference \nfor transect medians",
                           c(-limit, limit),
                           difference = TRUE)

values <- transect_long_df[transect_long_df$metric == "uncertainty_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(transect_long_df,
                           "uncertainty_difference",
                           "CDF of uncertainty difference \nfor transect medians",
                           c(-limit, limit),
                           difference = TRUE)









# Transect prediction decision tree






transect_habitat_variables <- unit_data[,c(natura_habitat_variables, 
                                           corine_habitat_variables, 
                                           other_variables)]
transect_habitat_variables$transect_name <- unit_data$Transect
other_columns <- setdiff(colnames(transect_habitat_variables), "transect_name")
transect_habitat_variables <- aggregate(transect_habitat_variables,
                                        cbind(Luonnonmetsät,
                                              Tunturikoivikot,
                                              Lehdot,
                                              Tulvametsät,
                                              NaturaPatchDensity,
                                              Havumetsät.kivennäismaalla,
                                              Sekametsät.kivennäismaalla,
                                              Sekametsät.turvemaalla,
                                              Lehtimetsät.kivennäismaalla,
                                              Havumetsät.kalliomaalla,
                                              CorinePatchDensity,
                                              Temperature,
                                              Rainfall) ~ transect_name,
                                        FUN = mean)
rownames(transect_habitat_variables) <- transect_habitat_variables$transect_name
transect_habitat_variables$transect_name <- NULL


long_data_error <- transect_long_df[transect_long_df$metric == "error_difference",]
decision_tree_data_error <- transect_habitat_variables
decision_tree_data_error$value <- rep(0, nrow(decision_tree_data_error))
for (transect in rownames(decision_tree_data_error)) {
    decision_tree_data_error[transect,]$value <- long_data_error[long_data_error$transect == transect,"value"]
}
error_decision_tree <- rpart(value ~., data = decision_tree_data_error)
rpart.plot(error_decision_tree)
plotcp(error_decision_tree)


long_data_error <- transect_long_df[transect_long_df$metric == "uncertainty_difference",]
decision_tree_data_error <- transect_habitat_variables
decision_tree_data_error$value <- rep(0, nrow(decision_tree_data_error))
for (transect in rownames(decision_tree_data_error)) {
    decision_tree_data_error[transect,]$value <- long_data_error[long_data_error$transect == transect,"value"]
}
error_decision_tree <- rpart(value ~., data = decision_tree_data_error)
rpart.plot(error_decision_tree)
plotcp(error_decision_tree)






# Model difference vs. transect environmental variables

#non_habitat_variables <- c("CorinePatchDensity", "NaturaPatchDensity", "Rainfall", "Temperature")
non_habitat_variables <- c("")
long_data_error <- transect_long_df[transect_long_df$metric == "error_difference",]
decision_tree_data_error <- transect_habitat_variables[,setdiff(colnames(transect_habitat_variables), non_habitat_variables)]
decision_tree_data_error$value <- rep("no difference", nrow(decision_tree_data_error))
for (transect in rownames(decision_tree_data_error)) {
    if (long_data_error[long_data_error$transect == transect,"value"] > 0) {
        decision_tree_data_error[transect,]$value <- "lc lower"
    } else if (long_data_error[long_data_error$transect == transect,"value"] < 0) {
        decision_tree_data_error[transect,]$value <- "hab lower"
    } 
}
decision_tree_data_error$value <- factor(decision_tree_data_error$value)
error_decision_tree <- rpart(value ~., data = decision_tree_data_error)

par(mfrow = c(1, 2))
rpart.plot(error_decision_tree)
plotcp(error_decision_tree)


#non_habitat_variables <- c("CorinePatchDensity", "NaturaPatchDensity", "Rainfall", "Temperature")
long_data_uncertainty <- transect_long_df[transect_long_df$metric == "uncertainty_difference",]
decision_tree_data_uncertainty <- transect_habitat_variables[,setdiff(colnames(transect_habitat_variables), non_habitat_variables)]
decision_tree_data_uncertainty$value <- rep("no difference", nrow(decision_tree_data_uncertainty))
for (transect in rownames(decision_tree_data_uncertainty)) {
    if (long_data_uncertainty[long_data_uncertainty$transect == transect,"value"] > 0) {
        decision_tree_data_uncertainty[transect,]$value <- "lc lower"
    } else if (long_data_uncertainty[long_data_uncertainty$transect == transect,"value"] < 0) {
        decision_tree_data_uncertainty[transect,]$value <- "hab lower"
    } 
}
decision_tree_data_uncertainty$value <- factor(decision_tree_data_uncertainty$value)
uncertainty_decision_tree <- rpart(value ~., data = decision_tree_data_uncertainty)
rpart.plot(uncertainty_decision_tree)
plotcp(uncertainty_decision_tree)









# SPECIESWISE COMPARISON


species_long_df <- aggregate(overall_long_df,
                             value ~ species + metric,
                             FUN = median)









par(mfrow = c(2, 3))


plot_model_evaluation_ecdf(species_long_df, 
                           c("corine_occurrence_prob", "natura_occurrence_prob"), 
                           "CDF of expected occurrence probability E",
                           c(0, 1))

plot_model_evaluation_ecdf(species_long_df,
                           c("corine_error", "natura_error"), 
                           "CDF of mean absolute error MAE",
                           c(0, 1))

plot_model_evaluation_ecdf(species_long_df, 
                           c("corine_uncertainty", "natura_uncertainty"),
                           "CDF of uncertainty UC",
                           c(0, 1))


values <- species_long_df[species_long_df$metric == "occurrence_prob_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(species_long_df,
                           "occurrence_prob_difference",
                           "CDF of expected occurrence probablity difference \nfor species medians",
                           c(-limit, limit),
                           difference = TRUE)


values <- species_long_df[species_long_df$metric == "error_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(species_long_df,
                           "error_difference",
                           "CDF of mean absolute error difference \nfor species medians",
                           c(-limit, limit),
                           difference = TRUE)

values <- species_long_df[species_long_df$metric == "uncertainty_difference",]$value
limit <- max(abs(values))
plot_model_evaluation_ecdf(species_long_df,
                           "uncertainty_difference",
                           "CDF of uncertainty difference \nfor species medians",
                           c(-limit, limit),
                           difference = TRUE)













load(file = file.path(dir_data, "trait_data.RData"))
load(file = file.path(dir_data, "species_prevalences.RData"))
species_list <- species_long_df[species_long_df$metric == "error_difference",]$species

species_traits <- data.frame(species = species_list,
                             feeding = trait_data[species_list,]$Feeding,
                             mass = trait_data[species_list,]$Mass,
                             transect_prevalence = species_prevalences[species_list,]$transects)

rownames(species_traits) <- species_traits$species






long_data_error <- species_long_df[species_long_df$metric == "error_difference",]
decision_tree_data_error <- species_traits[,c("feeding", "mass", "transect_prevalence")]
decision_tree_data_error$value <- rep(0, nrow(decision_tree_data_error))
for (species in rownames(decision_tree_data_error)) {
    decision_tree_data_error[species,]$value <- long_data_error[long_data_error$species == species,"value"]
}
error_decision_tree <- rpart(value ~., data = decision_tree_data_error)
rpart.plot(error_decision_tree)
plotcp(error_decision_tree)


long_data_uncertainty <- species_long_df[species_long_df$metric == "uncertainty_difference",]
decision_tree_data_uncertainty <- species_traits[,c("feeding", "mass", "transect_prevalence")]
decision_tree_data_uncertainty$value <- rep(0, nrow(decision_tree_data_uncertainty))
for (species in rownames(decision_tree_data_uncertainty)) {
    decision_tree_data_uncertainty[species,]$value <- long_data_uncertainty[long_data_uncertainty$species == species,"value"]
}
uncertainty_decision_tree <- rpart(value ~., data = decision_tree_data_uncertainty)
rpart.plot(uncertainty_decision_tree)
plotcp(uncertainty_decision_tree)




# Species prediction decision tree


# Model difference vs. species traits

long_data_error <- species_long_df[species_long_df$metric == "error_difference",]
decision_tree_data_error <- species_traits[,c("feeding", "mass", "transect_prevalence")]
decision_tree_data_error$value <- rep("no difference", nrow(decision_tree_data_error))
for (species in rownames(decision_tree_data_error)) {
    if (long_data_error[long_data_error$species == species,"value"] > 0) {
        decision_tree_data_error[species,]$value <- "lc lower"
    } else if (long_data_error[long_data_error$species == species, "value"] < 0) {
        decision_tree_data_error[species,]$value <- "hab lower"
    } 
}
decision_tree_data_error$value <- factor(decision_tree_data_error$value)
error_decision_tree <- rpart(value ~., data = decision_tree_data_error)

par(mfrow = c(1, 2))
rpart.plot(error_decision_tree)
plotcp(error_decision_tree)



long_data_error <- species_long_df[species_long_df$metric == "uncertainty_difference",]
decision_tree_data_error <- species_traits[,c("feeding", "mass", "transect_prevalence")]
decision_tree_data_error$value <- rep("no difference", nrow(decision_tree_data_error))
for (species in rownames(decision_tree_data_error)) {
    if (long_data_error[long_data_error$species == species,"value"] > 0) {
        decision_tree_data_error[species,]$value <- "lc lower"
    } else if (long_data_error[long_data_error$species == species, "value"] < 0) {
        decision_tree_data_error[species,]$value <- "hab lower"
    } 
}
decision_tree_data_error$value <- factor(decision_tree_data_error$value)
error_decision_tree <- rpart(value ~., data = decision_tree_data_error)

par(mfrow = c(1, 2))
rpart.plot(error_decision_tree)
plotcp(error_decision_tree)
















