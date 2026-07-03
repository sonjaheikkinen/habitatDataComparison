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






# Differences between groups

# Divide transects into groups based on which model has lower error for them

transects_natura <- c()
transects_corine <- c()
transects_no_difference <- c()
limit <- c(-0.005, 0.005)
for (transect in unique(transect_long_df$transect)) {
    data_for_transect <- transect_long_df[transect_long_df$transect == transect,]
    error_difference <- data_for_transect[data_for_transect$metric == "error_difference",]$value
    if (error_difference < limit[1]) {
        transects_natura <- c(transects_natura, transect)
    } 
    if (error_difference > limit[2]) {
        transects_corine <- c(transects_corine, transect)
    }
    if (error_difference >= limit[1] & error_difference <= limit[2]) {
        transects_no_difference <- c(transects_no_difference, transect)
    }
}


# Compare environmental variables between groups

par(mfrow = c(2, 3))
for (variable in colnames(transect_habitat_variables)) {
    plot(density(transect_habitat_variables[transects_natura,][,variable]),
         main = "natura",
         xlab = variable)
    plot(density(transect_habitat_variables[transects_corine,][,variable]),
         main = "corine",
         xlab = variable)
    plot(density(transect_habitat_variables[transects_no_difference,][,variable]),
         main = "no_difference",
         xlab = variable)
    plot(ecdf(transect_habitat_variables[transects_natura,][,variable]),
         main = "natura",
         xlab = variable,
         verticals = TRUE,
         do.points = FALSE)
    abline(h = 0.5, col = "grey")
    plot(ecdf(transect_habitat_variables[transects_corine,][,variable]),
         main = "corine",
         xlab = variable,
         verticals = TRUE,
         do.points = FALSE)
    abline(h = 0.5, col = "grey")
    plot(ecdf(transect_habitat_variables[transects_no_difference,][,variable]),
         main = "no difference",
         xlab = variable,
         verticals = TRUE,
         do.points = FALSE)
    abline(h = 0.5, col = "grey")
}



# Plot environmental variable comparison
variable_data <- data.frame(variable = rep(colnames(transect_habitat_variables), each = 3),
                            group = rep(times = length(colnames(transect_habitat_variables)), 
                                        c("Natura smaller error", "Corine smaller error", "Error difference < 0.005")),
                            value = 0)
values <- c()
for (variable in colnames(transect_habitat_variables)) {
    natura_group_median <- median(transect_habitat_variables[transects_natura,][,variable])
    corine_group_median <- median(transect_habitat_variables[transects_corine,][,variable])
    no_difference_group_median <- median(transect_habitat_variables[transects_no_difference,][,variable])
    values <- c(values, natura_group_median, corine_group_median, no_difference_group_median)
}

variable_data$value <- values

ggplot(data = variable_data[variable_data$variable %in% setdiff(unique(variable_data$variable), 
                                                                c("Temperature",
                                                                "Rainfall",
                                                                "NaturaPatchDensity",
                                                                "CorinePatchDensity")),], 
       aes(x = value, 
           y = variable)) +
    geom_col(aes(fill = group),
             position = "dodge") +
    xlim(0, 0.6) +
    ggtitle("Median values for habitat variables in each group") +
    theme_minimal()









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
                             transect_prevalence = species_prevalences[species_list,]$transects,
                             mig = trait_data[species_list,]$Mig) 

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





# Differences between groups

# Divide species into groups based on which model has lower error for them

species_natura <- c()
species_corine <- c()
species_no_difference <- c()
limit <- c(-0.005, 0.005)
for (species in unique(species_long_df$species)) {
    data_for_species <- species_long_df[species_long_df$species == species,]
    error_difference <- data_for_species[data_for_species$metric == "error_difference",]$value
    print("")
    print(species)
    print(sprintf("Error difference: %s", error_difference))
    print(sprintf("Natura uncertainty: %s", data_for_species[data_for_species$metric == "natura_uncertainty",]$value))
    print(sprintf("Corine uncertainty: %s", data_for_species[data_for_species$metric == "corine_uncertainty",]$value))
    if (error_difference < limit[1]) {
        species_natura <- c(species_natura, species)
    } 
    if (error_difference > limit[2]) {
        species_corine <- c(species_corine, species)
    }
    if (error_difference >= limit[1] & error_difference <= limit[2]) {
        print("SMALL")
        species_no_difference <- c(species_no_difference, species)
    }
}


# Compare environmental variables between groups

par(mfrow = c(2, 3))


plot(density(species_traits[species_natura,][,"transect_prevalence"]),
     main = "natura",
     xlab = "transect prevalence")
plot(density(species_traits[species_corine,][,"transect_prevalence"]),
     main = "corine",
     xlab = "transect prevalence")
plot(density(species_traits[species_no_difference,][,"transect_prevalence"]),
     main = "no difference",
     xlab = "transect prevalence")
plot(ecdf(species_traits[species_natura,][,"transect_prevalence"]),
     main = "natura",
     xlab = "transect prevalence",
     verticals = TRUE,
     do.points = FALSE)
abline(h = 0.5, col = "grey")
plot(ecdf(species_traits[species_corine,][,"transect_prevalence"]),
     main = "corine",
     xlab = "transcet prevalence",
     verticals = TRUE,
     do.points = FALSE)
abline(h = 0.5, col = "grey")
plot(ecdf(species_traits[species_no_difference,][,"transect_prevalence"]),
     main = "no difference",
     xlab = "transect prevalence",
     verticals = TRUE,
     do.points = FALSE)
abline(h = 0.5, col = "grey")


par(mfrow = c(1, 3))

barplot(table(species_traits[species_natura,][,"feeding"]),
        main = "natura")
barplot(table(species_traits[species_corine,][,"feeding"]),
        main = "corine")
barplot(table(species_traits[species_no_difference,][,"feeding"]),
        main = "no difference")


barplot(table(species_traits[species_natura,][,"mig"]),
        main = "natura")
barplot(table(species_traits[species_corine,][,"mig"]),
        main = "corine")
barplot(table(species_traits[species_no_difference,][,"mig"]),
        main = "no difference")




# Plot trait comparison
trait_names <- c("feeding I", "feeding H", "feeding C", "feeding O", "feeding M",
                 "mig L", "mig S", "mig R",
                 "Transect prevalence")
variable_data <- data.frame(variable = rep(trait_names, each = 3),
                            group = rep(times = length(trait_names), 
                                        c("Natura smaller error", "Corine smaller error", "Error difference < 0.005")),
                            value = 0)
values <- c()
for (variable in trait_names) {
    if (variable == "Transect prevalence") {
        natura_group_value <- median(species_traits[species_natura,][,"transect_prevalence"])
        corine_group_value <- median(species_traits[species_corine,][,"transect_prevalence"])
        no_difference_group_value <- median(species_traits[species_no_difference,][,"transect_prevalence"])
    } else {
        variable_name <- strsplit(variable, " ")[[1]][1]
        variable_value <- strsplit(variable, " ")[[1]][2]
        data_for_variable <- species_traits[species_traits[,variable_name] == variable_value,]
        natura_group_value <- sum(species_natura %in% rownames(data_for_variable))
        corine_group_value <- sum(species_corine %in% rownames(data_for_variable))
        no_difference_group_value <- sum(species_no_difference %in% rownames(data_for_variable))
    }
    values <- c(values, natura_group_value, corine_group_value, no_difference_group_value)
}

variable_data$value <- values

ggplot(data = variable_data[variable_data$variable %in% setdiff(unique(variable_data$variable), 
                                                                c("Transect prevalence")),], 
       aes(x = value, 
           y = variable)) +
    geom_col(aes(fill = group),
             position = "dodge") +
    ggtitle("Number of species having specific traits in each group") +
    labs(x = "Number of species") +
    theme_minimal()




ggplot(data = variable_data[variable_data$variable == "Transect prevalence",], 
       aes(x = value, 
           y = variable)) +
    geom_col(aes(fill = group),
             position = "dodge") +
    ggtitle("Median transect prevalence in each group") +
    labs(x = "Transect prevalence") +
    theme_minimal()













