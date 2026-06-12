
# FUNCTIONS | PREDICTION


# THIS FUNCTION GENERATES A PREDICTION FOR GIVEN SPATIAL AND TEMPORAL COORDINATES
make_prediction <- function(prediction_units, variable_names) {
        
    # Set random noise to coordinates because 
    # predicting to the same exact coordinates does not work 
    new_coords <- data.frame(x = prediction_units$x + rnorm(nrow(prediction_units), mean = 0, sd = 0.01),
                             y = prediction_units$y + rnorm(nrow(prediction_units), mean = 0, sd = 0.01))
    
    # Similarly set random noise to years
    new_years <- data.frame(year = prediction_units$Year + rnorm(nrow(prediction_units), 
                                                                 mean = 0, 
                                                                 sd = 0.01))

    # Generate sample numbers
    new_sample <- data.frame(sample = 1:nrow(prediction_units))
    
    # Get environmental variables
    new_env_vars <- prediction_units[,variable_names, drop = FALSE]
    
    # Prepare prediction gradient over the generated transects
    gradient <- prepareGradient(model, 
                                XDataNew = new_env_vars, 
                                sDataNew = list(Sample = new_sample,
                                                Transect = new_coords,
                                                Year = new_years))
    
    # Predict occurrence over given gradient
    prediction <- predict(model,
                          Gradient = gradient,
                          expected = TRUE,
                          nParallel = 1)
    
    return(prediction)
    
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


# Download fitted models
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

# Save variable scales for variable effect predictions
variable_scales <- list()

# Save model names for later use
model_names <- c(strsplit(basename(fitted_models[1]), "\\.")[[1]][1],
                 strsplit(basename(fitted_models[2]), "\\.")[[1]][1])




# MAKE PREDICTIONS

for (model_number in 1:length(fitted_models)) {
        
    # Load fitted model and extract its name
    load(fitted_models[model_number])
    model_name <- model_names[model_number]
    model <- fitted_model
    

    # Load original environmental data and coordinates from the model
    env_vars <- data.frame(model$XData, check.names = FALSE)
    load(file = file.path(dir_data, "spatiotemporal_context.RData"))
    coords <- spatiotemporal_context
    unit_data <- env_vars
    unit_data$x <- coords$x
    unit_data$y <- coords$y
    unit_data$Year <- coords$Year
    unit_data$Transect <- coords$Transect
    
    
    # Make a prediction over the same coordinates as original data
    prediction <- make_prediction(unit_data, colnames(env_vars))
    save(prediction, file = file.path(dir_results, sprintf("prediction_%s.RData", model_name)))
    save(unit_data, file = file.path(dir_results, sprintf("unit_data_%s.RData", model_name)))
    
    
    # Predict variable effects
    variable_prediction_transect <- "000"
    variable_prediction_x <- mean(unit_data$x)
    variable_prediction_y <- mean(unit_data$y)
    variable_prediction_year <- round(mean(unit_data$Year))
    
    
    # Marginal effect - all else in their mean
    predicted_variable_effects <- list()
    for (focal_variable in colnames(env_vars)) {
        
        focal_variable_values <- unit_data[,focal_variable]
        
        if (!(focal_variable == "Temperature") 
            & !(focal_variable == "Rainfall")
            & !(focal_variable == "PatchDensity")) {
            focal_variable_values <- seq(0,
                                         1,
                                         length.out = 20)
        } else {
            focal_variable_values <- seq(min(focal_variable_values),
                                         max(focal_variable_values),
                                         length.out = 20)
        }
        data_length <- length(focal_variable_values)
        variable_scales[[model_name]][[focal_variable]] <- focal_variable_values
        
        unit_data_for_focal <- data.frame(Transect = rep(variable_prediction_transect, data_length),
                                          x = rep(variable_prediction_x, data_length),
                                          y = rep(variable_prediction_y, data_length),
                                          Year = rep(variable_prediction_year, data_length))
        unit_data_for_focal[,focal_variable] <- focal_variable_values
        
        other_variables <- setdiff(colnames(env_vars), focal_variable)
        
        for (non_focal_variable in other_variables) {
            values <- unit_data[,non_focal_variable]
            unit_data_for_focal[,non_focal_variable] <- rep(mean(values), data_length)
        }
        
        prediction <- make_prediction(unit_data_for_focal, colnames(env_vars))
        predicted_variable_effects[[focal_variable]] <- prediction
    }
    save(predicted_variable_effects, file = file.path(dir_results, sprintf("predicted_variable_effects_%s.RData", model_name)))
    
    
    # All else in their most likely value given the value of focal
    predicted_variable_effects_non_marginal <- list()
    for (focal_variable in colnames(env_vars)) {
        
        focal_variable_values <- unit_data[,focal_variable]
        
        if (!(focal_variable == "Temperature") 
            & !(focal_variable == "Rainfall")
            & !(focal_variable == "PatchDensity")) {
            focal_variable_values <- seq(0,
                                         1,
                                         length.out = 20)
        } else {
            focal_variable_values <- seq(min(focal_variable_values),
                                         max(focal_variable_values),
                                         length.out = 20)
        }
        data_length <- length(focal_variable_values)
        variable_scales[[model_name]][[focal_variable]] <- focal_variable_values
        
        unit_data_for_focal <- data.frame(Transect = rep(variable_prediction_transect, data_length),
                                          x = rep(variable_prediction_x, data_length),
                                          y = rep(variable_prediction_y, data_length),
                                          Year = rep(variable_prediction_year, data_length))
        unit_data_for_focal[,focal_variable] <- focal_variable_values
        
        other_variables <- setdiff(colnames(env_vars), focal_variable)
        
        for (non_focal_variable in other_variables) {
            
            values_focal <- unit_data[, focal_variable]
            values_nonfocal <- unit_data[, non_focal_variable]
            fitted_lm <- lm(values_nonfocal ~ values_focal)
            predicted_non_focal <- predict(fitted_lm,
                                           newdata = data.frame(values_focal = focal_variable_values))
            unit_data_for_focal[, non_focal_variable] <- predicted_non_focal
            
        }
        
        prediction <- make_prediction(unit_data_for_focal, colnames(env_vars))
        predicted_variable_effects_non_marginal[[focal_variable]] <- prediction
    }
    save(predicted_variable_effects_non_marginal, file = file.path(dir_results, sprintf("predicted_variable_effects_non_marginal_%s.RData", model_name)))
    save(variable_scales, file = file.path(dir_results, sprintf("variable_scales")))
}



# CALCULATE EXPECTED VALUES

# Load variable_scales
load(file.path(dir_results, "variable_scales"))


expected_values_list <- list()
expected_values_over_variables <- list()
expected_values_over_variables_non_marginal <- list()
predictions <- list()



# Calculate expected values from the predictions
for (model_number in 1:length(fitted_models)) {
    
    # Get model name
    model_name <- model_names[model_number]
    
    # Load prediction
    load(file.path(dir_results, sprintf("prediction_%s.RData", model_name)))
    predictions[[model_number]] <- prediction
    # Load predicted_variable_effects
    load(file.path(dir_results, sprintf("predicted_variable_effects_%s.RData", model_name)))
    # Load predicted_variable_effects_non_marginal
    load(file.path(dir_results, sprintf("predicted_variable_effects_non_marginal_%s.RData", model_name)))
    
    
        
    
    # The prediction is a list of 2000 samples (500 from each chain)
    # Each sample is matrix of occurrence probabilities 
    # and has dimensions sampling units x species
    
    # CALCULATE EXPECTED OCCURRENCE PROBABILITIES FOR EACH SPECIES AT EACH SAMPLING UNIT
    

    # Calculate averages by summing over all corresponding values in each matrix, 
    # then dividing by the number of matrices
    expected_values <- Reduce("+", prediction) / length(prediction)
    expected_values_list[[model_name]] <- expected_values
    
    
    expected_values <- list()
    for (variable in names(predicted_variable_effects)) {
        prediction_for_variable <- predicted_variable_effects[[variable]]
        expected_values[[variable]] <- Reduce("+", prediction_for_variable) / length(prediction_for_variable)
    }
    expected_values_over_variables[[model_name]] <- expected_values
    
    expected_values <- list()
    for (variable in names(predicted_variable_effects_non_marginal)) {
        prediction_for_variable <- predicted_variable_effects_non_marginal[[variable]]
        expected_values[[variable]] <- Reduce("+", prediction_for_variable) / length(prediction_for_variable)
    }
    expected_values_over_variables_non_marginal[[model_name]] <- expected_values

}


# FORMAT PREDICTION INFORMATION FOR COMPARISON


# Load data
load(fitted_models[[1]])
occurrence <- fitted_model$Y
load(file = file.path(dir_data, "spatiotemporal_context.RData"))
load(file.path(dir_data, "species_prevalences.RData"))
load(file.path(dir_results, sprintf("unit_data_%s.RData", model_names[[1]])))
corine_units <- unit_data
load(file.path(dir_results, sprintf("unit_data_%s.RData", model_names[[2]])))
natura_units <- unit_data



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





# Format one common unit data combined from the unit data of both models
unit_data <- natura_units
unit_data$NaturaPatchDensity <- natura_units$PatchDensity
unit_data$CorinePatchDensity <- corine_units$PatchDensity
unit_data$PatchDensity <- NULL
for (variable in setdiff(corine_habitat_variables, "CorinePatchDensity")) {
    unit_data[,variable] <- corine_units[,variable]
}




# Make one big dataframe for all prediction information
number_of_species <- ncol(occurrence)
number_of_samples <- nrow(spatiotemporal_context)
prediction_dataframe <- data.frame(transect = rep(spatiotemporal_context$Transect, number_of_species),
                                   x = rep(corine_units$x, number_of_species),
                                   y = rep(spatiotemporal_context$y, number_of_species),
                                   year = rep(spatiotemporal_context$Year, number_of_species),
                                   sample_number = rep(1:nrow(spatiotemporal_context), number_of_species),
                                   species = rep(colnames(occurrence), each = number_of_samples),
                                   true_value = as.vector(occurrence),
                                   corine_occurrence_prob = as.vector(expected_values_list[[1]]),
                                   natura_occurrence_prob = as.vector(expected_values_list[[2]]))
prediction_dataframe$occurrence_prob_difference <- prediction_dataframe$natura_occurrence_prob - prediction_dataframe$corine_occurrence_prob
prediction_dataframe$average_occurrence_prob <- rowMeans(prediction_dataframe[,c("natura_occurrence_prob", "corine_occurrence_prob")])
prediction_dataframe$occurrence_prob_relative_difference <- prediction_dataframe$occurrence_prob_difference / prediction_dataframe$average_occurrence_prob





# Add model uncertainty to the dataframe
uncertainties <- list()
for (model_number in 1:length(predictions)) {
    model_predictions <- predictions[[model_number]]
    model_predictions_array <- simplify2array(model_predictions)
    credible_interval_lower_bound <- apply(model_predictions_array, 
                                           c(1, 2), 
                                           quantile, 
                                           probs = 0.025)
    credible_interval_upper_bound <- apply(model_predictions_array, 
                                           c(1, 2), 
                                           quantile, 
                                           probs = 0.975)
    uncertainty <- credible_interval_upper_bound - credible_interval_lower_bound
    uncertainties[[model_number]] <- uncertainty
}
names(uncertainties) <- c("corine", "natura")
prediction_dataframe$corine_uncertainty <- rep(0, nrow(prediction_dataframe))
prediction_dataframe$natura_uncertainty <- rep(0, nrow(prediction_dataframe))
for (row in 1:nrow(prediction_dataframe)) {
    row_species <- prediction_dataframe[row,]$species
    row_sample <- prediction_dataframe[row,]$sample_number
    prediction_dataframe[row,]$corine_uncertainty <- uncertainties$corine[row_sample,row_species]
    prediction_dataframe[row,]$natura_uncertainty <- uncertainties$natura[row_sample,row_species]
}
prediction_dataframe$uncertainty_difference <- prediction_dataframe$natura_uncertainty - prediction_dataframe$corine_uncertainty
prediction_dataframe$uncertainty_relative_difference <- prediction_dataframe$uncertainty_difference / prediction_dataframe$average_occurrence_prob
prediction_dataframe$average_uncertainty <- rowMeans(prediction_dataframe[,c("natura_uncertainty", "corine_uncertainty")])


# Add  model error to the dataframe
prediction_dataframe$corine_error <- abs(prediction_dataframe$true_value - prediction_dataframe$corine_occurrence_prob)
prediction_dataframe$natura_error <- abs(prediction_dataframe$true_value - prediction_dataframe$natura_occurrence_prob)
prediction_dataframe$error_difference <- prediction_dataframe$natura_error - prediction_dataframe$corine_error
prediction_dataframe$error_relative_difference <- prediction_dataframe$error_difference / prediction_dataframe$average_occurrence_prob
prediction_dataframe$average_error <- rowMeans(prediction_dataframe[,c("natura_error", "corine_error")])



# Average dataframe data over years to get rid of temporal replicates
averaged <- aggregate(prediction_dataframe,
                      cbind(x,
                            y, 
                            corine_occurrence_prob,
                            natura_occurrence_prob,
                            occurrence_prob_difference,
                            occurrence_prob_relative_difference,
                            average_occurrence_prob,
                            corine_uncertainty,
                            natura_uncertainty,
                            uncertainty_difference,
                            average_uncertainty,
                            uncertainty_relative_difference,
                            corine_error,
                            natura_error,
                            error_difference,
                            error_relative_difference,
                            average_error) ~ species + transect,
                            FUN = mean)


# Reshape data to long format
overall_long_df <- reshape(averaged,
                           varying = setdiff(names(averaged), c("species", 
                                                                "transect", 
                                                                "x", 
                                                                "y")),
                           v.names = "value",
                           timevar = "metric",
                           times = setdiff(names(averaged), c("species", 
                                                              "transect", 
                                                              "x", 
                                                              "y")),
                           idvar = c("species", "transect"),
                           direction = "long")





# OVERALL MODEL COMPARISON



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











    


















