
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



# SCRIPT STARTS


# Download fitted models
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

# Save variable scales for variable effect predictions
variable_scales <- list()

# Save model names for later use
model_names <- c(strsplit(basename(fitted_models[1]), "\\.")[[1]][1],
                 strsplit(basename(fitted_models[2]), "\\.")[[1]][1])


# Make predictions for both models
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
    
    
    # MAKE A PREDICTION OVER THE SAME COORDINATES AS ORIGINAL DATA
    prediction <- make_prediction(unit_data, colnames(env_vars))
    save(prediction, file = file.path(dir_results, sprintf("prediction_%s.RData", model_name)))
    save(unit_data, file = file.path(dir_results, sprintf("unit_data_%s.RData", model_name)))
    
    
    
    # PREDICT VARIABLE EFFECTS
    
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


# Load original occurrence data
load(fitted_models[[1]])
occurrence <- fitted_model$Y

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



# Load spatiotemporal context
load(file = file.path(dir_data, "spatiotemporal_context.RData"))


# Load trait data
load(file.path(dir_data, "species_prevalences.RData"))


# PREDICTIONS 

# Load unit_data for predictions
load(file.path(dir_results, sprintf("unit_data_%s.RData", model_names[[1]])))
corine_units <- unit_data
load(file.path(dir_results, sprintf("unit_data_%s.RData", model_names[[2]])))
natura_units <- unit_data

unit_data <- natura_units
unit_data$NaturaPatchDensity <- natura_units$PatchDensity
unit_data$CorinePatchDensity <- corine_units$PatchDensity
unit_data$PatchDensity <- NULL
for (variable in setdiff(corine_habitat_variables, "CorinePatchDensity")) {
    unit_data[,variable] <- corine_units[,variable]
}

# MAKE ONE BIG DATAFRAME FOR ALL PREDICTION INFORMATION
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





# ANALYSE UNCERTAINTY AND CORRECTNESS OF PREDICTIONS TO THE DATAFRAME

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


prediction_dataframe$corine_error <- abs(prediction_dataframe$true_value - prediction_dataframe$corine_occurrence_prob)
prediction_dataframe$natura_error <- abs(prediction_dataframe$true_value - prediction_dataframe$natura_occurrence_prob)
prediction_dataframe$error_difference <- prediction_dataframe$natura_error - prediction_dataframe$corine_error
prediction_dataframe$error_relative_difference <- prediction_dataframe$error_difference / prediction_dataframe$average_occurrence_prob
prediction_dataframe$average_error <- rowMeans(prediction_dataframe[,c("natura_error", "corine_error")])



# AVERAGE OVER YEARS

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



plot(overall_long_df[overall_long_df$metric == "natura_occurrence_prob",]$value,
     overall_long_df[overall_long_df$metric == "natura_error",]$value,
     xlab = "Occurrence probability",
     ylab = "Mean absolute error",
     main = "E (hab) vs. MAE (hab)")

plot(overall_long_df[overall_long_df$metric == "natura_occurrence_prob",]$value,
     overall_long_df[overall_long_df$metric == "natura_uncertainty",]$value,
     xlab = "Occurrence probability",
     ylab = "Uncertainty",
     main = "E (hab) vs. UC (hab)")


plot(overall_long_df[overall_long_df$metric == "natura_uncertainty",]$value,
     overall_long_df[overall_long_df$metric == "natura_error",]$value,
     xlab = "Uncertainty",
     ylab = "Mean absolute error",
     main = "UC (hab) vs. MAE (hab)")


plot(overall_long_df[overall_long_df$metric == "corine_occurrence_prob",]$value,
     overall_long_df[overall_long_df$metric == "occurrence_prob_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Occurrence probability difference",
     main = "E (hab) vs. E diff")


plot(overall_long_df[overall_long_df$metric == "natura_occurrence_prob",]$value,
     overall_long_df[overall_long_df$metric == "error_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Mean absolute error difference",
     main = "E (hab) vs. MAE diff")

plot(overall_long_df[overall_long_df$metric == "natura_occurrence_prob",]$value,
     overall_long_df[overall_long_df$metric == "uncertainty_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Uncertainty difference",
     main = "E (hab) vs. UC diff")

plot(overall_long_df[overall_long_df$metric == "uncertainty_difference",]$value,
     overall_long_df[overall_long_df$metric == "error_difference",]$value,
     xlab = "Uncertainty difference",
     ylab = "Mean absolute error difference",
     main = "UC diff vs. MAE diff")


cor(overall_long_df[overall_long_df$metric == "uncertainty_difference",]$value,
    overall_long_df[overall_long_df$metric == "error_difference",]$value)


par(mfrow = c(2, 3))


plot(ecdf(overall_long_df[overall_long_df$metric == "corine_occurrence_prob",]$value),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of expected occurrence probability E",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(overall_long_df[overall_long_df$metric == "natura_occurrence_prob",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
     lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("E (hab)", "E (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)


plot(ecdf(overall_long_df[overall_long_df$metric == "corine_error",]$value),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of mean absolute error MAE",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(overall_long_df[overall_long_df$metric == "natura_error",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
      lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("MAE (hab)", "MAE (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)


plot(ecdf(overall_long_df[overall_long_df$metric == "corine_uncertainty",]$value),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of uncertainty UC",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(overall_long_df[overall_long_df$metric == "natura_uncertainty",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
      lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("UC (hab)", "UC (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)



values <- overall_long_df[overall_long_df$metric == "occurrence_prob_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of E difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))


values <- overall_long_df[overall_long_df$metric == "error_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of MAE difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))


values <- overall_long_df[overall_long_df$metric == "uncertainty_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of UC difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))












par(mfrow=c(1, 3))
hist(overall_long_df[overall_long_df$metric == "occurrence_prob_relative_difference",]$value,
     main = "E difference scaled with E",
     xlab = "Value",
     col = "forestgreen")
hist(overall_long_df[overall_long_df$metric == "error_relative_difference",]$value,
     main = "MAE difference scaled with E",
     xlab = "Value",
     col = "red")
hist(overall_long_df[overall_long_df$metric == "uncertainty_relative_difference",]$value,
     main = "UC difference scaled with E",
     xlab = "Value",
     col = "blue")


# TRANSECTWISE ANALYSIS

transect_long_df <- aggregate(overall_long_df,
                              cbind(value, x, y) ~ transect + metric,
                              FUN = median)


par(mfrow = c(2, 4))

plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "corine_occurrence_prob",]$value,
     xlab = "Occurrence probability (hab)",
     ylab = "Occurrence probability (lc)",
     main = "E (hab) vs. E (lc)")


plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "natura_error",]$value,
     xlab = "Occurrence probability",
     ylab = "Mean absolute error",
     main = "E (hab) vs. MAE (hab)")

plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "natura_uncertainty",]$value,
     xlab = "Occurrence probability",
     ylab = "Uncertainty",
     main = "E (hab) vs. UC (hab)")


plot(transect_long_df[transect_long_df$metric == "natura_uncertainty",]$value,
     transect_long_df[transect_long_df$metric == "natura_error",]$value,
     xlab = "Uncertainty",
     ylab = "Mean absolute error",
     main = "UC (hab) vs. MAE (hab)")


plot(transect_long_df[transect_long_df$metric == "corine_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "occurrence_prob_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Occurrence probability difference",
     main = "E (hab) vs. E diff")


plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "error_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Mean absolute error difference",
     main = "E (hab) vs. MAE diff")

plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "uncertainty_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Uncertainty difference",
     main = "E (hab) vs. UC diff")

plot(transect_long_df[transect_long_df$metric == "uncertainty_difference",]$value,
     transect_long_df[transect_long_df$metric == "error_difference",]$value,
     xlab = "Uncertainty difference",
     ylab = "Mean absolute error difference",
     main = "UC diff vs. MAE diff")






par(mfrow = c(2, 3))


plot(ecdf(transect_long_df[transect_long_df$metric == "corine_occurrence_prob",]$value),
     col = "blue",
     lwd = 3,
     yaxt = "n",
     verticals = TRUE,
     do.points = FALSE,
     main = "CDF of expected occurrence probability E",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
      lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("E (hab)", "E (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)


plot(ecdf(transect_long_df[transect_long_df$metric == "corine_error",]$value),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of mean absolute error MAE",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(transect_long_df[transect_long_df$metric == "natura_error",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
      lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("MAE (hab)", "MAE (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)


plot(ecdf(transect_long_df[transect_long_df$metric == "corine_uncertainty",]$value),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of uncertainty UC",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(transect_long_df[transect_long_df$metric == "natura_uncertainty",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
      lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("UC (hab)", "UC (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)



values <- transect_long_df[transect_long_df$metric == "occurrence_prob_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of E difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))


values <- transect_long_df[transect_long_df$metric == "error_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of MAE difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))


values <- transect_long_df[transect_long_df$metric == "uncertainty_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of UC difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))










par(mfrow = c(1, 2))

plot(ecdf(transect_long_df[transect_long_df$metric == "corine_occurrence_prob",]$value),
     col = "forestgreen",
     do.points = FALSE,
     verticals = TRUE,
     lwd = 3,
     xlim = c(0, 1),
     yaxt = "n",
     main = "CDF of E, MAE and UC for both models")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value),
      col = "lightgreen",
      lwd = 3,
      do.points = FALSE,
      verticals = TRUE)
lines(ecdf(transect_long_df[transect_long_df$metric == "corine_error",]$value),
      col = "red",
      lwd = 3,
      do.points = FALSE,
      verticals = TRUE)
lines(ecdf(transect_long_df[transect_long_df$metric == "natura_error",]$value),
      col = "pink",
      lwd = 3,
      do.points = FALSE,
      verticals = TRUE)
lines(ecdf(transect_long_df[transect_long_df$metric == "corine_uncertainty",]$value),
      col = "blue",
      lwd = 3,
      do.points = FALSE,
      verticals = TRUE)
lines(ecdf(transect_long_df[transect_long_df$metric == "natura_uncertainty",]$value),
      col = "lightblue",
      lwd = 3,
      do.points = FALSE,
      verticals = TRUE)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("E (hab)", "E (lc)",
                  "MAE (hab)", "MAE (lc)",
                  "UC (hab)", "UC (lc"),
       col = c("lightgreen", "forestgreen",
               "pink", "red",
               "lightblue", "blue"),
       lty = 1)



par(mfrow = c(1, 3))


plot(ecdf(transect_long_df[transect_long_df$metric == "occurrence_prob_difference",]$value),
     col = "blue",
     lwd = 3,
     do.points = FALSE,
     verticals = TRUE,
     yaxt = "n",
     main = "CDF of expected occurrence probablity difference \nfor transect medians",
     xlab = "Difference",
     ylab = "Cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))

plot(ecdf(transect_long_df[transect_long_df$metric == "error_difference",]$value),
     col = "blue",
     lwd = 3,
     do.points = FALSE,
     verticals = TRUE,
     yaxt = "n",
     main = "CDF of mean absolute error difference \nfor transect medians",
     xlab = "Difference",
     ylab = "Cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))


plot(ecdf(transect_long_df[transect_long_df$metric == "uncertainty_difference",]$value),
     col = "blue",
     lwd = 3,
     do.points = FALSE,
     verticals = TRUE,
     yaxt = "n",
     main = "CDF of uncertainty difference \nfor transect medians",
     xlab = "Difference",
     ylab = "Cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))










plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "natura_error",]$value,
     xlab = "Occurrence probability",
     ylab = "Error",
     main = "E (hab) vs. MAE (hab)")


plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "error_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Error difference",
     main = "E (hab) vs. MAE diff")





transect_list <- transect_long_df[transect_long_df$metric == "natura_error",]$transect
x_list <- transect_long_df[transect_long_df$metric == "natura_error",]$x
y_list <- transect_long_df[transect_long_df$metric == "natura_error",]$y


calculate_iqr_median <- function(values) {
    iqr <- quantile(values, 0.75) - quantile(values, 0.25)
    iqr_median <- iqr / median(values)
    return(iqr_median)
}

natura_trrv_E <- c()
corine_trrv_E <- c()
natura_trrv_MAE <- c()
corine_trrv_MAE <- c()
natura_trrv_UC <- c()
corine_trrv_UC <- c()
for (transect in transect_list) {
    
    natura_probs <- prediction_dataframe[prediction_dataframe$transect == transect,]$natura_occurrence_prob
    natura_trrv_E <- c(natura_trrv_E, calculate_iqr_median(natura_probs))
    corine_probs <- prediction_dataframe[prediction_dataframe$transect == transect,]$corine_occurrence_prob
    corine_trrv_E <- c(corine_trrv_E, calculate_iqr_median(corine_probs))
    
    natura_error <- prediction_dataframe[prediction_dataframe$transect == transect,]$natura_error
    natura_trrv_MAE <- c(natura_trrv_MAE, calculate_iqr_median(natura_error))
    corine_error <- prediction_dataframe[prediction_dataframe$transect == transect,]$corine_error
    corine_trrv_MAE <- c(corine_trrv_MAE, calculate_iqr_median(corine_error))
    
    natura_uc <- prediction_dataframe[prediction_dataframe$transect == transect,]$natura_uncertainty
    natura_trrv_UC <- c(natura_trrv_UC, calculate_iqr_median(natura_uc))
    corine_uc <- prediction_dataframe[prediction_dataframe$transect == transect,]$corine_uncertainty
    corine_trrv_UC <- c(corine_trrv_UC, calculate_iqr_median(corine_uc))
    
    
}



transect_long_df <- rbind(transect_long_df, data.frame(transect = rep(transect_list, 6),
                                     metric = c(rep("natura_trrv_E", length(transect_list)),
                                                rep("corine_trrv_E", length(transect_list)),
                                                rep("natura_trrv_UC", length(transect_list)),
                                                rep("corine_trrv_UC", length(transect_list)),
                                                rep("natura_trrv_MAE", length(transect_list)),
                                                rep("corine_trrv_MAE", length(transect_list))),
                                     value = c(natura_trrv_E, corine_trrv_E,
                                               natura_trrv_UC, corine_trrv_UC,
                                               natura_trrv_MAE, corine_trrv_MAE),
                                     x = rep(x_list, 6),
                                     y = rep(y_list, 6)))



par(mfrow = c(1, 2))

plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "natura_error",]$value,
     xlab = "Occurrence probability",
     ylab = "Mean absolute error",
     main = "mtr_E vs. mtr_MAE (hab)")

plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "natura_uncertainty",]$value,
     xlab = "Occurrence probability",
     ylab = "Uncertainty",
     main = "mtr_E vs. mtr_UC (hab)")

plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "error_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Mean absolute error difference",
     main = "E (hab) vs. MAE diff")

plot(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "uncertainty_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Uncertainty difference",
     main = "E (hab) vs. UC diff")



cor(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "natura_error",]$value)


cor(transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]$value,
     transect_long_df[transect_long_df$metric == "natura_uncertainty",]$value)





create_transect_difference_plot <- function(df,
                                            type,
                                            title,
                                            relative_difference) {
    relative_string <- ""
    if (relative_difference == TRUE) {
        title <- sprintf("%s RELATIVE", title)
        relative_string <- "_relative"
    }
    
    title <- sprintf("%s DIFFERENCE", title)
    
    metric_name <- sprintf("%s%s_difference", type, relative_string)
    metric_values <- df[df$metric == metric_name,]$value
    
    plot <- ggplot(df[df$metric == metric_name,], 
                   aes(x = x, 
                       y = y, 
                       color = factor(value > 0, levels = c(FALSE, TRUE)))) +
        geom_point(size = 2 + 3 * abs((metric_values / max(metric_values)))) +
        scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                           labels = c("Land cover", "Habitat type"),
                           name = sprintf("Model with higher %s", type),
                           drop = FALSE) +
        labs(title = sprintf("%s\n\n%s/%s negative vs. positive differences | min %s | max %s",
                             title,
                             sum(metric_values < 0),
                             sum(metric_values > 0),
                             round(min(metric_values), digits = 3),
                             round(max(metric_values), digits = 3))) +
        theme_minimal() +
        theme(plot.title = element_text(size = 10)) +
        coord_equal() 
    
    return(plot)
}

median_difference_plot <- create_transect_difference_plot(transect_long_df, "occurrence_prob", "OCCURRENCE PROBABILITY", FALSE)
uncertainty_difference_plot <- create_transect_difference_plot(transect_long_df, "uncertainty", "UNCERTAINTY", FALSE)
error_difference_plot <- create_transect_difference_plot(transect_long_df, "error", "ERROR", FALSE)

grid.arrange(median_difference_plot,
             uncertainty_difference_plot,
             error_difference_plot,
             nrow = 2)


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


par(mfrow = c(1, 2))

spatial_plot_data <- transect_long_df[transect_long_df$metric == "natura_occurrence_prob",]
plot(spatial_plot_data$x,
     spatial_plot_data$y,
     cex = 1 + 10 * spatial_plot_data$value,
     pch = 20,
     col = "forestgreen",
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "Occurrence probability mtr_E (hab)")
spatial_plot_data <- transect_long_df[transect_long_df$metric == "natura_trrv_E",]
plot(spatial_plot_data$x,
     spatial_plot_data$y,
     cex = spatial_plot_data$value,
     pch = 20,
     col = "forestgreen",
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "Occurrence probability relative variation imtr_E (hab)")



par(mfrow = c(1, 3))
spatial_plot_data <- transect_long_df[transect_long_df$metric == "occurrence_prob_difference",]
plot(spatial_plot_data$x,
     spatial_plot_data$y,
     pch = 20,
     col = "black",
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "Occurrence probability difference")
pos_data <- spatial_plot_data[spatial_plot_data$value > 0,]
points(pos_data$x,
       pos_data$y,
       cex = 2 + 50 * pos_data$value,
       pch = 20,
       col = "red")
neg_data <- spatial_plot_data[spatial_plot_data$value < 0,]
points(neg_data$x,
       neg_data$y,
       cex = 2 + 50 * abs(neg_data$value),
       pch = 20,
       col = "blue")
legend("bottomleft",
       legend = c("hab",
                  "lc"),
       col = c("red",
               "blue"),
       pch = 20,
       title = "Bigger value",
       bty = "n")

spatial_plot_data <- transect_long_df[transect_long_df$metric == "error_difference",]
plot(spatial_plot_data$x,
     spatial_plot_data$y,
     pch = 20,
     col = "black",
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "Mean absolute error difference")
pos_data <- spatial_plot_data[spatial_plot_data$value > 0,]
points(pos_data$x,
       pos_data$y,
       cex = 2 + 80 * pos_data$value,
       pch = 20,
       col = "red")
neg_data <- spatial_plot_data[spatial_plot_data$value < 0,]
points(neg_data$x,
       neg_data$y,
       cex = 2 + 80 * abs(neg_data$value),
       pch = 20,
       col = "blue")
legend("bottomleft",
       legend = c("hab",
                  "lc"),
       col = c("red",
               "blue"),
       pch = 20,
       title = "Bigger value",
       bty = "n")

spatial_plot_data <- transect_long_df[transect_long_df$metric == "uncertainty_difference",]
plot(spatial_plot_data$x,
     spatial_plot_data$y,
     pch = 20,
     col = "black",
     axes = FALSE,
     xlab = "",
     ylab = "",
     main = "Uncertainty difference")
pos_data <- spatial_plot_data[spatial_plot_data$value > 0,]
points(pos_data$x,
       pos_data$y,
       cex = 2 + 20 * pos_data$value,
       pch = 20,
       col = "red")
neg_data <- spatial_plot_data[spatial_plot_data$value < 0,]
points(neg_data$x,
       neg_data$y,
       cex = 2 + 20 * abs(neg_data$value),
       pch = 20,
       col = "blue")
legend("bottomleft",
       legend = c("hab",
                  "lc"),
       col = c("red",
               "blue"),
       pch = 20,
       title = "Bigger value",
       bty = "n")





natura_forest_types <- c("Luonnonmetsät", 
                         "Tunturikoivikot", 
                         "Lehdot", 
                         "Tulvametsät")

corine_forest_types <- c("Havumetsät.kivennäismaalla", 
                         "Sekametsät.kivennäismaalla", 
                         "Sekametsät.turvemaalla", 
                         "Lehtimetsät.kivennäismaalla", 
                         "Havumetsät.kalliomaalla")

other_types <- c("CorinePatchDensity", "NaturaPatchDensity",
                 "Temperature", "Rainfall")


metrics_to_plot <- c("occurrence_prob_difference",
                       "uncertainty_difference",
                       "error_difference")

variables <- c(natura_forest_types, corine_forest_types, other_types)

p_values_df <- data.frame(occurrence_prob_difference = rep(0, length(variables)),
                          uncertainty_difference = rep(0, length(variables)),
                          error_difference = rep(0, length(variables)))
estimates_df <- data.frame(occurrence_prob_difference = rep(0, length(variables)),
                          uncertainty_difference = rep(0, length(variables)),
                          error_difference = rep(0, length(variables)))



for (metric in metrics_to_plot) {
    
    p_values <- c()
    estimates <- c()
    
    
    for (variable in c(natura_forest_types, corine_forest_types, other_types)) {
        
        transect_habitat_values <- transect_habitat_variables[transect_list, variable]
        transect_metric_values <- transect_long_df[transect_long_df$metric == metric, c("value", "transect")]
        rownames(transect_metric_values) <- transect_metric_values$transect
        transect_metric_values <- transect_metric_values[transect_list,]
        
        
        model <- lm(transect_metric_values$value ~ transect_habitat_values)
        p_values <- c(p_values, summary(model)$coefficients[2,4])
        estimates <- c(estimates, summary(model)$coefficients[2,1])
    }
    
    names(p_values) <- c(natura_forest_types, corine_forest_types, other_types)
    p_values_df[,metric] <- p_values
    
    names(estimates) <- c(natura_forest_types, corine_forest_types, other_types)
    estimates_df[,metric] <- estimates
}

rownames(p_values_df) <- variables
colnames(p_values_df) <- c("Occurrence probability difference",
                           "Uncertainty difference",
                           "Error difference")

rownames(estimates_df) <- variables
colnames(estimates_df) <- c("Occurrence probability difference",
                           "Uncertainty difference",
                           "Error difference")


p_threshold <- 0.05 

heatmap_df <- expand.grid(variable = rownames(p_values_df), 
                          metric = colnames(p_values_df))

heatmap_df$p_value <- round(as.vector(as.matrix(p_values_df)), digits = 3)
heatmap_df$estimate <- round(as.vector(as.matrix(estimates_df)), digits = 3)


color_non_significant <- "grey"
color_positive <- "red"
color_negative <- "blue"


heatmap_df$fill_color <- ifelse(heatmap_df$p_value > 1, 
                                color_non_significant, 
                                ifelse(heatmap_df$estimate > 0, 
                                       color_positive,
                                       color_negative))



heatmap_df$fill_color <- factor(heatmap_df$fill_color,
                                levels = c(color_non_significant, 
                                           color_positive,
                                           color_negative))

heatmap_df$alpha <- ifelse(heatmap_df$p_value < 0.05,
                           1,
                           0.01)


ggplot(heatmap_df, aes(x = metric, y = variable)) +
    geom_tile(aes(fill = fill_color,
                  alpha = alpha)) +
    geom_text(aes(label = sprintf("P: %s, E: %s", p_value, estimate)), 
              size = 3,
              color = "gray30") +
    scale_fill_manual(values = c(grey = color_non_significant, 
                                 red = color_positive,
                                 blue = color_negative), 
                      name = "Direction", 
                      labels = c("Positive estimate", 
                                 "Negative estimate")) + 
    guides(alpha = "none") +
    labs(x = "",
         y = "",
         title = "Effects of transect environmental variables on difference metrics from linear models\n- p-value and estimate") +
    theme_minimal()



p_values_long <- data.frame(metric = rep(colnames(p_values_df), each = length(variables)),
                            variable = rep(variables, 3),
                            value = c(p_values_df[,1], p_values_df[,2], p_values_df[,3]))


ggplot(p_values_long, aes(x = metric, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    ylim(0, 1) +
    coord_flip() + 
    #geom_text(aes(label = round(value, 3)), 
    #          hjust = -0.2,
    #          size = 3,
    #          group = metric, 
    #          position = position_dodge(width = .9)) +
    geom_hline(yintercept = 0.05) +
    labs(title = "P-values",
         x = "",
         y = "",
         fill = "Variable"
    ) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme_minimal()





for (variable in c(natura_forest_types, corine_forest_types, other_types)) {
    habitat_values <- transect_habitat_variables[transect_list,variable]
    values <- list()
    for (name in variables_to_plot) {
        diff <- transect_long_df[transect_long_df$metric == name, c("value", "transect")]
        rownames(diff) <- diff$transect
        values[[name]] <- diff[transect_list,]$value  
    }
    plot(habitat_values,
         values[[1]],
         col = alpha("black", 0.5),
         pch = 3,
         ylim = c(min(unlist(values)), max(unlist(values))),
         xlab = "Amount of variable",
         ylab = "Difference",
         main = variable)
    abline(lm(values[[1]] ~ habitat_values), col = "black")
    points(habitat_values,
           values[[2]],
           pch = 20,
           cex = 2,
           col = alpha("dodgerblue", 0.25))
    abline(lm(values[[2]] ~ habitat_values), col = "dodgerblue")
    points(habitat_values,
           values[[3]],
           pch = 17,
           col = alpha("red", 0.25))
    abline(lm(values[[3]] ~ habitat_values), col = "red")
    abline(0, 0, col = "black", lty = "dotted")
    
}



# Transect prediction decision tree


long_data_error <- transect_long_df[transect_long_df$metric == "error_difference",]
decision_tree_data_error <- transect_habitat_variables[,setdiff(colnames(transect_habitat_variables), non_habitat_variables)]
decision_tree_data_error$value <- rep(0, nrow(decision_tree_data_error))
for (transect in rownames(decision_tree_data_error)) {
    decision_tree_data_error[transect,]$value <- long_data_error[long_data_error$transect == transect,"value"]
}
error_decision_tree <- rpart(value ~., data = decision_tree_data_error)
rpart.plot(error_decision_tree)
plotcp(error_decision_tree)


long_data_error <- transect_long_df[transect_long_df$metric == "uncertainty_difference",]
decision_tree_data_error <- transect_habitat_variables[,setdiff(colnames(transect_habitat_variables), non_habitat_variables)]
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


plot(ecdf(species_long_df[species_long_df$metric == "corine_occurrence_prob",]$value),
     col = "blue",
     lwd = 3,
     yaxt = "n",
     verticals = TRUE,
     do.points = FALSE,
     main = "CDF of expected occurrence probability E",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(species_long_df[species_long_df$metric == "natura_occurrence_prob",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
      lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("E (hab)", "E (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)


plot(ecdf(species_long_df[species_long_df$metric == "corine_error",]$value),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of mean absolute error MAE",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(species_long_df[species_long_df$metric == "natura_error",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
      lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("MAE (hab)", "MAE (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)


plot(ecdf(species_long_df[species_long_df$metric == "corine_uncertainty",]$value),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of uncertainty UC",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
lines(ecdf(species_long_df[species_long_df$metric == "natura_uncertainty",]$value),
      verticals = TRUE,
      do.points = FALSE,
      col = "lightblue",
      lwd = 3)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
legend("bottomright",
       legend = c("UC (hab)", "UC (lc)"),
       col = c("lightblue", "blue"),
       lty = 1)



values <- species_long_df[species_long_df$metric == "occurrence_prob_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of E difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))


values <- species_long_df[species_long_df$metric == "error_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of MAE difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))


values <- species_long_df[species_long_df$metric == "uncertainty_difference",]$value
limit <- max(abs(values))
plot(ecdf(values),
     verticals = TRUE,
     do.points = FALSE,
     col = "blue",
     lwd = 3,
     yaxt = "n",
     main = "CDF of UC difference",
     xlab = "value",
     ylab = "cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))











par(mfrow = c(1, 3))



plot(ecdf(species_long_df[species_long_df$metric == "occurrence_prob_difference",]$value),
     col = "blue",
     lwd = 3,
     do.points = FALSE,
     verticals = TRUE,
     yaxt = "n",
     main = "CDF of expected occurrence probablity difference \nfor species medians",
     xlab = "Difference",
     ylab = "Cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))

plot(ecdf(species_long_df[species_long_df$metric == "error_difference",]$value),
     col = "blue",
     lwd = 3,
     do.points = FALSE,
     verticals = TRUE,
     yaxt = "n",
     main = "CDF of mean absolute error difference \nfor species medians",
     xlab = "Difference",
     ylab = "Cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))


plot(ecdf(species_long_df[species_long_df$metric == "uncertainty_difference",]$value),
     col = "blue",
     lwd = 3,
     do.points = FALSE,
     verticals = TRUE,
     yaxt = "n",
     main = "CDF of uncertainty difference \nfor species medians",
     xlab = "Difference",
     ylab = "Cumulative probability")
axis (2, at = c(0, 0.25, 0.5, 0.75, 1), 
      labels = c(0, 0.25, 0.5, 0.75, 1),
      las = 2)
abline(0.5, 0, col = alpha("black", 0.25))
abline(0.25, 0, col = alpha("black", 0.25))
abline(0.75, 0, col = alpha("black", 0.25))
abline(v = 0, col = alpha("black", 0.25))









plot(species_long_df[species_long_df$metric == "natura_occurrence_prob",]$value,
     species_long_df[species_long_df$metric == "corine_occurrence_prob",]$value,
     xlab = "Occurrence probability (hab)",
     ylab = "Occurrence probability (lc)",
     main = "E (hab) vs. E (lc)")




par(mfrow = c(2, 4))

plot(species_long_df[species_long_df$metric == "natura_occurrence_prob",]$value,
     species_long_df[species_long_df$metric == "corine_occurrence_prob",]$value,
     xlab = "Occurrence probability (hab)",
     ylab = "Occurrence probability (lc)",
     main = "E (hab) vs. E (lc)")


plot(species_long_df[species_long_df$metric == "natura_occurrence_prob",]$value,
     species_long_df[species_long_df$metric == "natura_error",]$value,
     xlab = "Occurrence probability",
     ylab = "Mean absolute error",
     main = "E (hab) vs. MAE (hab)")

plot(species_long_df[species_long_df$metric == "natura_occurrence_prob",]$value,
     species_long_df[species_long_df$metric == "natura_uncertainty",]$value,
     xlab = "Occurrence probability",
     ylab = "Uncertainty",
     main = "E (hab) vs. UC (hab)")


plot(species_long_df[species_long_df$metric == "natura_uncertainty",]$value,
     species_long_df[species_long_df$metric == "natura_error",]$value,
     xlab = "Uncertainty",
     ylab = "Mean absolute error",
     main = "UC (hab) vs. MAE (hab)")


plot(species_long_df[species_long_df$metric == "corine_occurrence_prob",]$value,
     species_long_df[species_long_df$metric == "occurrence_prob_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Occurrence probability difference",
     main = "E (hab) vs. E diff")


plot(species_long_df[species_long_df$metric == "natura_occurrence_prob",]$value,
     species_long_df[species_long_df$metric == "error_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Mean absolute error difference",
     main = "E (hab) vs. MAE diff")

plot(species_long_df[species_long_df$metric == "natura_occurrence_prob",]$value,
     species_long_df[species_long_df$metric == "uncertainty_difference",]$value,
     xlab = "Occurrence probability",
     ylab = "Uncertainty difference",
     main = "E (hab) vs. UC diff")

plot(species_long_df[species_long_df$metric == "uncertainty_difference",]$value,
     species_long_df[species_long_df$metric == "error_difference",]$value,
     xlab = "Uncertainty difference",
     ylab = "Mean absolute error difference",
     main = "UC diff vs. MAE diff")



load(file = file.path(dir_data, "trait_data.RData"))
load(file = file.path(dir_data, "species_prevalences.RData"))
species_list <- species_long_df[species_long_df$metric == "error_difference",]$species

species_traits <- data.frame(species = species_list,
                             feeding = trait_data[species_list,]$Feeding,
                             mass = trait_data[species_list,]$Mass,
                             transect_prevalence = species_prevalences[species_list,]$transects)


par(mfrow = c(1, 2))
for (variable in c("transect_prevalence", "mass")) {
    variable_values <- species_traits[,variable]
    variables_to_plot <- c("occurrence_prob_difference",
                           "uncertainty_difference",
                           "error_difference")
    values <- list()
    for (name in variables_to_plot) {
        diff <- species_long_df[species_long_df$metric == name, c("value", "species")]
        rownames(diff) <- diff$species
        values[[name]] <- diff[species_list,]$value  
    }
    plot(variable_values,
         values[[1]],
         col = alpha("black", 0.5),
         pch = 3,
         ylim = c(min(unlist(values)), max(unlist(values))),
         xlab = "Amount of variable",
         ylab = "Difference",
         main = variable)
    abline(lm(values[[1]] ~ variable_values), col = "black")
    points(variable_values,
           values[[2]],
           pch = 20,
           cex = 2,
           col = alpha("dodgerblue", 0.25))
    abline(lm(values[[2]] ~ variable_values), col = "dodgerblue")
    points(variable_values,
           values[[3]],
           pch = 17,
           col = alpha("red", 0.25))
    abline(lm(values[[3]] ~ variable_values), col = "red")
    abline(0, 0, col = "black", lty = "dotted")
}


species_long_df$feeding <- trait_data[species_long_df$species,]$Feeding

par(mfrow = c(1, 3))

boxplot(value ~ feeding,
        species_long_df[species_long_df$metric == "occurrence_prob_difference",],
        main = "E difference")
abline(0, 0, col = "blue")
boxplot(value ~ feeding,
        species_long_df[species_long_df$metric == "error_difference",],
        main = "MAE difference")
abline(0, 0, col = "blue")
boxplot(value ~ feeding,
        species_long_df[species_long_df$metric == "uncertainty_difference",],
        main = "UC difference")
abline(0, 0, col = "blue")





species_list <- unique(prediction_dataframe$species)


metrics_to_plot <- c("occurrence_prob_difference",
                     "uncertainty_difference",
                     "error_difference")

variables <- c("mass", "transect_prevalence")

p_values_df <- data.frame(occurrence_prob_difference = rep(0, length(variables)),
                          uncertainty_difference = rep(0, length(variables)),
                          error_difference = rep(0, length(variables)))
estimates_df <- data.frame(occurrence_prob_difference = rep(0, length(variables)),
                           uncertainty_difference = rep(0, length(variables)),
                           error_difference = rep(0, length(variables)))


rownames(species_traits) <- species_traits$species


for (metric in metrics_to_plot) {
    
    p_values <- c()
    estimates <- c()
    
    
    for (variable in variables) {
        
        species_trait_values <- species_traits[species_list, variable]
        species_metric_values <- species_long_df[species_long_df$metric == metric, c("value", "species")]
        rownames(species_metric_values) <- species_metric_values$species
        species_metric_values <- species_metric_values[species_list,]
        
        
        model <- lm(species_metric_values$value ~ species_trait_values)
        p_values <- c(p_values, summary(model)$coefficients[2,4])
        estimates <- c(estimates, summary(model)$coefficients[2,1])
    }
    
    names(p_values) <- variables
    p_values_df[,metric] <- p_values
    
    names(estimates) <- variables
    estimates_df[,metric] <- estimates
}

rownames(p_values_df) <- variables
colnames(p_values_df) <- c("Occurrence probability difference",
                           "Uncertainty difference",
                           "Error difference")

rownames(estimates_df) <- variables
colnames(estimates_df) <- c("Occurrence probability difference",
                            "Uncertainty difference",
                            "Error difference")


p_threshold <- 0.05 

heatmap_df <- expand.grid(variable = rownames(p_values_df), 
                          metric = colnames(p_values_df))


color_positive <- "red"
color_negative <- "blue"

heatmap_df$p_value <- as.vector(as.matrix(p_values_df))
heatmap_df$estimate <- as.vector(as.matrix(estimates_df))



heatmap_df$fill_color <- ifelse(heatmap_df$p_value > 1, 
                                color_non_significant, 
                                ifelse(heatmap_df$estimate > 0, 
                                       color_positive,
                                       color_negative))

heatmap_df$p_value <- round(as.vector(as.matrix(p_values_df)), digits = 3)
heatmap_df$estimate <- round(as.vector(as.matrix(estimates_df)), digits = 3)



heatmap_df$fill_color <- factor(heatmap_df$fill_color,
                                levels = c(color_non_significant, 
                                           color_positive,
                                           color_negative))

heatmap_df$alpha <- ifelse(heatmap_df$p_value < 0.05,
                           1,
                           0.01)


ggplot(heatmap_df, aes(x = metric, y = variable)) +
    geom_tile(aes(fill = fill_color,
                  alpha = alpha)) +
    geom_text(aes(label = sprintf("P: %s, E: %s", p_value, estimate)), 
              size = 3,
              color = "gray30") +
    scale_fill_manual(values = c(grey = color_non_significant, 
                                 red = color_positive,
                                 blue = color_negative), 
                      name = "Direction", 
                      labels = c("Positive estimate", 
                                 "Negative estimate")) + 
    guides(alpha = "none") +
    labs(x = "",
         y = "",
         title = "Effects of transect environmental variables on difference metrics from linear models\n- p-value and estimate") +
    theme_minimal()




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







# INTERESTING SPECIES AND TRANSECTS

natura_species <- c()
corine_species <- c()

for (species in species_list) {
    
    data <- species_long_df[species_long_df$species == species,]
    overall_data <- overall_long_df[overall_long_df$species == species,]
    
    plot(overall_data[overall_data$metric == "occurrence_prob_relative_difference",]$value,
         main = species)
    text(overall_data[overall_data$metric == "occurrence_prob_relative_difference",]$value,
         labels = overall_data[overall_data$metric == "occurrence_prob_relative_difference",]$transect)
    
    prob_diff <- data[data$metric == "occurrence_prob_relative_difference",]$value
    error_diff <- data[data$metric == "error_relative_difference",]$value
    uc_diff <- data[data$metric == "uncertainty_relative_difference",]$value
    
    if (prob_diff > 0 & (error_diff < 0 | uc_diff < 0)) {
        natura_species <- c(natura_species, species)
    }
    if (prob_diff < 0 & (error_diff > 0 | uc_diff > 0)) {
        corine_species <- c(corine_species, species)
    }
}


load(file = file.path(dir_data, "trait_data.RData"))
View(trait_data[natura_species,])
View(trait_data[corine_species,])




natura_transects <- c()
corine_transects <- c()

for (transect in unique(prediction_dataframe$transect)) {
    
    data <- transect_long_df[transect_long_df$transect == transect,]
    
    prob_diff <- data[data$metric == "occurrence_prob_relative_difference",]$value
    error_diff <- data[data$metric == "error_relative_difference",]$value
    uc_diff <- data[data$metric == "uncertainty_relative_difference",]$value
    
    if (prob_diff > 0 & (error_diff < 0 | uc_diff < 0)) {
        natura_transects <- c(natura_transects, transect)
    }
    if (prob_diff < 0 & (error_diff > 0 | uc_diff > 0)) {
        corine_transects <- c(corine_transects, transect)
    }
}

View(unit_data[unit_data$Transect == natura_transects,])
View(unit_data[unit_data$Transect == corine_transects,])


plot(unit_data$x,
     unit_data$y)
points(unit_data[unit_data$Transect == natura_transects,]$x,
     unit_data[unit_data$Transect == natura_transects,]$y,
     col = "red",
     pch = 20,
     cex = 4)
points(unit_data[unit_data$Transect == corine_transects,]$x,
       unit_data[unit_data$Transect == corine_transects,]$y,
       col = "blue",
       pch = 20,
       cex = 4)
text(unit_data$x,
     unit_data$y,
     labels = unit_data$Transect)








# ANALYSE DISTRIBUTION OF PREDICTIONS





plot(ecdf(prediction_dataframe$posterior_mean_relative_difference),
     xlab = "Difference between predicted values",
     ylab = "Probability",
     main = "CDF of value difference")
lines(ecdf(prediction_dataframe$posterior_mean_relative_difference))


plot(prediction_dataframe$natura_posterior_mean, 
     prediction_dataframe$natura_uncertainty)

plot(prediction_dataframe$corine_posterior_mean,
     (prediction_dataframe$corine_error - prediction_dataframe$corine_posterior_mean) / prediction_dataframe$corine_posterior_mean)

plot(prediction_dataframe$natura_error,
     prediction_dataframe$natura_uncertainty)


par(mfrow = c(3, 3))

boxplot(prediction_dataframe[,c("corine_posterior_mean", 
                                "natura_posterior_mean")],
        main = sprintf("Posterior means from both models, \n medians: corine %s, natura %s",
                       round(median(prediction_dataframe$corine_posterior_mean), digits = 3),
                       round(median(prediction_dataframe$natura_posterior_mean), digits = 3)))
boxplot(prediction_dataframe$posterior_mean_difference,
        main = sprintf("Difference as natura - corine, \n median %s",
                       round(median(prediction_dataframe$posterior_mean_difference), digits = 3)))  
boxplot(prediction_dataframe$posterior_mean_relative_difference,
        main = sprintf("Relative difference of natura compared to corine, \n median %s",
                       round(median(prediction_dataframe$posterior_mean_relative_difference), digits = 3))) 

boxplot(prediction_dataframe[,c("corine_uncertainty", 
                                "natura_uncertainty")],
        main = sprintf("Posterior uncertainty from both models, \n medians: corine %s, natura %s",
                       round(median(prediction_dataframe$corine_uncertainty), digits = 3),
                       round(median(prediction_dataframe$natura_uncertainty), digits = 3)))
boxplot(prediction_dataframe$uncertainty_difference,
        main = sprintf("Difference as natura - corine, \n median %s",
                       round(median(prediction_dataframe$uncertainty_difference), digits = 3)))  
boxplot(prediction_dataframe$uncertainty_relative_difference,
        main = sprintf("Relative difference of natura compared to corine, \n median %s",
                       round(median(prediction_dataframe$uncertainty_relative_difference), digits = 3))) 


boxplot(prediction_dataframe[,c("corine_error", 
                                "natura_error")],
        main = sprintf("Posterior error from both models, \n medians: corine %s, natura %s",
                       round(median(prediction_dataframe$corine_error), digits = 3),
                       round(median(prediction_dataframe$natura_error), digits = 3)))
boxplot(prediction_dataframe$error_difference,
        main = sprintf("Difference as natura - corine, \n median %s",
                       round(median(prediction_dataframe$error_difference), digits = 3)))  
boxplot(prediction_dataframe$error_relative_difference,
        main = sprintf("Relative difference of natura compared to corine, \n median %s",
                       round(median(prediction_dataframe$error_relative_difference), digits = 3))) 
 

species_list <- unique(prediction_dataframe$species)
species_postdis_dataframe <- data.frame(corine_median = rep(0, length(species_list)),
                                        natura_median = rep(0, length(species_list)),
                                        corine_iqr =rep(0, length(species_list)),
                                        natura_iqr = rep(0, length(species_list)),
                                        corine_iqr_median = rep(0, length(species_list)),
                                        natura_iqr_median = rep(0, length(species_list)),
                                        corine_uncertainty = rep(0, length(species_list)),
                                        natura_uncertainty = rep(0, length(species_list)),
                                        corine_error = rep(0, length(species_list)),
                                        natura_error = rep(0, length(species_list)),
                                        median_diff = rep(0, length(species_list)),
                                        iqr_diff = rep(0, length(species_list)),
                                        iqr_median_diff = rep(0, length(species_list)),
                                        uncertainty_diff = rep(0, length(species_list)),
                                        error_diff = rep(0, length(species_list)),
                                        median_relative_diff = rep(0, length(species_list)),
                                        iqr_relative_diff = rep(0, length(species_list)),
                                        iqr_median_relative_diff = rep(0, length(species_list)),
                                        uncertainty_relative_diff = rep(0, length(species_list)),
                                        error_relative_diff = rep(0, length(species_list)))
rownames(species_postdis_dataframe) <- species_list

for (species in unique(prediction_dataframe$species)) {
    data <- prediction_dataframe[prediction_dataframe$species == species,]
    corine_median <- median(data$corine_posterior_mean)
    natura_median <- median(data$natura_posterior_mean)
    corine_iqr <- quantile(data$corine_posterior_mean, 0.75) - quantile(data$corine_posterior_mean, 0.25)
    natura_iqr <- quantile(data$natura_posterior_mean, 0.75) - quantile(data$natura_posterior_mean, 0.25)
    corine_iqr_median <- corine_iqr / corine_median
    natura_iqr_median <- natura_iqr / natura_median
    corine_uncertainty <- median(data$corine_uncertainty)
    natura_uncertainty <- median(data$natura_uncertainty)
    corine_error <- median(data$corine_error)
    natura_error <- median(data$natura_error)
    median_diff <- natura_median - corine_median
    iqr_diff <- natura_iqr - corine_iqr
    iqr_median_diff <- natura_iqr_median - corine_iqr_median
    uncertainty_diff <- natura_uncertainty - corine_uncertainty
    error_diff <- natura_error - corine_error
    median_relative_diff <- (natura_median - corine_median) / corine_median
    iqr_relative_diff <- (natura_iqr - corine_iqr) / corine_iqr
    iqr_median_relative_diff <- (natura_iqr_median - corine_iqr_median) / corine_iqr_median
    uncertainty_relative_diff <- (natura_uncertainty - corine_uncertainty) / corine_uncertainty
    error_relative_diff <- (natura_error - corine_error) / corine_error
    species_postdis_dataframe[species,] <- c(corine_median, natura_median,
                                             corine_iqr, natura_iqr,
                                             corine_iqr_median, natura_iqr_median,
                                             corine_uncertainty, natura_uncertainty,
                                             corine_error, natura_error,
                                             median_diff, iqr_diff, iqr_median_diff, uncertainty_diff, error_diff,
                                             median_relative_diff, iqr_relative_diff, iqr_median_relative_diff, uncertainty_relative_diff, error_relative_diff)
}

species_postdis_dataframe$species <- rownames(species_postdis_dataframe)

plot(species_prevalences[species_list[1],]$transects,
     species_postdis_dataframe[species_list[1],]$natura_error,
     ylim = c(0, 1),
     xlim = c(0, 50))
for (species in species_list[2:length(species_list)]) {
    points(species_prevalences[species,]$transects,
           species_postdis_dataframe[species,]$natura_error)
}



true_data_rare <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0)
true_data_medium <- c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
true_data_common <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0)
prob_rare <- c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
prob_medium <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
prob_common <- c(0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 0.9)

res <- list()
res[["rare_rare"]] <- abs(true_data_rare - prob_rare)
res[["rare_medium"]] <- abs(true_data_rare - prob_medium)
res[["rare_common"]] <- abs(true_data_rare - prob_common)
res[["medium_rare"]] <- abs(true_data_medium - prob_rare)
res[["medium_medium"]] <- abs(true_data_medium- prob_medium)
res[["medium_common"]] <- abs(true_data_medium - prob_common)
res[["common_rare"]] <- abs(true_data_common - prob_rare)
res[["common_medium"]] <- abs(true_data_common - prob_medium)
res[["common_common"]] <- abs(true_data_common - prob_common)

res_avg <- lapply(res, FUN = mean)
res_median <- lapply(res, FUN = median)






long_df <- reshape(species_postdis_dataframe,
                   varying = names(species_postdis_dataframe)[names(species_postdis_dataframe) != "species"],
                   v.names = "value",
                   timevar = "metric",
                   times = names(species_postdis_dataframe)[names(species_postdis_dataframe) != "species"],
                   idvar = "species",
                   direction = "long")


ordering_rows <- long_df[long_df$metric == "natura_median",]
species_order <- ordering_rows[order(ordering_rows$value),]$species
long_df$species <- factor(long_df$species, levels = species_order)


create_specieswise_prediction_plot <- function(long_df, type, title, hide_species = FALSE) {
    
    relative_diffs <- long_df[long_df$metric == sprintf("%s_relative_diff", type),]$value
    abs_min_relative_diff <- abs(min(relative_diffs))
    max_relative_diff <- max(relative_diffs)
    
    plot <- ggplot(long_df[long_df$metric %in% c(sprintf("natura_%s", type), 
                                                 sprintf("corine_%s", type)),], 
                   aes(x = species, y = value, fill = metric)) +
            geom_col(position = "identity", alpha = 0.5) +
            geom_point(data = long_df[long_df$metric == sprintf("%s_relative_diff", type),],
                       aes(x = species,
                           y = (value + abs_min_relative_diff) / (max_relative_diff + abs_min_relative_diff))) +
            geom_hline(yintercept = abs_min_relative_diff / (max_relative_diff + abs_min_relative_diff)) +
            geom_text(data = long_df[long_df$metric == sprintf("%s_relative_diff", type),],
                      size = 3,
                      aes(x = species,
                          y = 1.2,
                          label = round(value, 3))) +
            ggtitle(sprintf("%s\n\nDifference | %s/%s negative vs. positive | min %s | max %s\nRelative difference | min % s | max %s",
                            title,
                            sum(long_df[long_df$metric == sprintf("%s_diff", type),]$value < 0),
                            sum(long_df[long_df$metric == sprintf("%s_diff", type),]$value > 0),
                            round(min(long_df[long_df$metric == sprintf("%s_diff", type),]$value), digits = 3),
                            round(max(long_df[long_df$metric == sprintf("%s_diff", type),]$value), digits = 3),
                            round(min(long_df[long_df$metric == sprintf("%s_relative_diff", type),]$value), digits = 3),
                            round(max(long_df[long_df$metric == sprintf("%s_relative_diff", type),]$value), digits = 3))) +
            theme_minimal() +
            theme(plot.title = element_text(size = 8),
                  axis.text.y= if (hide_species) {
                                    element_blank() 
                                } else {
                                    element_text()
                                }) +
            ylim(0, 1.3) +
            coord_flip()
    
    return(plot)
}


median_ratio_plot <- create_specieswise_prediction_plot(long_df, "median", "MEDIAN")
iqr_ratio_plot <- create_specieswise_prediction_plot(long_df, "iqr", "IQR", hide_species = TRUE)

iqr_median_ratio_plot <- ggplot(long_df[long_df$metric %in% c("natura_iqr_median", "corine_iqr_median"),], 
                                aes(x = species, y = value, fill = metric)) +
    geom_col(position = "identity", alpha = 0.5) +
    theme_minimal() +
    theme(plot.title = element_text(size = 8),
          axis.text.y=element_blank()) +
    ggtitle(sprintf("IQR / MEDIAN\n\n%s/%s negative vs. positive differences | min %s | max %s",
                    sum(long_df[long_df$metric == "iqr_median_diff",]$value < 0),
                    sum(long_df[long_df$metric == "iqr_median_diff",]$value > 0),
                    round(min(long_df[long_df$metric == "iqr_median_diff",]$value), digits = 3),
                    round(max(long_df[long_df$metric == "iqr_median_diff",]$value), digits = 3))) +
    coord_flip()


uncertainty_ratio_plot <- create_specieswise_prediction_plot(long_df, "uncertainty", "UNCERTAINTY")
error_ratio_plot <- create_specieswise_prediction_plot(long_df, "error", "ERROR", hide_species = TRUE)

grid.arrange(median_plot, iqr_ratio_plot, iqr_median_ratio_plot,
             uncertainty_ratio_plot, error_ratio_plot,
             ncol = 3, 
             widths = c(1.4, 1, 1))







transect_list <- unique(prediction_dataframe$transect)
transect_postdis_dataframe <- data.frame(corine_median = rep(0, length(transect_list)),
                                         natura_median = rep(0, length(transect_list)),
                                         corine_iqr = rep(0, length(transect_list)),
                                         natura_iqr = rep(0, length(transect_list)),
                                         corine_iqr_median = rep(0, length(transect_list)),
                                         natura_iqr_median = rep(0, length(transect_list)),
                                         corine_uncertainty = rep(0, length(transect_list)),
                                         natura_uncertainty = rep(0, length(transect_list)),
                                         corine_error = rep(0, length(transect_list)),
                                         natura_error = rep(0, length(transect_list)),
                                         median_diff = rep(0, length(transect_list)),
                                         iqr_diff = rep(0, length(transect_list)),
                                         iqr_median_diff = rep(0, length(transect_list)),
                                         uncertainty_diff = rep(0, length(transect_list)),
                                         error_diff = rep(0, length(transect_list)),
                                         median_relative_diff = rep(0, length(transect_list)),
                                         iqr_relative_diff = rep(0, length(transect_list)),
                                         iqr_median_relative_diff = rep(0, length(transect_list)),
                                         uncertainty_relative_diff = rep(0, length(transect_list)),
                                         error_relative_diff = rep(0, length(transect_list)),
                                         x = rep(0, length(transect_list)),
                                         y = rep(0, length(transect_list)))
rownames(transect_postdis_dataframe) <- transect_list



for (transect in unique(prediction_dataframe$transect)) {
    data <- prediction_dataframe[prediction_dataframe$transect == transect,]
    corine_median <- median(data$corine_posterior_mean)
    natura_median <- median(data$natura_posterior_mean)
    corine_iqr <- quantile(data$corine_posterior_mean, 0.75) - quantile(data$corine_posterior_mean, 0.25)
    natura_iqr <- quantile(data$natura_posterior_mean, 0.75) - quantile(data$natura_posterior_mean, 0.25)
    corine_iqr_median <- corine_iqr / corine_median
    natura_iqr_median <- natura_iqr / natura_median
    corine_uncertainty <- median(data$corine_uncertainty)
    natura_uncertainty <- median(data$natura_uncertainty)
    corine_error <- median(data$corine_error)
    natura_error <- median(data$natura_error)
    median_diff <- natura_median - corine_median
    iqr_diff <- natura_iqr - corine_iqr
    iqr_median_diff <- natura_iqr_median - corine_iqr_median
    uncertainty_diff <- natura_uncertainty - corine_uncertainty
    error_diff <- natura_error - corine_error
    median_relative_diff <- (natura_median - corine_median) / corine_median
    iqr_relative_diff <- (natura_iqr - corine_iqr) / corine_iqr
    iqr_median_relative_diff <- (natura_iqr_median - corine_iqr_median) / corine_iqr_median
    uncertainty_relative_diff <- (natura_uncertainty - corine_uncertainty) / corine_uncertainty
    error_relative_diff <- (natura_error - corine_error) / corine_error
    x <- data[1,]$x
    y <- data[1,]$y
    transect_postdis_dataframe[transect,] <- c(corine_median, natura_median,
                                               corine_iqr, natura_iqr,
                                               corine_iqr_median, natura_iqr_median,
                                               corine_uncertainty, natura_uncertainty,
                                               corine_error, natura_error,
                                               median_diff, iqr_diff, iqr_median_diff, uncertainty_diff, error_diff,
                                               median_relative_diff, iqr_relative_diff, iqr_median_relative_diff, uncertainty_relative_diff, error_relative_diff,
                                               x, y)
    
}

transect_postdis_dataframe$transect <- rownames(transect_postdis_dataframe)


long_df <- reshape(transect_postdis_dataframe,
                   varying = names(transect_postdis_dataframe)[names(transect_postdis_dataframe) != "transect"],
                   v.names = "value",
                   timevar = "metric",
                   times = names(transect_postdis_dataframe)[names(transect_postdis_dataframe) != "transect"],
                   idvar = "transect",
                   direction = "long")





median_natura_plot <- ggplot(transect_postdis_dataframe, aes(x, y)) +
    geom_point(size = 1 + 10 * transect_postdis_dataframe$natura_median, 
               color = "forestgreen")+
    coord_equal() +
    labs(title = "Median (Habitat type)") +
    theme_minimal()
iqr_natura_plot <- ggplot(transect_postdis_dataframe, aes(x, y)) +
    geom_point(size = 1 + 10 * transect_postdis_dataframe$natura_iqr, 
               color = "forestgreen") +
    coord_equal() +
    labs(title = "IQR (Habitat type)") +
    theme_minimal()
iqr_median_natura_plot <- ggplot(transect_postdis_dataframe, aes(x, y)) +
    geom_point(size = 1 + 0.5 * (transect_postdis_dataframe$natura_iqr_median), 
               color = "forestgreen") +
    coord_equal() +
    labs(title = "IQR/Median (Habitat type)") +
    theme_minimal()
uncertainty_natura_plot <- ggplot(transect_postdis_dataframe, aes(x, y)) +
    geom_point(size = 1 + 10 * transect_postdis_dataframe$natura_uncertainty, 
               color = "forestgreen") +
    coord_equal() +
    labs(title = "Uncertainty (Habitat type)") +
    theme_minimal()
error_natura_plot <- ggplot(transect_postdis_dataframe, aes(x, y)) +
    geom_point(size = 1 + 15 * transect_postdis_dataframe$corine_error, 
               color = "forestgreen") +
    coord_equal() +
    labs(title = "Error (Habitat type)") +
    theme_minimal()



grid.arrange(median_natura_plot, 
             iqr_natura_plot,
             iqr_median_natura_plot,
             uncertainty_natura_plot,
             error_natura_plot,
             ncol = 5)



create_transect_difference_plot <- function(transect_postdis_dataframe,
                                            type,
                                            title,
                                            relative_diff) {
    relative_string <- ""
    if (relative_diff == TRUE) {
        title <- sprintf("%s RELATIVE", title)
        relative_string <- "_relative"
    }
    
    title <- sprintf("%s DIFFERENCE", title)
    
    colname <- sprintf("%s%s_diff", type, relative_string)
    values <- transect_postdis_dataframe[,colname]
    
    plot <- ggplot(transect_postdis_dataframe, 
                   aes(x = x, 
                       y = y, 
                       color = factor(!!sym(colname) > 0, levels = c(FALSE, TRUE)))) +
            geom_point(size = 2 + 3 * abs((values / max(values)))) +
            scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                               labels = c("Land cover", "Habitat type"),
                               name = sprintf("Model with higher %s", type),
                               drop = FALSE) +
            labs(title = sprintf("%s\n\n%s/%s negative vs. positive differences | min %s | max %s",
                                 title,
                                 sum(values < 0),
                                 sum(values > 0),
                                 round(min(values), digits = 3),
                                 round(max(values), digits = 3))) +
            theme_minimal() +
            theme(plot.title = element_text(size = 10))
            coord_equal() 
            
    return(plot)
}


median_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "median", "MEDIAN", FALSE)
iqr_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "iqr", "IQR", FALSE)
iqr_median_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "iqr_median", "IQR / MEDIAN", FALSE)
uncertainty_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "uncertainty", "UNCERTAINTY", FALSE)
error_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "error", "ERROR", FALSE)
median_relative_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "median", "MEDIAN", TRUE)
iqr_relative_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "iqr", "IQR", TRUE)
iqr_median_relative_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "iqr_median", "IQR / MEDIAN", TRUE)
uncertainty_relative_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "uncertainty", "UNCERTAINTY", TRUE)
error_relative_difference_plot <- create_transect_difference_plot(transect_postdis_dataframe, "error", "ERROR", TRUE)


grid.arrange(median_relative_difference_plot,
             iqr_relative_difference_plot,
             iqr_median_relative_difference_plot,
             uncertainty_relative_difference_plot,
             error_relative_difference_plot,
             ncol = 3)





















#### OLD ANALYSIS #####



# Plot predictions against true values (an average over 2000 prediction samples)
for (species in colnames(occurrence)) {
    
    # First get true values and predictions
    data_for_species <- prediction_dataframe[prediction_dataframe$species == species,]
    true_values <- data_for_species$true_value
    predictions_corine <- data_for_species$corine_estimate
    predictions_natura <- data_for_species$natura_estimate
    
    prediction_info <- data.frame(x = rep(data_for_species$x, 3),
                                  y = rep(data_for_species$y, 3),
                                  year = rep(data_for_species$year, 3),
                                  Data = factor(rep(c("True", "Corine", "Natura"),
                                                    each = length(true_values)),
                                                levels = c("True", "Corine", "Natura")),
                                  Value = c(true_values,
                                            predictions_corine,
                                            predictions_natura))
    
    spatial_df <- aggregate(Value ~ x + y + Data,
                            data = prediction_info,
                            FUN = mean,
                            na.rm = TRUE)
    
    spatial_plot <- ggplot(spatial_df, aes(x, y, color = Value)) +
        geom_point(size = 1.5) +
        facet_wrap(~ Data, nrow = 1) +
        scale_color_gradient(low = "red",
                             high = "blue",
                             limits = c(0, 1)) +
        coord_equal() +
        labs(title = paste(species, "- Transect averages over years"),
            color = "Mean occurrence") +
        theme_minimal()
    

    temporal_df <- aggregate(Value ~ year + Data,
        data = prediction_info,
        FUN = mean,
        na.rm = TRUE)
    
    temporal_plot <- ggplot(temporal_df,
                         aes(year, Value, color = Data)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2) +
        scale_y_continuous(limits = c(0, 1)) +
        labs(title = paste(species, "- Yearly averages over transects"),
            y = "Mean occurrence",
            x = "Year",
            color = "Data") +
        theme_minimal()
    
    print(grid.arrange(spatial_plot, temporal_plot))
    
}


for (variable in c(natura_habitat_variables)) {
    if (variable == "NaturaPatchDensity") {
        prediction_dataframe[,variable] <- rep(natura_units[,"PatchDensity"], number_of_species)
    } else {
        prediction_dataframe[,variable] <- rep(natura_units[,variable], number_of_species)
    }
}

for (variable in c(corine_habitat_variables, other_variables)) {
    if (variable == "CorinePatchDensity") {
        prediction_dataframe[,variable] <- rep(corine_units[,"PatchDensity"], number_of_species)
    } else {
        prediction_dataframe[,variable] <- rep(corine_units[,variable], number_of_species)
    }
}


prediction_dataframe$estimate_difference <- prediction_dataframe$natura_estimate - prediction_dataframe$corine_estimate

prediction_dataframe$corine_error <- prediction_dataframe$corine_estimate - prediction_dataframe$true_value
prediction_dataframe$natura_error <- prediction_dataframe$natura_estimate - prediction_dataframe$true_value

prediction_dataframe$abs_error_difference <- abs(prediction_dataframe$natura_error) - abs(prediction_dataframe$corine_error)




plot_spatial_average <- function(data, column_name) {
    
    spatial_df <- aggregate(as.formula(paste(column_name, "~ x + y")),
                            data = data,
                            FUN = mean,
                            na.rm = TRUE)
    
    spatial_plot <- ggplot(spatial_df, aes(x, y, color = factor(.data[[column_name]] > 0, levels = c(FALSE, TRUE)))) +
        geom_point(size = 2 + 10 * abs(spatial_df[,column_name])) +
        scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red"),
                           labels = c("Negative", "Positive"),
                           name = "Sign",
                           drop = FALSE) +
        coord_equal() +
        labs(title = sprintf("%s", column_name),
             color = "Value") +
        theme_minimal()
    
    print(spatial_plot)
    
}

plot_temporal_average <- function(data, columns) {
    
    temporal_df <- aggregate(data[columns], 
                             by = list(year = data$year),
                             FUN = mean,
                             na.rm = TRUE)
    
    temporal_long <- reshape(temporal_df,
                             varying = columns,
                             v.names = "value",
                             timevar = "variable",
                             times = columns,
                             direction = "long")
    
    temporal_plot <- ggplot(temporal_long, aes(x = year, y = value, color = variable)) +
                        geom_line(linewidth = 1) +
                        geom_point(size = 2) +
                        geom_hline(yintercept = 0) +
                        labs(title = "Temporal averages", 
                             x = "Year",
                             y = "Mean value",
                             color = "Variable") +
                        theme_minimal()
    
    print(temporal_plot)

}

difference_names <- c("estimate_difference", "corine_error", "natura_error", "abs_error_difference")


for (name in difference_names) {
    plot_spatial_average(prediction_averages_over_species, name)
}

plot_temporal_average(prediction_averages_over_species, difference_names)



plot_against_environmental <- function(data, column_name) {
    
    old_par <- par(no.readonly = TRUE) 
    par(mfrow = c(3, 5), 
        mai = c(0.55, 0.3, 0.2, 0.01))
    for (variable in c(natura_habitat_variables)) {
        plot(prediction_dataframe[,variable], 
             prediction_dataframe[,column_name],
             xlab = variable,
             ylab = "",
             main = column_name,
             col = "blue")
        abline(0, 0)
        
    }
    for (variable in c(corine_habitat_variables, other_variables)) {
        plot(prediction_dataframe[,variable], 
             prediction_dataframe[,column_name],
             xlab = variable,
             ylab = "",
             main = column_name,
             col = "blue")
        abline(0, 0)
        
    }
    par(old_par)
    
}

for (name in difference_names) {
    plot_against_environmental(prediction_dataframe, name)
}


# PREDICTION CORRELATION AGAINST SPECIES TRAITS

pearson_correlation_species <- data.frame(species = unique(prediction_dataframe$species),
                                           correlation = 0)


for (species in unique(prediction_dataframe$species)) {
    data_for_species <- prediction_dataframe[prediction_dataframe$species == species,]
    correlation <- cor(data_for_species$corine_estimate, 
                       data_for_species$natura_estimate)
    pearson_correlation_species[pearson_correlation_species == species,]$correlation <- correlation
}

species_average_error_difference <- aggregate(abs_error_difference ~ species,
                                              data = prediction_dataframe,
                                              FUN = mean)
rownames(species_average_error_difference) <- species_average_error_difference$species
pearson_correlation_species$average_abs_error_difference <- species_average_error_difference[pearson_correlation_species$species,]$abs_error_difference



ggplot(pearson_correlation_species, 
       aes(x=reorder(species, correlation), 
           y=correlation,
           fill=factor(average_abs_error_difference > 0, 
                        levels = c(FALSE, TRUE)) )) + 
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "red"),
                       labels = c("TRUE" = "Corine",
                                  "FALSE" = "Natura"),
                       name = "Better model \n(smaller error)",
                       drop = FALSE) +
    geom_bar(stat = "identity") +
    labs(title = sprintf("Species correlation for predicted estimates \nmean %.3f, sd %.3f", 
                         mean(pearson_correlation_species$correlation),
                         sd(pearson_correlation_species$correlation)))




# load trait_data
load(file = file.path(dir_data, "trait_data.RData"))
load(file.path(dir_data, "species_prevalences.RData"))

pearson_correlation_species$feeding <- trait_data[trait_data$Species == pearson_correlation_species$species,]$Feeding
pearson_correlation_species$mig <- trait_data[trait_data$Species == pearson_correlation_species$species,]$Mig
pearson_correlation_species$mass <- trait_data[trait_data$Species == pearson_correlation_species$species,]$Mass
pearson_correlation_species$order <- trait_data[trait_data$Species == pearson_correlation_species$species,]$Order
pearson_correlation_species$sample_prevalence <- species_prevalences[pearson_correlation_species$species,]$samples
pearson_correlation_species$transect_prevalence <- species_prevalences[pearson_correlation_species$species,]$transects
pearson_correlation_species$year_prevalence <- species_prevalences[pearson_correlation_species$species,]$years




par(mfrow = c(1, 2))
boxplot(correlation ~ feeding, data = pearson_correlation_species,
        main = "Feeding type")
boxplot(correlation ~ mig, data = pearson_correlation_species,
        main = "Migration type")

par(mfrow = c(1,3))
for (trait in c("mass", "year_prevalence", "transect_prevalence")) {
    plot(pearson_correlation_species[,trait],
         pearson_correlation_species$correlation,
         main = trait,
         xlab = trait)
}



# PREDICTION CORRELATION AGAINST TRANSECT

unique_transects <- spatiotemporal_context[!duplicated(spatiotemporal_context$Transect), ]
pearson_correlation_transect <- data.frame(transect = unique_transects$Transect,
                                           x = unique_transects$x,
                                           y = unique_transects$y,
                                           correlation = 0
)
for (transect in unique(prediction_dataframe$transect)) {
    data_for_transect <- prediction_dataframe[prediction_dataframe$transect == transect,]
    correlation <- cor(data_for_transect$corine_estimate, 
                       data_for_transect$natura_estimate)
    pearson_correlation_transect[pearson_correlation_transect == transect,]$correlation <- correlation
}


ggplot(pearson_correlation_transect, 
       aes(x=reorder(transect, correlation), 
           y=correlation)) + 
    coord_flip() +
    geom_bar(stat = "identity", fill = "dodgerblue") +
    geom_text(aes(label = round(correlation, 3)), hjust = -0.2, size = 3) +
    labs(title = sprintf("Transect correlation for predicted estimates \nmean %.3f, sd %.3f", 
                         mean(pearson_correlation_transect$correlation),
                         sd(pearson_correlation_transect$correlation)))

transect_correlation_plot <- ggplot(pearson_correlation_transect, 
       aes(x, y, color = correlation)) +
    geom_point(size = 3 * pearson_correlation_transect$correlation) +
    scale_color_gradient(low = "red",
                         high = "blue") +
    coord_equal() +
    labs(title = sprintf("Transect correlation for predicted estimates \nmean %.3f, sd %.3f", 
                         mean(pearson_correlation_transect$correlation),
                         sd(pearson_correlation_transect$correlation))) +
    theme_minimal()


transect_average_error_difference <- aggregate(abs_error_difference ~ transect,
                                               data = prediction_dataframe,
                                               FUN = mean)
rownames(transect_average_error_difference) <- transect_average_error_difference$transect
pearson_correlation_transect$average_abs_error_difference <- transect_average_error_difference[pearson_correlation_transect$transect,]$abs_error_difference





transect_error_plot <- ggplot(pearson_correlation_transect, 
                       aes(x, 
                           y, 
                           color = factor(average_abs_error_difference > 0, 
                                          levels = c(FALSE, TRUE)))) +
    geom_point(size = 2 + 10 * abs(pearson_correlation_transect[,"average_abs_error_difference"])) +
    scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"),
                       labels = c("TRUE" = "Corine",
                                  "FALSE" = "Natura"),
                       name = "Better model \n(smaller error)",
                       drop = FALSE) +
    coord_equal() +
    labs(title = "Average abs error difference Natura - Corine",
         color = "Value") +
    theme_minimal()




grid.arrange(transect_correlation_plot, 
             transect_error_plot,
             nrow = 1)



pearson_correlation_transect$separator <- "High correlation"
pearson_correlation_transect[pearson_correlation_transect$correlation < 0.97,]$separator <- "Low correlation"
natura_better <- pearson_correlation_transect$average_abs_error_difference <= 0 & pearson_correlation_transect$separator == "Low correlation"
corine_better <- pearson_correlation_transect$average_abs_error_difference > 0 & pearson_correlation_transect$separator == "Low correlation"
pearson_correlation_transect[natura_better,]$separator <- "Natura better"
pearson_correlation_transect[corine_better,]$separator <- "Corine better"



correlation_long <- reshape(pearson_correlation_transect, 
                            direction = "long",
                            varying = c(natura_habitat_variables, 
                                        corine_habitat_variables, 
                                        other_variables),
                            v.names = "value",
                            timevar = "variable",
                            times = c(natura_habitat_variables, 
                                      corine_habitat_variables, 
                                      other_variables))

not_habitat_types <- c("Rainfall", "Temperature", "CorinePatchDensity", "NaturaPatchDensity")

ggplot(correlation_long[!correlation_long$variable %in% not_habitat_types,], aes(x = variable, y = value, fill = separator)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() + 
    labs(title = "Environmental variable comparison",
         x = "Environmental variable",
         y = "Amount",
         fill = "Separator"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("High correlation" = "purple", 
                                 "Natura better" = "red",
                                 "Corine better" = "blue"))







# PREDICTION CORRELATION AGAINST YEAR

pearson_correlation_year <- data.frame(year = unique(prediction_dataframe$year),
                                       correlation = 0)

for (year in unique(prediction_dataframe$year)) {
    data_for_year <- prediction_dataframe[prediction_dataframe$year == year,]
    correlation <- cor(data_for_year$corine_estimate, 
                       data_for_year$natura_estimate)
    pearson_correlation_year[pearson_correlation_year == year,]$correlation <- correlation
}

year_correlation_plot <- ggplot(pearson_correlation_year, 
       aes(x = year, y = correlation)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    labs(title = sprintf("Year correlation for predicted estimates \nmean %.3f, sd %.3f", 
                         mean(pearson_correlation_year$correlation),
                         sd(pearson_correlation_year$correlation))) +
    theme_minimal()


grid.arrange(transect_correlation_plot, 
             year_correlation_plot,
             nrow = 1)














# PREDICTION CORRELATION AGAINST ENVIRONMENTAL VARIABLES



old_par <- par(no.readonly = TRUE) 
par(mfrow = c(3, 5), 
    mai = c(0.55, 0.3, 0.2, 0.01))
for (variable in c(natura_habitat_variables, corine_habitat_variables, other_variables)) {
    
    pearson_correlation_transect[,variable] <- 0
    for (transect in pearson_correlation_transect$transect) {
        pearson_correlation_transect[pearson_correlation_transect == transect,
                                     variable] <- mean(prediction_dataframe[prediction_dataframe$transect == transect, 
                                                                                                           variable])
    }
    
    plot(pearson_correlation_transect[,variable],
         pearson_correlation_transect$correlation,
         main = variable,
         xlab = variable,
         ylab = "correlation")
    
}
par <- old_par

ggplot(pearson_correlation_env, aes(x = group, y = correlation, color = variable)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~variable)
    labs(title = "Correlations for environmental variables", 
         x = "Group",
         y = "Correlation",
         color = "Variable") +
    theme_minimal()
    


















