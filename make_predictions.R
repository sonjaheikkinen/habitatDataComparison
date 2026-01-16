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





fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

expected_values_over_transects <- list()
expected_values_over_years <- list()
expected_values_over_variables <- list()
expected_values_over_variables_non_marginal <- list()
variable_scales <- list()
prediction_units_list <- list()

for (model_number in 1:length(fitted_models)) {
        
        # Load fitted models into "models"
        load(fitted_models[model_number])
        model_name <- strsplit(basename(fitted_models[model_number]), "\\.")[[1]][1]
        model <- fitted_model
        

        # Load original unit data
        env_vars <- data.frame(model$XData, check.names = FALSE)
        load(file = file.path(dir_data, "spatiotemporal_context.RData"))
        coords <- spatiotemporal_context
        unit_data <- env_vars
        unit_data$x <- coords$x
        unit_data$y <- coords$y
        unit_data$Year <- coords$Year
        unit_data$Transect <- coords$Transect
        
        # Randomly select some units for prediction
        #number_of_prediction_units <- 100
        #prediction_unit_rows <- sample(1:nrow(unit_data), number_of_prediction_units)
        #prediction_units <- unit_data[prediction_unit_rows, ]
        
        
        # PREDICT ALL TRANSECTS FOR ONE YEAR
        
        # First find the newest year with most transect visits
        transect_visits_per_year <- table(unit_data$Year)
        maximum_number_of_visits <- max(transect_visits_per_year)
        years_with_maximum_number_of_visits <- as.numeric(names(transect_visits_per_year[transect_visits_per_year == maximum_number_of_visits]))
        prediction_year <- max(years_with_maximum_number_of_visits)
        
        # Select prediction units from this year
        prediction_units <- unit_data[unit_data$Year == prediction_year,]
        prediction_unit_list_name <- sprintf("over_transects_%s", model_name)
        prediction_units_list[[prediction_unit_list_name]] <- prediction_units
        
        # Predict over transects
        prediction_over_transects <- make_prediction(prediction_units, colnames(env_vars))
        prediction_over_transects <- lapply(prediction_over_transects, function(x) {
                                                                        rownames(x) <- prediction_units$Transect
                                                                        return (x)})
        
        save(prediction_over_transects, file = file.path(dir_results, "predictions_over_transects.RData"))
        
        
        # PREDICT ALL YEARS FOR ONE TRANSCET
        
        number_of_years_transect_was_visited <- table(unit_data$Transect)
        prediction_transect <- names(which.max(number_of_years_transect_was_visited))
        
        prediction_units <- unit_data[unit_data$Transect == prediction_transect,]
        prediction_unit_list_name <- sprintf("over_years_%s", model_name)
        prediction_units_list[[prediction_unit_list_name]] <- prediction_units
        prediction_over_years <- make_prediction(prediction_units, colnames(env_vars))
        prediction_over_years <- lapply(prediction_over_years, function(x) {
                                                                rownames(x) <- prediction_units$Year
                                                                return (x)
                                                                })
        
        save(prediction_over_years, file = file.path(dir_results, "predictions_over_transects.RData"))
        
        
        
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
        
        
        
        
        # The prediction is a list of 2000 samples (500 from each chain)
        # Each sample is matrix of occurrence probabilities 
        # and has dimensions sampling units x species
        
        # CALCULATE EXPECTED OCCURRENCE PROBABILITIES FOR EACH SPECIES AT EACH SAMPLING UNIT
        

        # Calculate averages by summing over all corresponding values in each matrix, then dividing by the number of matrices
        expected_values <- Reduce("+", prediction_over_transects) / length(prediction_over_transects)
        expected_values_over_transects[[model_name]] <- expected_values
        
        expected_values <- Reduce("+", prediction_over_years) / length(prediction_over_years)
        expected_values_over_years[[model_name]] <- expected_values
        
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



load(fitted_models[[1]])
occurrences <- fitted_model$Y

natura_habitat_variables <- c("Luonnonmetsät",
                              "Tunturikoivikot",
                              "Lehdot",
                              "Tulvametsät",
                              "PatchDensity")
corine_habitat_variables <- c("Havumetsät.kivennäismaalla",
                              "Sekametsät.kivennäismaalla",
                              "Sekametsät.turvemaalla",
                              "Lehtimetsät.kivennäismaalla",
                              "Havumetsät.kalliomaalla",
                              "PatchDensity")
other_variables <- c("Temperature", "Rainfall")

# LOAD SPATIOTEMPORAL CONTEXT
spatiotemporal_context <- model$studyDesign


# PREDICTIONS OVER YEARS
rows_for_prediction_transect <- spatiotemporal_context$Transect == prediction_transect

for (species in colnames(occurrences)) {
    true_values <- occurrences[rows_for_prediction_transect,species]
    predictions_corine <- expected_values_over_years[[1]][,species]
    predictions_natura <- expected_values_over_years[[2]][,species]
    
    ymin <- min(c(true_values, predictions_corine, predictions_natura))
    ymax <- max(c(true_values, predictions_corine, predictions_natura), ymin + 0.05)
    
    plot(spatiotemporal_context[rows_for_prediction_transect,]$Year, 
         true_values,
         ylim = c(0, 1),
         col = "black",
         main = species)
    points(spatiotemporal_context[rows_for_prediction_transect,]$Year, 
           predictions_corine,
           col = "blue")
    points(spatiotemporal_context[rows_for_prediction_transect,]$Year, 
           predictions_natura,
           col = "red")
}


# PREDICTIONS OVER TRANSECTS
rows_for_prediction_year <- spatiotemporal_context$Year == prediction_year
transect_problems <- data.frame(transect = spatiotemporal_context[rows_for_prediction_year,]$Transect,
                                problems_natura = 0,
                                problems_corine = 0)
rownames(transect_problems) <- transect_problems$transect


prediction_differences <- data.frame(transect = prediction_units_list$over_transects_3_probit_natura_forest_thin_100_samples_500_fitted$Transect,
                                     estimate_difference = 0,
                                     error_difference = 0)
rownames(prediction_differences) <- prediction_differences$transect
prediction_differences$transect <- NULL


for (species in colnames(occurrences)) {
    
    true_values <- occurrences[rows_for_prediction_year,species]
    names(true_values) <- spatiotemporal_context[rows_for_prediction_year,]$Transect
    predictions_corine <- expected_values_over_transects[[1]][,species]
    predictions_natura <- expected_values_over_transects[[2]][,species]
    for (transect in rownames(prediction_differences)) {
        estimate_difference <- predictions_natura[transect] - predictions_corine[transect]
        corine_abs_error <- abs(predictions_corine[transect] - true_values[transect])
        natura_abs_error <- abs(predictions_natura[transect] - true_values[transect])
        abs_error_difference <- natura_abs_error - corine_abs_error
        prediction_differences[transect,]$estimate_difference <- prediction_differences[transect,]$estimate_difference + estimate_difference
        prediction_differences[transect,]$error_difference <- prediction_differences[transect,]$error_difference + abs_error_difference
    }

    
    plot(spatiotemporal_context[rows_for_prediction_year,]$Transect, 
         true_values,
         ylim = c(0, 1),
         col = "black",
         main = species)
    points(spatiotemporal_context[rows_for_prediction_year,]$Transect, 
           predictions_corine,
           col = "blue")
    points(spatiotemporal_context[rows_for_prediction_year,]$Transect, 
           predictions_natura,
           col = "red")
}

# Calculate mean difference and mean absolute error difference
prediction_differences$estimate_difference <- prediction_differences$estimate_difference / nrow(prediction_differences)
prediction_differences$error_difference <- prediction_differences$error_difference / nrow(prediction_differences)

prediction_units_natura <- prediction_units_list$over_transects_3_probit_natura_forest_thin_100_samples_500_fitted
rownames(prediction_units_natura) <- prediction_units_natura$Transect
prediction_units_corine <- prediction_units_list$over_transects_3_probit_corine_forest_thin_100_samples_500_fitted
rownames(prediction_units_corine) <- prediction_units_corine$Transect

old_par <- par(no.readonly = TRUE) 
par(mfrow = c(3, 5), 
    mai = c(0.55, 0.3, 0.2, 0.01))
for (variable in c(natura_habitat_variables)) {
    prediction_differences[,variable] <- prediction_units_natura[rownames(prediction_differences),][,variable]
    plot(prediction_differences[,variable], prediction_differences$estimate_difference,
         xlab = variable,
         ylab = "",
         col = "blue")
    abline(0, 0)
    
}
for (variable in c(corine_habitat_variables, other_variables)) {
    prediction_differences[,variable] <- prediction_units_corine[rownames(prediction_differences),][,variable]
    plot(prediction_differences[,variable], prediction_differences$estimate_difference,
         xlab = variable,
         ylab = "",
         col = "blue")
    abline(0, 0)
    
}
par(old_par)

values <- prediction_differences$error_difference
limit <- max(abs(min(values)), abs(max(values)))

old_par <- par(no.readonly = TRUE) 
par(mfrow = c(3, 5), 
    mai = c(0.55, 0.3, 0.2, 0.01))
for (variable in c(natura_habitat_variables)) {
    prediction_differences[,variable] <- prediction_units_natura[rownames(prediction_differences),][,variable]
    plot(prediction_differences[,variable], prediction_differences$error_difference,
         xlab = variable,
         ylab = "",
         ylim = c(-limit, limit),
         col = "blue")
    abline(0, 0)
    
}
for (variable in c(corine_habitat_variables, other_variables)) {
    prediction_differences[,variable] <- prediction_units_corine[rownames(prediction_differences),][,variable]
    plot(prediction_differences[,variable], prediction_differences$error_difference,
         xlab = variable,
         ylab = "",
         ylim = c(-limit, limit),
         col = "blue")
    abline(0, 0)
    
}
par(old_par)


# PREDICTIONS OVER EXPECTED VALUES 

par(mfrow = c(2, 2))
comparable_variables <- c("Temperature", "Rainfall")
for (variable in comparable_variables) {
    
    prediction_for_corine <- expected_values_over_variables[[1]][[variable]]
    prediction_for_natura <- expected_values_over_variables[[2]][[variable]]
    prediction_for_corine_non_marginal <- expected_values_over_variables_non_marginal[[1]][[variable]]
    prediction_for_natura_non_marginal <- expected_values_over_variables_non_marginal[[2]][[variable]]
    colors <- rainbow(ncol(occurrences))
    
    plot(variable_scales[[1]][[variable]],
         prediction_for_corine[,1],
         col = colors[1],
         xlab = variable,
         ylab = "Predicted occurrence",
         ylim = c(0, 1),
         type = "l",
         main = sprintf("Corine %s", variable))
    for (column in 2:ncol(occurrences)) {
        lines(variable_scales[[1]][[variable]],
             prediction_for_corine[,column],
             col = colors[column])
    }
    mean_prediction <- rowMeans(prediction_for_corine, na.rm = TRUE)
    lines(variable_scales[[1]][[variable]],
          mean_prediction,
          col = "black",
          lwd = 3,          
          lty = 2) 
    
    plot(variable_scales[[2]][[variable]],
         prediction_for_natura[,1],
         col = colors[1],
         xlab = variable,
         ylab = "Predicted occurrence",
         ylim = c(0, 1),
         type = "l",
         main = sprintf("Natura %s", variable))
    for (column in 2:ncol(occurrences)) {
        lines(variable_scales[[2]][[variable]],
             prediction_for_natura[,column],
             col = colors[column])
    }
    mean_prediction <- rowMeans(prediction_for_natura, na.rm = TRUE)
    lines(variable_scales[[2]][[variable]],
          mean_prediction,
          col = "black",
          lwd = 3,          
          lty = 2)
    
    
    plot(variable_scales[[1]][[variable]],
         prediction_for_corine_non_marginal[,1],
         col = colors[1],
         xlab = variable,
         ylab = "Predicted occurrence",
         ylim = c(0, 1),
         type = "l",
         main = sprintf("Non marginal Corine %s", variable))
    for (column in 2:ncol(occurrences)) {
        lines(variable_scales[[1]][[variable]],
              prediction_for_corine_non_marginal[,column],
              col = colors[column])
    }
    mean_prediction <- rowMeans(prediction_for_corine_non_marginal, na.rm = TRUE)
    lines(variable_scales[[1]][[variable]],
          mean_prediction,
          col = "black",
          lwd = 3,          
          lty = 2)   
    
    plot(variable_scales[[2]][[variable]],
         prediction_for_natura_non_marginal[,1],
         col = colors[1],
         xlab = variable,
         ylab = "Predicted occurrence",
         ylim = c(0, 1),
         type = "l",
         main = sprintf("Non marginal Natura %s", variable))
    for (column in 2:ncol(occurrences)) {
        lines(variable_scales[[2]][[variable]],
              prediction_for_natura_non_marginal[,column],
              col = colors[column],
              main = sprintf("Natura %s", variable))
    }
    mean_prediction <- rowMeans(prediction_for_natura_non_marginal, na.rm = TRUE)
    lines(variable_scales[[2]][[variable]],
          mean_prediction,
          col = "black",
          lwd = 3,          
          lty = 2)
    
}

par(mfrow = c(1, 2))


for (variable in natura_habitat_variables) {
    
    prediction_for_natura <- expected_values_over_variables[[2]][[variable]]
    prediction_for_natura_non_marginal <- expected_values_over_variables_non_marginal[[2]][[variable]]
    variable_values <- variable_scales[[2]][[variable]]
    colors <- rainbow(ncol(occurrences))
    
    plot(variable_values,
         prediction_for_natura[,1],
         col = colors[1],
         xlab = variable,
         ylab = "Predicted occurrence",
         ylim = c(0, 1),
         xlim = c(min(variable_values), max(variable_values)),
         type = "l",
         main = sprintf("Natura %s", variable))
    for (column in 2:ncol(occurrences)) {
        lines(variable_values,
              prediction_for_natura[,column],
              col = colors[column])
    }
    mean_prediction <- rowMeans(prediction_for_natura, na.rm = TRUE)
    lines(variable_values,
          mean_prediction,
          col = "black",
          lwd = 3,          
          lty = 2)  
    
    
    plot(variable_values,
         prediction_for_natura_non_marginal[,1],
         col = colors[1],
         xlab = variable,
         ylab = "Predicted occurrence",
         ylim = c(0, 1),
         xlim = c(min(variable_values), max(variable_values)),
         type = "l",
         main = sprintf("Non marginal Natura %s", variable))
    for (column in 2:ncol(occurrences)) {
        lines(variable_values,
              prediction_for_natura_non_marginal[,column],
              col = colors[column])
    }
    mean_prediction <- rowMeans(prediction_for_natura_non_marginal, na.rm = TRUE)
    lines(variable_values,
          mean_prediction,
          col = "black",
          lwd = 3,          
          lty = 2)   
}


for (variable in corine_habitat_variables) {
    
    prediction_for_corine <- expected_values_over_variables[[1]][[variable]]
    prediction_for_corine_non_marginal <- expected_values_over_variables_non_marginal[[1]][[variable]]
    variable_values <- variable_scales[[1]][[variable]]
    colors <- rainbow(ncol(occurrences))
    
    plot(variable_values,
         prediction_for_corine[,1],
         col = colors[1],
         xlab = variable,
         ylab = "Predicted occurrence",
         ylim = c(0, 1),
         xlim = c(min(variable_values), max(variable_values)),
         type = "l",
         main = sprintf("Corine %s", variable))
    for (column in 2:ncol(occurrences)) {
        lines(variable_values,
              prediction_for_corine[,column],
              col = colors[column])
    }
    mean_prediction <- rowMeans(prediction_for_corine, na.rm = TRUE)
    lines(variable_values,
          mean_prediction,
          col = "black",
          lwd = 3,          
          lty = 2)     
    
    plot(variable_values,
         prediction_for_corine_non_marginal[,1],
         col = colors[1],
         xlab = variable,
         ylab = "Predicted occurrence",
         ylim = c(0, 1),
         xlim = c(min(variable_values), max(variable_values)),
         type = "l",
         main = sprintf("Non-marginal Corine %s", variable))
    for (column in 2:ncol(occurrences)) {
        lines(variable_values,
              prediction_for_corine_non_marginal[,column],
              col = colors[column])
    }
    mean_prediction <- rowMeans(prediction_for_corine_non_marginal, na.rm = TRUE)
    lines(variable_values,
          mean_prediction,
          col = "black",
          lwd = 3,          
          lty = 2)   
    
    
}
















