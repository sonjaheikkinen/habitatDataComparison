make_prediction <- function(prediction_units, env_vars) {
    
    # Set random noise to coordinates because 
    # predicting to the same exact coordinates does not work 
    new_coords <- data.frame(
        x = prediction_units$x + rnorm(nrow(prediction_units), mean = 0, sd = 0.01),
        y = prediction_units$y + rnorm(nrow(prediction_units), mean = 0, sd = 0.01)
    )
    
    # Similarly set random noise to years
    new_years <- data.frame(year = prediction_units$Year + rnorm(nrow(prediction_units), 
                                                                 mean = 0, 
                                                                 sd = 0.01))
    # Generate sample numbers
    new_sample <- data.frame(sample = 1:nrow(prediction_units))
    
    # Get environmental variables
    new_env_vars <- prediction_units[,env_vars, drop = FALSE]
    
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





generate_spatial_predictions <- TRUE

fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

expected_values_over_transects <- list()
expected_values_over_years <- list()
prediction_transect <- c()
prediction_year <- c()


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
        
        # Predict over transects
        prediction_over_transects <- make_prediction(prediction_units, colnames(env_vars))
        
        save(prediction_over_transects, file = file.path(dir_results, "predictions_over_transects.RData"))
        
        
        # PREDICT ALL YEARS FOR ONE TRANSCET
        
        number_of_years_transect_was_visited <- table(unit_data$Transect)
        prediction_transect <- names(which.max(number_of_years_transect_was_visited))
        
        prediction_units <- unit_data[unit_data$Transect == prediction_transect,]
        prediction_over_years <- make_prediction(prediction_units, colnames(env_vars))
        
        save(prediction_over_years, file = file.path(dir_results, "predictions_over_transects.RData"))
        
        
        # The prediction is a list of 2000 samples (500 from each chain)
        # Each sample is matrix of occurrence probabilities 
        # and has dimensions sampling units x species
        
        # CALCULATE EXPECTED OCCURRENCE PROBABILITIES FOR EACH SPECIES AT EACH SAMPLING UNIT
        

        # Calculate averages by summing over all corresponding values in each matrix, then dividing by the number of matrices
        expected_values <- Reduce("+", prediction_over_transects) / length(prediction_over_transects)
        expected_values_over_transects[[model_name]] <- expected_values
        
        expected_values <- Reduce("+", prediction_over_years) / length(prediction_over_years)
        expected_values_over_years[[model_name]] <- expected_values

}



load(fitted_models[[1]])
occurrences <- fitted_model$Y

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

for (species in colnames(occurrences)) {
    true_values <- occurrences[rows_for_prediction_year,species]
    predictions_corine <- expected_values_over_transects[[1]][,species]
    predictions_natura <- expected_values_over_transects[[2]][,species]
    
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






