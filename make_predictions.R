generate_spatial_predictions <- TRUE

fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)


for (model_number in 1:length(fitted_models)) {
        
        # Load fitted models into "models"
        load(fitted_models[model_number])
        model_name <- strsplit(basename(fitted_models[model_number]), "\\.")[[1]][1]
        
        
        # Open up pdf to save prediction results
        pdf(file = file.path(dir_results, sprintf("%s_predictions.pdf", model_name)))
        

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
        number_of_prediction_units <- 100
        prediction_unit_rows <- sample(1:nrow(unit_data), number_of_prediction_units)
        prediction_units <- unit_data[prediction_unit_rows, ]
        
        # Set random noise to coordinates because 
        # predicting to the same exact coordinates does not work (see notes)
        new_coords <- data.frame(
            x = prediction_units$x + rnorm(nrow(prediction_units), mean = 0, sd = 0.01),
            y = prediction_units$y + rnorm(nrow(prediction_units), mean = 0, sd = 0.01)
        )
        
        # Similarly set random noise to years
        new_years <- data.frame(year = prediction_units$Year + rnorm(nrow(prediction_units), 
                                                                     mean = 0, 
                                                                     sd = 0.01))
        # Generate sample numbers
        new_sample <- data.frame(sample = 1:number_of_prediction_units)
        
        # Set the cluster values for prediction units
        new_env_vars <- prediction_units[,colnames(env_vars), drop = FALSE]
        
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
        
        # The prediction is a list of 1000 samples (250 from each chain)
        # Each sample is matrix of occurrence probabilities 
        # and has dimensions sampling units x species
        
        # Calculate the expected occurrence probability 
        # for each sampling unit / species -combination
        # Calculate averages by summing over all corresponding values in each matrix
        # Then dividing by the number of matrices
        expected_values <- Reduce("+", prediction) / length(prediction)
        save(expected_values, file = file.path(dir_results, "predictions.RData"))
            
        
        load(file.path(dir_results, "predictions.RData"))
        
        # Sum the expected occurrence probability for each species over each unit
        expected_species_richness <- rowSums(expected_values)
        
        # Community weighted mean traits
        # First get the expected occurrence probabilities (unit x species)
        # And the species traits (species x traits)
        # Now, for each sampling unit:
        # Take the species occurrence probabilities for each species
        # Take one trait from the trait matrix
        # Multiply the trait value with each of the species occurrence probabilities
        # Sum the results for that trait
        # Do this for all traits and all sampling units
        # Final matrix is community weighted trait values for each unit (units x traits)
        # Finally scale them all to same scale by dividing with species richness
        cwm <- (expected_values %*% model$Tr) /  matrix(rep(expected_species_richness, model$nt),
                                                        ncol = model$nt)
        
        # Example species (Fringilla montifringilla = Järripeippo)
        example_species <- expected_values[,"Fringilla_montifringilla"]
        
        temperature <- new_env_vars$Temperature
        
        prediction_data <- data.frame(new_coords, 
                                      new_years,
                                      temperature,
                                      expected_species_richness,
                                      example_species,
                                      cwm,
                                      stringsAsFactors = TRUE)
        
        # Plot clusters on map
        print(ggplot(data = prediction_data, 
                     aes(x = x, 
                         y = y, 
                         color = new_env_vars$Temperature)) + 
                  geom_point(size = 2) + 
                  ggtitle("Temperature") + 
                  coord_equal())
        
        # Plot prediction for example species
        print(ggplot(data = prediction_data, 
                     aes(x = x, 
                         y = y, 
                         color = example_species)) + 
                  geom_point(size = 2) + 
                  ggtitle(expression(italic("Fringilla montifringilla"))) + 
                  scale_color_gradient(low = "blue", high = "red") + 
                  coord_equal())
        
        # Plot predicted species richness
        print(ggplot(data = prediction_data, 
                     aes(x = x, 
                         y = y, 
                         color = expected_species_richness)) + 
                  geom_point(size = 2) + 
                  ggtitle("Species richness") + 
                  scale_color_gradient(low = "blue", high = "red") + 
                  coord_equal())
        
            
        
        
        dev.off()
        
    
}
