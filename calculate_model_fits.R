# SCRIPT FOR CALCULATING FIT VALUES FOR EACH MODEL

overwrite_modelfits <- TRUE


# FUNCTIONS FOR THE SCRIPT



calculate_modelfits <- function(folds, file, model, model_name, partition_transect, partition_year) {
    
    print(sprintf("Calculating fits for %s", model_name))
    print(sprintf("Calculation started %s", date()))
    
    # CALCULATE EXPLANATORY POWER
    
    # Compute predicted values for species occurrence/abundance (matrix Y)
    # based on model fitted using the whole data
    predicted_values_whole_data <- computePredictedValues(model)
    # Compute measures of model fit based on the predicted values
    # Each modelfit measure is a vector with one value for each species
    # The species are in same order as the columns in matrix Y (based on source code)
    explanatory_power <- evaluateModelFit(hM = model, predY = predicted_values_whole_data)
    # Set species names as rownames
    explanatory_power <- as.data.frame(explanatory_power)
    rownames(explanatory_power) <- colnames(model$Y)
    
    
    # CALCULATE PREDICTIVE POWER BY TRANSECT
    # Compute predicted values for test set, based on data fitted with training set
    #predicted_values_test_set <- computePredictedValues(model, 
    #                                                    partition = partition_transect,
    #                                                    nParallel = 4)
    # Compute measures of model fit based on the predicted values for test set
    #predictive_power_transect <- evaluateModelFit(hM = model, predY = predicted_values_test_set)
    # Set species names as rownames
    #predictive_power_transect <- as.data.frame(predictive_power_transect)
    #rownames(predictive_power_transect) <- colnames(model$Y)
    
    
    # CALCULATE PREDICTIVE POWER BY YEAR
    # Compute predicted values for test set, based on data fitted with training set
    predicted_values_test_set <- computePredictedValues(model, 
                                                        partition = partition_year,
                                                        nParallel = 4)
    # Compute measures of model fit based on the predicted values for test set
    predictive_power_year <- evaluateModelFit(hM = model, predY = predicted_values_test_set)
    # Set species names as rownames
    predictive_power_year <- as.data.frame(predictive_power_year)
    rownames(predictive_power_year) <- colnames(model$Y)
    
    
    
    # PREDICTIVE POWER BY  FOCAL TRANSECT
    # TO DO: Refactor so that focal transects and years are the same for every model
    #number_of_focal_transects <- 1
    #focal_transect_indices <- sample(1:length(unique(model$studyDesign$Transect)), number_of_focal_transects)
    #focal_transects <- unique(model$studyDesign$Transect)[focal_transect_indices]
    #partition <- (model$studyDesign$Transect == focal_transects[1]) * 1
    #partition[partition == 0] <- 2
    #predicted_values_test_set <- computePredictedValues(model,
    #                                                    partition = partition,
    #                                                    nParallel = 4)
    #predictive_power_transect <- evaluateModelFit(hM = model, predY = predicted_values_test_set)
    #predictive_power_transect <- as.data.frame(predictive_power_transect)
    
    
    
    
    # CALCULATE WAIC
    waic <-  computeWAIC(model)
    
    
    
    print(sprintf("Calculation ended %s", date()))
    print("")
    
    # SAVE RESULTS TO FILE
    # TO DO: Combine all in same files, then edit saving so that everything goes to one file
    save(explanatory_power, predictive_power_transect, waic, file = file)
    save(predictive_power_year, file = file.path(dir_modelfits, sprintf("modelfit_year_%s.RData", model_name)))
    
}



# SCRIPT STARTS

# List filenames of fitted models
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

partition_transect <- c()
partition_year <- c()

for (model_number in 1:length(fitted_models)) {

    # GET MODEL INFORMATION
    load(fitted_models[model_number]) # Load file into fitted_model
    model <- fitted_model
    model_name <- strsplit(basename(fitted_models[model_number]), "\\.")[[1]][1]
    modelfit_file <- file.path(dir_modelfits, sprintf("modelfit_%s.RData", model_name))
    
    # Partition the data to training and test sets
    # use transect column as the partition column
    # this ensures, that when predicting based on the habitat types of transect, 
    # the habitat types and occurrences for that transect have not yet been seen
    if (model_number <- 1) {
        partition_transect <- createPartition(model, nfolds = modelfit_folds, column = "Transect")
        partition_year <- createPartition(model, nfolds = modelfit_folds, column = "Year")
    }
    
    # CALCULATE FIT IF NEEDED
    if (file.exists(modelfit_file) & !overwrite_modelfits) {
        print(sprintf("Modelfits already calculated for %s", model_name))
    } else {
        calculate_modelfits(modelfit_folds, 
                            modelfit_file, 
                            model, 
                            model_name, 
                            partition_transect, 
                            partition_year)
    }
    
    
}



