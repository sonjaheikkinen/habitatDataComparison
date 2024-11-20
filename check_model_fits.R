# SCRIPT FOR CALCULATING FIT VALUES FOR EACH MODEL

# FUNCTIONS FOR THE SCRIPT

calculate_modelfits <- function(folds, file, model, model_name) {
    
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
    explanatory_power <- as.data.frame(explanatory_power)
    rownames(explanatory_power) <- colnames(model$Y)
    
    
    # CALCULATE PREDICTIVE POWER
    
    # Partition the data to training and test sets
    # use transect column as the partition column
    # this ensures, that when predicting based on the habitat types of transect, 
    # the habitat types and occurrences for that transect have not yet been seen
    partition <- createPartition(model, nfolds = folds, column = "Transect")
    # Compute predicted values for test set, based on data fitted with training set
    predicted_values_test_set <- pcomputePredictedValues(model, 
                                                         partition = partition,
                                                         nParallel = 4)
    # Compute measured of model fit based on the predicted values for test set
    predictive_power <- evaluateModelFit(hM = model, predY = predicted_values_test_set)
    
    
    print(sprintf("Calculation ended %s", date()))
    print("")
    
    # SAVE RESULTS TO FILE
    save(explanatory_power, predictive_power, file = file)
    
}


# SCRIPT STARTS

# List filenames of fitted models
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

# Will old fits be overwritten?

for (model_number in 1:length(fitted_models)) {

    # GET MODEL INFORMATION
    load(fitted_models[model_number]) # Load file into fitted_model
    model <- fitted_model
    model_name <- strsplit(basename(fitted_models[model_number]), "\\.")[[1]][1]
    thinning_value <- strsplit(model_name, "_")[[1]][4]
    modelfit_file <- file.path(dir_modelfits, sprintf("modelfit_%s", model_name))
    
    # CALCULATE FIT IF NEEDED
    if (file.exists(modelfit_file) & !overwrite_modelfits) {
        print(sprintf("Modelfits already calculated for %s", model_name))
    } else {
        calculate_modelfits(as.numeric(thinning_value), modelfit_folds, modelfit_file, model, model_name)
    }
    
    # CREATE PLOTS OF FIT VALUES
    
    
}