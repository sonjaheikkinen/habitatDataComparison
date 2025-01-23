# SCRIPT FOR CALCULATING FIT VALUES FOR EACH MODEL

# FUNCTIONS FOR THE SCRIPT

extract_thinning_value <- function(filename) {
    parts <- strsplit(filename, "_")[[1]] 
    thinning_value <- as.numeric(parts[length(parts) - 1])
    print(thinning_value)
    return(thinning_value)
}

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
    # Set species names as rownames
    explanatory_power <- as.data.frame(explanatory_power)
    rownames(explanatory_power) <- colnames(model$Y)
    
    
    # CALCULATE PREDICTIVE POWER
    
    # Partition the data to training and test sets
    # use transect column as the partition column
    # this ensures, that when predicting based on the habitat types of transect, 
    # the habitat types and occurrences for that transect have not yet been seen
    partition <- createPartition(model, nfolds = folds, column = "Transect")
    # Compute predicted values for test set, based on data fitted with training set
    predicted_values_test_set <- computePredictedValues(model, 
                                                        partition = partition,
                                                        nParallel = 4)
    # Compute measures of model fit based on the predicted values for test set
    predictive_power <- evaluateModelFit(hM = model, predY = predicted_values_test_set)
    # Set species names as rownames
    predictive_power <- as.data.frame(predictive_power)
    rownames(predictive_power) <- colnames(model$Y)
    
    
    print(sprintf("Calculation ended %s", date()))
    print("")
    
    # SAVE RESULTS TO FILE
    save(explanatory_power, predictive_power, file = file)
    
}


create_modelfit_plot <- function(explanatory_power, predictive_power, type, model_name, thinning_value) {
    
    if (!type %in% colnames(explanatory_power)) {
        return("Type not yet calculated or does not apply for this model")
    }
    if (!is.null(explanatory_power[,type])) {
        plot(explanatory_power[,type],
             predictive_power[,type],
             xlim = c(-1,1),
             ylim = c(-1,1),
             xlab = "explanatory power (MF)",
             ylab = "predictive power (MFCV)",
             main = sprintf("%s\n thin = %s: %s. \nmean(MF) = %s, mean(MFCV) = %s",
                            model_name,
                            as.character(thinning_value),
                            type,
                            as.character(mean(explanatory_power[,type], na.rm = TRUE)),
                            as.character(mean(predictive_power[,type], na.rm = TRUE))))
        abline(h = 0)
        abline(v = 0)
    }
    
}


# SCRIPT STARTS

# List filenames of fitted models
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

fitted_models <- fitted_models[7]

for (model_number in 1:length(fitted_models)) {

    # GET MODEL INFORMATION
    load(fitted_models[model_number]) # Load file into fitted_model
    model <- fitted_model
    model_name <- strsplit(basename(fitted_models[model_number]), "\\.")[[1]][1]
    thinning_value <- extract_thinning_value(model_name)
    modelfit_file <- file.path(dir_modelfits, sprintf("modelfit_%s", model_name))
    
    # CREATE PDF FOR PLOTTING THE FIT VALUES
    pdf(file = file.path(dir_results, sprintf("modelfit_results_%s.pdf", model_name)))
    
    # CALCULATE FIT IF NEEDED
    if (file.exists(modelfit_file) & !overwrite_modelfits) {
        print(sprintf("Modelfits already calculated for %s", model_name))
    } else {
        calculate_modelfits(modelfit_folds, modelfit_file, model, model_name)
    }
    
    # CREATE PLOTS OF FIT VALUES
    if (file.exists(modelfit_file)) {
        load(modelfit_file)
        create_modelfit_plot(explanatory_power, predictive_power, "TjurR2", model_name, thinning_value)
        create_modelfit_plot(explanatory_power, predictive_power, "R2", model_name, thinning_value)
        create_modelfit_plot(explanatory_power, predictive_power, "AUC", model_name, thinning_value)
        create_modelfit_plot(explanatory_power, predictive_power, "SR2", model_name, thinning_value)
    } else {
        print(sprintf("Modelfit results not found for %s", model_name))
    }
    
}

# Close pdf
dev.off()
