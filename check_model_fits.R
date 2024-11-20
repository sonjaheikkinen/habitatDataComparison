# SCRIPT FOR CALCULATING FIT VALUES FOR EACH MODEL

# FUNCTIONS FOR THE SCRIPT

calculate_modelfits <- function(thinning_value, folds, file, model, model_name) {
    
    print(sprintf("Calculating fits for %s", model_name))
    print(sprintf("Calculation started %s", date()))
    
    
    print(sprintf("Calculation ended %s", date()))
    print("")
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
    
    # SAVE FIT VALUES TO FILE
    
    
}