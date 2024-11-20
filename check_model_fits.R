# SCRIPT FOR CALCULATING FIT VALUES FOR EACH MODEL

# FUNCTIONS FOR THE SCRIPT


# SCRIPT STARTS

# List filenames of fitted models
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

for (model_number in 1:length(fitted_models)) {

    # GET MODEL INFORMATION
    load(fitted_models[model_number]) # Load file into fitted_model
    model <- fitted_model
    model_name <- strsplit(basename(fitted_models[model_number]), "\\.")[[1]][1]
    thinning_value <- strsplit(model_name, "_")[[1]][4]
    modelfit_file <- file.path(dir_modelfits, sprintf("modelfit_%s", model_name))
    print(modelfit_file)
    
    # CALCULATE FIT IF NEEDED
    
    
    # SAVE FIT VALUES TO FILE
    
    
}