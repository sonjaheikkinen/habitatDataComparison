# LOAD DATA
load(file = file.path(dir_models, "models_unfitted.RData"))

overwrite <- FALSE


for (i in 1:length(model_list)) {
    
    model_to_fit <- model_list[[i]]
    model_name <- names(model_list[i])
    
    print(sprintf("Fitting model: %s", model_name))
    print(model_to_fit)
    print(sprintf("Fitting started %s", date()))
    
    
    for (j in 1:length(thinning_values)) {
        
        thinning_value <- thinning_values[j]
        model_file <- file.path(dir_fitted, sprintf("%s_thin_%s_fitted.Rdata", 
                                                        model_name,
                                                        thinning_value))
        found <- file.exists(model_file)
        
        if (found & !overwrite) {
            print(sprintf("Model %s with thinning value %s already fitted, not overwriting",
                          model_name,
                          thinning_value))
        } else {
            print(sprintf("Thinning value: %s, started %s", thinning_value, date()))
            
            fitted_model <- sampleMcmc(model_to_fit, 
                                       samples = number_of_samples, 
                                       thin = thinning_value, 
                                       transient = round(0.2 * number_of_samples * thinning_value),
                                       nChains = 4, 
                                       #nParallel = 4, #If parallel is set, verbose prints do not show
                                       verbose = 10)
            
            print(sprintf("Thinning value: %s, ended %s", thinning_value, date()))
            
            save(fitted_model, file = model_file)
        }
    }
    
    print(sprintf("Fitting ended %s", date()))
    print("")
    
}
