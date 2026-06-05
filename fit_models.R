# LOAD DATA
load(file = file.path(dir_models, "models_unfitted.RData"))

overwrite <- TRUE

sink(file.path(dir_results, "output.txt"), split = TRUE, append = TRUE)


for (model_number in 1:length(model_list)) {
    
    model_to_fit <- model_list[[model_number]]
    model_name <- names(model_list[model_number])
    
    print(sprintf("Fitting model: %s", model_name))
    print(model_to_fit)
    print(sprintf("Fitting started %s", date()))
    
    
    for (thinning_value_number in 1:length(thinning_values)) {
        
        thinning_value <- thinning_values[thinning_value_number]
        
        
        
        model_file <- file.path(dir_fitted, sprintf("3_%s_thin_%s_samples_%s_fitted.RData", 
                                                        model_name,
                                                        thinning_value,
                                                        number_of_samples))
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
                                       transient = ceiling(0.5 * number_of_samples * thinning_value),
                                       nChains = 4,
                                       nParallel = 4, 
                                       verbose = 10)
            
            print(sprintf("Thinning value: %s, ended %s", thinning_value, date()))
            
            save(fitted_model, file = model_file)
        }
    }
    
    print(sprintf("Fitting ended %s", date()))
    print("")
}

sink()
