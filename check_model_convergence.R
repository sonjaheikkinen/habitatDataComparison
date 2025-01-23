# SCRIPT FOR CHECKING MODEL CONVERGENCE

# CONVERGENCE FUNCTIONS
plot_convergence <- function(variable, model_name, psrf_point_estimates) {
    par(mfrow=c(1,2))
    vioplot(psrf_point_estimates, 
            col = "dodgerblue",
            ylim = c(min(psrf_point_estimates) - (max(psrf_point_estimates) - min(psrf_point_estimates)), 
                     max(psrf_point_estimates)))
    vioplot(psrf_point_estimates,
            col = "dodgerblue",
            ylim = c(0.9,1.1))
    title(sprintf("PSRF (%s): %s", variable, model_name), 
          outer = TRUE,
          line = -3)
} 

save_convergence_info <- function(samples, variable, model_name, convergence_file) {
    psrf <- gelman.diag(samples, multivariate = FALSE)$psrf
    append_to_file(sprintf("%s\n", summary(psrf)[,1]), file = convergence_file)
    append_to_file("\n\n", file = convergence_file)
    psrf_point_estimates <- psrf[,1]
    plot_convergence(variable, model_name, psrf_point_estimates)
}

extract_thinning_value <- function(filename) {
    parts <- strsplit(filename, "_")[[1]] 
    thinning_value <- as.numeric(parts[length(parts) - 1])
    print(thinning_value)
    return(thinning_value)
}


# SCRIPT STARTS

# List filenames of fitted models
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)
model_thinning_values <- sapply(fitted_models, extract_thinning_value)
fitted_models <- fitted_models[order(model_thinning_values)]

# Create file for convergence plots
pdf(file = file.path(dir_results,"convergence_results.pdf"))

# Create text file for convergence results
convergence_file <- file.path(dir_results, "convergence_results.txt")

for (model_number in 1:length(fitted_models)) {
    
    # GET MODEL INFORMATION
    load(fitted_models[model_number]) # Load file into fitted_model
    model <- fitted_model
    model_name <- strsplit(basename(fitted_models[model_number]), "\\.")[[1]][1]
    append_to_file(sprintf("%s\n\n", model_name), convergence_file)
    thinning_value <- strsplit(model_name, "_")[[1]][4]
    number_of_random_levels <- model$nr
    
    print(sprintf("Calculating convergence for %s", model_name))
    
    # GET THE POSTERIOR SAMPLES OF PARAMETERS AS CODA OBJECT
    # Get only names (not numbers) of species and covariates
    posterior <- convertToCodaObject(model, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    print("Converted to coda object")
    
    # CALCULATE THE CONVERGENCE FOR EACH VARIABLE
    # First get posterior samples for each variable
    # The posterior sample object is a list of 4, containing the samples (250) for each of the 4 chains
    # The samples themselves are vectorized matrices, with variable dimensions depending on the variable
    # First calculate the potential scale reduction factors
    # Then write the point estimate summaries (min, max, mean, median, quantiles) to file
    # Finally, create column in the corresponding dataframe to save the actual estimates for plotting
    
    
    # CONVERGENCE FOR SPECIES NICHES (BETA)
    append_to_file("BETA\n\n", file = convergence_file)
    samples_beta <- posterior$Beta # 
    save_convergence_info(samples_beta, "beta", model_name, convergence_file)
    print("Calculated convergence for beta")
    
    # CONVERGENCE FOR INFLUENCE OF TRAITS (GAMMA)
    append_to_file("GAMMA\n\n", file = convergence_file)
    samples_gamma <- posterior$Gamma
    save_convergence_info(samples_gamma, "gamma", model_name, convergence_file)
    print("Calculated convergence for gamma")
    
    # CONVERGENCE FOR INFLUENCE OF PHYLOGENY (RHO)
    if (!is.null(posterior$Rho)) {
        append_to_file("RHO\n\n", file = convergence_file)
        samples_rho <- posterior$Rho
        save_convergence_info(samples_rho, "rho", model_name, convergence_file)
        print("Calculated convergence for rho")
    }
    
    if (number_of_random_levels > 0) {
        
        # CONVERGENCE FOR RANDOM VARIATION IN CO-OCCURRENCE (OMEGA)
        append_to_file("OMEGA\n\n", file = convergence_file)
        samples_omega <- posterior$Omega
        for (randomlevel_number in 1:number_of_random_levels) {
            randomlevel_name <- names(model$ranLevels)[randomlevel_number]
            append_to_file(sprintf("%s\n\n", randomlevel_name), file = convergence_file)
            samples_omega_randomlevel <- samples_omega[[randomlevel_number]]
            randomlevel_length <- dim(samples_omega_randomlevel[[1]])[2]
            if (randomlevel_length > 1000) {
                random_sample_indices <- sample(1:randomlevel_length, size = 1000)
                for (chain in 1:length(samples_omega_randomlevel)) {
                    samples_omega_randomlevel[[chain]] <- samples_omega_randomlevel[[chain]][,random_sample_indices]
                }
            }
            save_convergence_info(samples_omega_randomlevel, 
                                  sprintf("omega, %s", randomlevel_name), 
                                  model_name, 
                                  convergence_file)
            print(sprintf("Calculated convergence for omega %s", randomlevel_name))
        }
        
        
        # CONVERGENCE FOR SPATIAL SCALE ALPHA
        append_to_file("ALPHA\n\n", file = convergence_file)
        samples_alpha <- posterior$Alpha
        for (randomlevel_number in 1:number_of_random_levels) {
            randomlevel_name <- names(model$ranLevels)[randomlevel_number]
            append_to_file(sprintf("%s\n\n", randomlevel_name), file = convergence_file)
            randomlevel <- model$ranLevels[[randomlevel_number]]
            number_of_spatial_dimensions <- randomlevel$sDim
            if (number_of_spatial_dimensions > 0) {
                samples_alpha_randomlevel <- samples_alpha[[randomlevel_number]] 
                save_convergence_info(samples_alpha_randomlevel, 
                                      sprintf("alpha, %s", randomlevel_name), 
                                      model_name, 
                                      convergence_file)
                print(sprintf("Calculated convergence for alpha %s", randomlevel_name))
            } else {
                append_to_file(sprintf("%s is not a spatial random level", randomlevel_name), convergence_file)
            }
        }
        
        
    }
}

# Close pdf
dev.off()
