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


# SCRIPT STARTS

# List filenames of fitted models
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

# Create file for convergence plots
pdf(file = file.path(dir_results,"convergence_results.pdf"))

# Create text file for convergence results
convergence_file <- file.path(dir_results, "convergence_results.txt")

for (i in 1:length(fitted_models)) {
    
    # GET MODEL
    load(fitted_models[i]) # Load file into fitted_model
    model <- fitted_model
    model_name <- strsplit(basename(fitted_models[i]), "\\.")[[1]][1]
    append_to_file(model_name, convergence_file, sep = "\n\n")
    thinning_value <- strsplit(model_name, "_")[[1]][4]
    
    # GET THE POSTERIOR SAMPLES OF PARAMETERS AS CODA OBJECT
    # Get only names (not numbers) of species and covariates
    posterior <- convertToCodaObject(model, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
    
    # CALCULATE THE CONVERGENCE FOR EACH VARIABLE
    # First calculate the potential scale reduction factors
    # Then write the point estimate summaries (min, max, mean, median, quantiles) to file
    # Finally, create column in the corresponding dataframe to save the actual estimates for plotting
    
    
    # CONVERGENCE OF SPECIES NICHES (BETA)
    psrf_beta <- gelman.diag(posterior$Beta, multivariate=FALSE)$psrf
    append_to_file("beta\n\n", file = convergence_file)
    append_to_file(sprintf("%s\n", summary(psrf_beta)[,1]), file = convergence_file)
    append_to_file("\n\n", file = convergence_file)
    beta_psrf_point_estimates <- psrf_beta[,1]
    plot_convergence("beta", model_name, beta_psrf_point_estimates)
    
    # CONVERGENCE OF INFLUENCE OF TRAITS (GAMMA)
    psrf_gamma <- gelman.diag(posterior$Gamma, multivariate=FALSE)$psrf
    append_to_file("gamma\n\n", file = convergence_file)
    append_to_file(sprintf("%s\n", summary(psrf_gamma)[,1]), file = convergence_file)
    append_to_file("\n\n", file = convergence_file)
    gamma_psrf_point_estimates <- psrf_gamma[,1]
    plot_convergence("gamma", model_name, gamma_psrf_point_estimates)    

}

# Close pdf
dev.off()
