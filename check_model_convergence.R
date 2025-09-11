# SCRIPT FOR CHECKING MODEL CONVERGENCE

# CONVERGENCE FUNCTIONS



plot_convergence <- function(variable, model_name, psrf_point_estimates) {
    par(mfrow=c(1,2))
    print("")
    print("min")
    print(min(psrf_point_estimates))
    print("max")
    print(max(psrf_point_estimates))
    print("")
    vioplot(psrf_point_estimates, 
            col = "dodgerblue",
            ylim = c(min(psrf_point_estimates) - (max(psrf_point_estimates) - min(psrf_point_estimates)), 
                     max(psrf_point_estimates)))
    vioplot(psrf_point_estimates,
            col = "dodgerblue",
            ylim = c(0.9,1.1))
    title(sprintf("PSRF (%s):\n %s\n[%s, %s], mean %s, median %s, sd %s", 
                  variable, 
                  model_name,
                  round(min(psrf_point_estimates), 3), 
                  round(max(psrf_point_estimates), 3),
                  round(mean(psrf_point_estimates), 3), 
                  round(median(psrf_point_estimates), 3),
                  round(sd(psrf_point_estimates), 3)), 
          outer = TRUE,
          line = -3)
    par(mfrow=c(1,1))
} 

ownEffectiveSize <- function(x)
{
    if (is.mcmc.list(x))
    {
        ##RGA changed to sum across all chains
        ess <- do.call("rbind",lapply(x,ownEffectiveSize))
        ans <- apply(ess,2,sum)
    }
    else
    {
        x <- as.mcmc(x)
        x <- as.matrix(x)
        spec <- test(x)$spec
        ans <- ifelse(spec==0, 0, nrow(x) * apply(x, 2, var)/spec)
    }
    return(ans)
}


test <- function(x)
{
    x <- as.matrix(x)
    v0 <- order <- numeric(ncol(x))
    names(v0) <- names(order) <- colnames(x)
    z <- 1:nrow(x)
    for (i in 1:ncol(x))
    {
        lm.out <- lm(x[,i] ~ z)
        if (identical(all(sd(residuals(lm.out)) ==  0), TRUE)) {
            v0[i] <- 0
            order[i] <- 0
        }
        else {
            ar.out <- ar(x[,i], aic=TRUE)
            v0[i] <- ar.out$var.pred/(1 - sum(ar.out$ar))^2
            order[i] <- ar.out$order
        }
    }
    return(list(spec=v0, order=order))
}

save_convergence_info <- function(samples, variable, model_name, convergence_file) {
    #samples <- filter_zero_columns(samples)
    psrf <- gelman.diag(samples, multivariate = FALSE, autoburnin = FALSE)$psrf
    print(psrf)
    psrf_problem_file <- file.path(dir_results, sprintf("psrf_issues/%s.txt", model_name))
    append_to_file(sprintf("\n\n%s\n", variable), psrf_problem_file)
    for (row_number in 1:nrow(psrf)) {
        row_name <- rownames(psrf)[row_number]
        row_values <- psrf[row_number,]
        psrf_point_estimate <- row_values[1]
        if (psrf_point_estimate > 1.1) {
            append_to_file(sprintf("%s %s\n", row_name, psrf_point_estimate), psrf_problem_file)
        }
    }
    append_to_file(sprintf("%s\n", summary(psrf)[,1]), file = convergence_file)
    append_to_file("\n\n", file = convergence_file)
    psrf_point_estimates <- psrf[,1]
    plot_convergence(variable, model_name, psrf_point_estimates)
    effective_sizes <- ownEffectiveSize(samples)
    effectiveSize_problem_file <- file.path(dir_results, sprintf("effectiveSize_issues/%s.txt", model_name))
    append_to_file(sprintf("\n\n%s\n", variable), effectiveSize_problem_file)
    for (index in 1:length(effective_sizes)) {
        name <- names(effective_sizes)[index]
        effective_size <- effective_sizes[index]
        if (effective_size < 1000) {
            append_to_file(sprintf("%s %s\n", name, effective_size), effectiveSize_problem_file)
        }
    }
    plot_histogram(effective_sizes,
                   sprintf("Effective size of %s", variable))
    combined_samples <- as.mcmc.list(samples)
    traceplot(combined_samples[,1],
              main = sprintf("Traceplot of %s, %s", 
                             variable,
                             colnames(combined_samples[[1]])[1]))
    # TO DO: Visualize traces for all variables, not just first
    par(mfrow = c(2, 2))
    acf_colors <- rainbow(ncol(as.data.frame(combined_samples[[1]])))
    acf_problem_file <- file.path(dir_results, sprintf("acf_issues/%s.txt", model_name))
    append_to_file(sprintf("\n\n%s\n", variable), acf_problem_file)
    for (chain in 1:4) {
        append_to_file(sprintf("\n\nChain %s\n", chain), acf_problem_file)
        combined_samples_chain <- as.data.frame(combined_samples[[chain]])
        for (column in 1:ncol(combined_samples_chain)) {
            acf <- acf(as.numeric(combined_samples_chain[,column]),
                       plot = FALSE)
            confidence_upper_limit <- qnorm((1 + (1 - 0.05))/2)/sqrt(acf$n.used)
            confidence_lower_limit <- -confidence_upper_limit
            if (column == 1) {
                plot(acf$lag, acf$acf, 
                     col = acf_colors[column],
                     ylim = c(-1, 1),
                     main = sprintf("Chain %s", chain),
                     xlab = "Lag", ylab = "ACF")
                abline(confidence_upper_limit, 0, col = "lightgrey")
                abline(confidence_lower_limit, 0, col = "lightgrey")
            } else {
                points(acf$lag, 
                       acf$acf,
                       col = acf_colors[column])
            }
            if (any(abs(acf$acf[2:length(acf$acf)]) > 0.2)) {
                column_name <- colnames(combined_samples_chain)[column]
                append_to_file(sprintf("\n%s", column_name), acf_problem_file)
            }
        }
    }
    title(main = sprintf("Autocorrelation of %s", variable),
          outer = TRUE,
          line = -3)
    
}

# Remove those columns from samples, that have been set to zero across all chains
filter_zero_columns <- function(samples) {
    sample_colsums <- lapply(samples, colSums)
    columns_zero_in_any_chain <- Reduce('|', lapply(sample_colsums, function(value) value == 0))
    filtered_samples <- lapply(samples, function(chain) chain[, !columns_zero_in_any_chain, drop = FALSE])
    return(filtered_samples)  
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
    samples_beta <- posterior$Beta
    save_convergence_info(samples_beta, "Beta - species niches", model_name, convergence_file)
    print("Calculated convergence for beta")
    
    # CONVERGENCE FOR INFLUENCE OF TRAITS (GAMMA)
    append_to_file("GAMMA\n\n", file = convergence_file)
    samples_gamma <- posterior$Gamma
    save_convergence_info(samples_gamma, "Gamma - influence of traits", model_name, convergence_file)
    print("Calculated convergence for gamma")
    
    # CONVERGENCE FOR UNEXPLAINED VARIATION IN SPECIES NICHES (V)
    append_to_file("V\n\n", file = convergence_file)
    samples_v <- posterior$V
    save_convergence_info(samples_v, "V - unexplained variation in species niches", model_name, convergence_file)
    print("Calculated convergence for V")
    
    # CONVERGENCE FOR RESIDUAL VARIANCE (SIGMA)
    #append_to_file("SIGMA\n\n", file = convergence_file)
    #samples_sigma <- posterior$Sigma
    #save_convergence_info(samples_sigma, "Sigma - residual variance", model_name, convergence_file)
    #print("Calculated convergence for Sigma")
    
    
    # CONVERGENCE FOR INFLUENCE OF PHYLOGENY (RHO)
    if (!is.null(posterior$Rho)) {
        append_to_file("RHO\n\n", file = convergence_file)
        samples_rho <- posterior$Rho
        save_convergence_info(samples_rho, "Rho - influence of phylogeny", model_name, convergence_file)
        print("Calculated convergence for rho")
    }
    
    if (number_of_random_levels > 0) {
        
        
        for (randomlevel_number in 1:number_of_random_levels) {
            
            randomlevel_name <- names(model$ranLevels)[[randomlevel_number]]
            
       
            # CONVERGENCE FOR SITE LOADINGS (ETA)
            append_to_file("ETA\n\n", file = convergence_file)
            samples_eta <- posterior$Eta
            append_to_file(sprintf("%s\n\n", randomlevel_name), file = convergence_file)
            samples_eta_randomlevel <- samples_eta[[randomlevel_number]]
            randomlevel_length <- dim(samples_eta_randomlevel[[1]])[2]
            if (randomlevel_length > 1000) {
                random_sample_indices <- sample(1:randomlevel_length, size = 1000)
                for (chain in 1:length(samples_eta_randomlevel)) {
                    samples_eta_randomlevel[[chain]] <- samples_eta_randomlevel[[chain]][,random_sample_indices]
                }
            }
            save_convergence_info(samples_eta_randomlevel, 
                                  sprintf("Eta - site loadings, %s", randomlevel_name), 
                                  model_name, 
                                  convergence_file)
            print(sprintf("Calculated convergence for Eta %s", randomlevel_name))
        
        
        
            # CONVERGENCE FOR SPECIES LOADINGS (LAMBDA)
            append_to_file("LAMBDA\n\n", file = convergence_file)
            samples_lambda <- posterior$Lambda
            append_to_file(sprintf("%s\n\n", randomlevel_name), file = convergence_file)
            samples_lambda_randomlevel <- samples_lambda[[randomlevel_number]]
            randomlevel_length <- dim(samples_lambda_randomlevel[[1]])[2]
            if (randomlevel_length > 1000) {
                random_sample_indices <- sample(1:randomlevel_length, size = 1000)
                for (chain in 1:length(samples_lambda_randomlevel)) {
                    samples_lambda_randomlevel[[chain]] <- samples_lambda_randomlevel[[chain]][,random_sample_indices]
                }
            }
            save_convergence_info(samples_lambda_randomlevel, 
                                  sprintf("Lambda - species loadings, %s", randomlevel_name), 
                                  model_name, 
                                  convergence_file)
            print(sprintf("Calculated convergence for Lambda %s", randomlevel_name))
        
        
        
            # CONVERGENCE FOR RANDOM VARIATION IN CO-OCCURRENCE (OMEGA)
            append_to_file("OMEGA\n\n", file = convergence_file)
            samples_omega <- posterior$Omega
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
                                  sprintf("Omega - random variation in co-occurrence, %s", randomlevel_name), 
                                  model_name, 
                                  convergence_file)
            print(sprintf("Calculated convergence for omega %s", randomlevel_name))
        
        
        
            # CONVERGENCE FOR SPATIAL SCALE ALPHA
            append_to_file("ALPHA\n\n", file = convergence_file)
            samples_alpha <- posterior$Alpha
            append_to_file(sprintf("%s\n\n", randomlevel_name), file = convergence_file)
            randomlevel <- model$ranLevels[[randomlevel_number]]
            number_of_spatial_dimensions <- randomlevel$sDim
            if (number_of_spatial_dimensions > 0) {
                samples_alpha_randomlevel <- samples_alpha[[randomlevel_number]] 
                save_convergence_info(samples_alpha_randomlevel, 
                                      sprintf("Alpha - spatial scale, %s", randomlevel_name), 
                                      model_name, 
                                      convergence_file)
                print(sprintf("Calculated convergence for alpha %s", randomlevel_name))
            } else {
                append_to_file(sprintf("%s is not a spatial random level", randomlevel_name), convergence_file)
            }
        
        
        
            # CONVERGENCE FOR LOCAL SHRINKAGE OF SPECIES LOADINGS (PSI)
            append_to_file("PSI\n\n", file = convergence_file)
            samples_psi <- posterior$Psi
            append_to_file(sprintf("%s\n\n", randomlevel_name), file = convergence_file)
            samples_psi_randomlevel <- samples_psi[[randomlevel_number]]
            randomlevel_length <- dim(samples_psi_randomlevel[[1]])[2]
            if (randomlevel_length > 1000) {
                random_sample_indices <- sample(1:randomlevel_length, size = 1000)
                for (chain in 1:length(samples_psi_randomlevel)) {
                    samples_psi_randomlevel[[chain]] <- samples_psi_randomlevel[[chain]][,random_sample_indices]
                }
            }
            save_convergence_info(samples_psi_randomlevel, 
                                  sprintf("Psi - local shrinkage of species loadings (?), %s", randomlevel_name), 
                                  model_name, 
                                  convergence_file)
            print(sprintf("Calculated convergence for Psi %s", randomlevel_name))
        
        
        
            # CONVERGENCE FOR GLOBAL SHRINKAGE OF SPECIES LOADINGS (DELTA)
            append_to_file("DELTA\n\n", file = convergence_file)
            samples_delta <- posterior$Delta
            append_to_file(sprintf("%s\n\n", randomlevel_name), file = convergence_file)
            samples_delta_randomlevel <- samples_delta[[randomlevel_number]]
            randomlevel_length <- dim(samples_delta_randomlevel[[1]])[2]
            if (randomlevel_length > 1000) {
                random_sample_indices <- sample(1:randomlevel_length, size = 1000)
                for (chain in 1:length(samples_delta_randomlevel)) {
                    samples_delta_randomlevel[[chain]] <- samples_delta_randomlevel[[chain]][,random_sample_indices]
                }
            }
            save_convergence_info(samples_delta_randomlevel, 
                                  sprintf("Delta - global shrinkage of species loadings, %s", randomlevel_name), 
                                  model_name, 
                                  convergence_file)
            print(sprintf("Calculated convergence for Delta %s", randomlevel_name))
            
            
        }
        
        
        
    }
}

# Close pdf
dev.off()
