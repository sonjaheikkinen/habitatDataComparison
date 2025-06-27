# SCRIPT FOR CALCULATING PARAMETER ESTIMATES AND VARIANCE PARTITIONING

support_level <- 0.95

file_parameter_estimates <- file.path(dir_results,"/parameter_estimates.txt")
append_to_file("Additional information regarding parameter estimates\n\n", file_parameter_estimates)


# GET MODELS WITH HIGHEST THINNING VALUE
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

# LOOP TROUGH MODELS
for (model_number in 1:length(fitted_models)) {
    
    # LOAD MODEL INFORMATION
    load(fitted_models[model_number])
    model <- fitted_model
    model_name <- strsplit(basename(fitted_models[model_number]), "\\.")[[1]][1]
    append_to_file(sprintf("\n%s\n\n", model_name), file_parameter_estimates)
    number_of_random_levels <- model$nr
    number_of_species <- model$ns
    number_of_covariates <- model$nc # Includes intercept, and breakdowns for categorical variables
    number_of_traits <- model$nt
    
    # OPEN PDF FOR GRAPHS
    pdf(file = file.path(dir_results, sprintf("parameter_estimates_%s.pdf", model_name)))
    
    # SET NUMBER OF MODELS TO 1 FOR LOOPS (LEFTOVER FROM WHEN THE STRUCTURE WAS DIFFERENT)
    number_of_models <- 1

        
    # EXTRACT NAMES FOR ENVIRONMENTAL VARIABLES
    # Both are needed, because first does not work with categorical variables,
    # second does not work with empty formula
    # Intercept is not included
    if (model$XFormula == "~.") {
        env_variable_names <- colnames(model$X)[-1]
    } else {
        env_variable_names <- attr(terms(model$XFormula),"term.labels")
    }
    number_of_environmental_variables <- length(env_variable_names)
    
    
    # GET VARIANCE PARTITIONING ONLY IF THERE IS ENOUGH VARIABLES
    if ((number_of_random_levels + number_of_environmental_variables) > 1 
        & number_of_species > 1) {
        
        # CALCULATE VARIANCE PARTITIONING
        variance_partitioning <- computeVariancePartitioning(model)
        
        # CALCULATE MODEL FIT
        predicted_values <- computePredictedValues(model)
        model_fit <- evaluateModelFit(hM = model, predY = predicted_values)

        # VARIANCE IN OBSERVATION DATA EXPLAINED BY ENVIRONMENTAL COVARIATES
        partitioning_values <- variance_partitioning$vals
        
        # CALCULATE R2 FOR EACH SPECIES AND BIND TO VARIANCE PARTITIONING
        R2 <- NULL
        if (!is.null(model_fit$TjurR2)) {
            R2 <- model_fit$TjurR2
            partitioning_values <- rbind(partitioning_values, R2)
        }
        if (!is.null(model_fit$R2)) {
            R2 <- model_fit$R2
            partitioning_values <- rbind(partitioning_values, R2)
        }
        if (!is.null(model_fit$SR2)) {
            R2 <- model_fit$SR2
            partitioning_values <- rbind(partitioning_values, R2)
        }
        
        # WRITE VARIANCE PARTITIONING OF COVARIATES PER SPECIES INTO A FILE
        file_variance_partitioning <- file.path(dir_results, 
                                                sprintf("parameter_estimates_VP_%s.csv", 
                                                        model_name))
        write.csv(partitioning_values, file = file_variance_partitioning)
        
        
        # VARIANCE IN SPECIES RESPONSES TO COVARIATES EXPLAINED BY TRAITS
        if (!is.null(variance_partitioning$R2T$Beta)) {
            file_vp_r2t_beta <- file.path(dir_results,
                                          sprintf("parameter_estimates_VP_R2T_Beta_%s.csv", 
                                                  model_name))
            write.csv(variance_partitioning$R2T$Beta, file = file_vp_r2t_beta)
        }
        
        # VARIANCE IN SPECIES OBSERVATION DATA EXPLAINED BY TRAITS
        if (!is.null(variance_partitioning$R2T$Y)) {
            file_vp_r2t_y <- file.path(dir_results, 
                                       sprintf("parameter_estimates_VP_R2T_Y_%s.csv",
                                               model_name))
            write.csv(variance_partitioning$R2T$Y, file = file_vp_r2t_y)
        }
        
        
        # SELECT THE ORDER IN WHICH TO SHOW VARIANCE PARTITIONING
        explained_variance_partition_order <- 1:number_of_species
        #explained_variance_partition_order <- order(R2, decreasing = TRUE)

        # ORDER VALUES IN THE SELECTED ORDER
        variance_partitioning$vals <- variance_partitioning$vals[,explained_variance_partition_order]
        
        # WRITE THE SPECIES NAMES IN CORRECT ORDER TO FILE
        append_to_file("\nexplained variance partitioning order\n\n", 
                       file_parameter_estimates, sep="\n")
        species_names_in_explained_variance_partitioning_order <- model$spNames[explained_variance_partition_order]
        append_to_file(sprintf("%s", species_names_in_explained_variance_partitioning_order), 
                       file_parameter_estimates)
        
        # PLOT VARIANCE PARTITIONING
        partitioning_colors <- rainbow(nrow(variance_partitioning$vals))
        plotVariancePartitioning(hM = model, 
                                 VP = variance_partitioning,
                                 main = sprintf("Proportion of explained variance, %s", model_name), 
                                 cex.main = 0.8, 
                                 cols = partitioning_colors, 
                                 args.leg = list(bg = "white",cex=0.7))
        
        
        # SAME FOR RAW VARIANCE
        raw_variance_partition_order <- 1:number_of_species
        #raw_variance_partition_order <- order(R2, decreasing = TRUE)

        # TRANSFORM VARIANCE PARTITIONING TO RAW VARIANCE PARTITIONING IF R2 EXISTS
        if (!is.null(R2)) {
            
            # INITIALIZE VARIABLE FOR RAW VARIANCES
            raw_variance_partitioning <- variance_partitioning
            
            # GET PROPORTION OF RAW VARIANCE FROM EXPLAINED VARIANCES
            for (k in 1:number_of_species) {
                raw_variance_partitioning$vals[,k] <- R2[k]*raw_variance_partitioning$vals[,k]
            }
            
            # ORDER RAW VARIANCE IN THE SELECTED ORDER
            raw_variance_partitioning$vals <- raw_variance_partitioning$vals[,raw_variance_partition_order]
            
            # APPEND SPECIES NAMES IN CORRECT ORDER TO FILE 
            append_to_file("\nraw variance partitioning order\n\n", file_parameter_estimates)
            append_to_file(sprintf("%s", model$spNames[raw_variance_partition_order]),
                           file_parameter_estimates, 
                           sep="\n")
            
            # PLOT RAW VARIANCE
            plotVariancePartitioning(hM = model, 
                                     VP = raw_variance_partitioning,
                                     main = sprintf("Proportion of raw variance, %s",model_name),
                                     cex.main = 0.8, 
                                     cols = partitioning_colors, 
                                     args.leg = list(bg = "white", cex = 0.7),
                                     ylim = c(0,1))
        }
    }

    
    # CREATE BETA PLOT (EFFECT OF ENVIRONMENTAL VARIABLES ON SPECIES)
    # AND SAVE THE SUPPORT VALUES FOR EACH ESTIMATE
    # IF THERE ARE OTHER COVARIATES BESIDES INTERCEPT
    if (number_of_covariates > 1) {
        
        # GET POSTERIOR ESTIMATES FOR ENVIRONMENTAL VARIABLE EFFECTS
        posterior_estimates_beta <- getPostEstimate(model, parName = "Beta")
        file_beta_estimates <- file.path(dir_results, sprintf("parameter_estimates_Beta_%s.xlsx", model_name))
        mean_estimates <- as.data.frame(t(posterior_estimates_beta$mean))
        mean_estimates <- cbind(model$spNames, mean_estimates)
        colnames(mean_estimates) <- c("Species", model$covNames)
        
        # GET SUPPORT VALUES FOR THE ESTIMATES
        support_estimates <- as.data.frame(t(posterior_estimates_beta$support))
        support_estimates <- cbind(model$spNames, support_estimates)
        colnames(support_estimates) <- c("Species", model$covNames)
        
        # GET NEGATIVE SUPPORT VALUES FOR THE ESTIMATES
        negative_support_estimates <- as.data.frame(t(posterior_estimates_beta$supportNeg))
        negative_support_estimates <- cbind(model$spNames, negative_support_estimates)
        colnames(negative_support_estimates) <- c("Species", model$covNames)
        values <- list("Posterior mean" = mean_estimates,
                       "Pr(x>0)" = support_estimates,
                       "Pr(x<0)" = negative_support_estimates)
        writexl::write_xlsx(values, path = file_beta_estimates)
        
        
        # SHOW SPECIES NAMES IN BETA PLOT IF THERE IS AT MOST 30 SPECIES
        show_species_names <- (is.null(model$phyloTree) && model$ns <= 30) 
        
        # CREATE BETA PLOT
        plot_tree <- !is.null(model$phyloTree)
        plotBeta(model, 
                 post = posterior_estimates_beta, 
                 supportLevel = support_level, 
                 param = "Sign",
                 plotTree = plot_tree,
                 covNamesNumbers = c(TRUE,FALSE),
                 spNamesNumbers = c(show_species_names, FALSE),
                 cex = c(0.6,0.6,0.8))
        title <- sprintf("BetaPlot, %s", model_name)
        if (!is.null(model$phyloTree)) {
            posterior <- convertToCodaObject(model)
            rho_values <- unlist(poolMcmcChains(posterior$Rho))
            title <- sprintf("%s, E[rho] = %f, Pr[rho>0] = %f",
                             title,
                             round(mean(rho_values), 2),
                             round(mean(rho_values > 0), 2))
        }
        title(main = title, line = 2.5, cex.main = 0.8)
    }

    
    # GAMMA PLOT (EFFECTS OF TRAITS ON SPECIES RESPONSES TO COVARIATES)
    # IF TRAITS AND COVARIATES INCLUDE OTHER THINGS BESIDES INTERCEPT
    if (number_of_traits > 1 & number_of_covariates > 1) {
        posterior_gamma <- getPostEstimate(model, parName = "Gamma")
        plotGamma(model, 
                  post = posterior_gamma, 
                  supportLevel = support_level, 
                  param = "Sign",
                  covNamesNumbers = c(TRUE,FALSE),
                  trNamesNumbers = c(number_of_traits < 21, FALSE),
                  cex = c(0.6,0.6,0.8))
        title(main = sprintf("GammaPlot %s", model_name), 
              line = 2.5,
              cex.main = 0.8)
    }
    
    # RANDOM LEVEL ASSOCIATIONS
    # IF THERE ARE RANDOM LEVELS AND MORE THAN 1 SPECIES
    if (number_of_random_levels > 0 & number_of_species > 1) {
        
        # CREATE ASSOCIATION MATRICES
        omega_matrices <- computeAssociations(model)
        
        # GO TROUGH ALL RANDOM LEVELS
        for (random_level_number in 1:number_of_random_levels) {
            
            random_level_support <- omega_matrices[[random_level_number]]$support
            random_level_mean <- omega_matrices[[random_level_number]]$mean
            omega_support_values <- (
                                        (random_level_support > support_level) 
                                        + (random_level_support < (1 - support_level)) > 0
                                     ) * sign(random_level_mean)
            
            show_species_names <- number_of_species <= 30
 
            if (!show_species_names) {
                colnames(omega_support_values) <- rep("", number_of_species)
                rownames(omega_support_values) <- rep("", number_of_species)
            }
            
            # SELECT ORDER IN WHICH TO PLOT THE ASSOCIATIONS
            #plot_order <- 1:number_of_species
            plot_order <- corrMatOrder(random_level_mean, order="AOE")

            # WRITE THE ORDER TO FILE
            append_to_file("\nomega order\n\n", file_parameter_estimates)
            append_to_file(sprintf("%s", model$spNames[plot_order]),
                           file_parameter_estimates,
                           sep = "\n")
            
            
            # CONSTRUCT TITLE FOR PLOT
            title <- sprintf("Associations, %s: %s", 
                             model_name,
                             names(model$ranLevels)[[random_level_number]])
            if (model$ranLevels[[random_level_number]]$sDim > 0) {
                posterior <- convertToCodaObject(model)
                alpha_values <- unlist(poolMcmcChains(posterior$Alpha[[random_level_number]][,1]))
                title <- sprintf("%s, E[alpha%s] = %s, Pr[alpha%s > 0] = %f",
                                 title,
                                 random_level_number,
                                 round(mean(alpha_values), 2),
                                 random_level_number,
                                 round(mean(alpha_values > 0), 2))
            }
            corrplot(omega_support_values[plot_order,plot_order], 
                     method = "color",
                     col = colorRampPalette(c("blue","white","red"))(3),
                     mar = c(0,0,1,0),
                     main = title,
                     cex.main = 0.8)
            
            # WRITE ESTIMATES TO FILE
            mean_estimates <- as.data.frame(omega_matrices[[random_level_number]]$mean)
            mean_estimates <- cbind(model$spNames, mean_estimates)
            colnames(mean_estimates)[1] <- ""
            support_estimates <- as.data.frame(omega_matrices[[random_level_number]]$support)
            support_estimates <- cbind(model$spNames, support_estimates)
            colnames(support_estimates)[1] = ""
            negative_support_estimates <- as.data.frame(1-omega_matrices[[random_level_number]]$support)
            negative_support_estimates <- cbind(model$spNames, negative_support_estimates)
            colnames(negative_support_estimates)[1] <- ""
            values <- list("Posterior mean" = mean_estimates,
                           "Pr(x>0)" = support_estimates,
                           "Pr(x<0)" = negative_support_estimates)
            file_omega = file.path(dir_results, 
                                   sprintf("parameter_estimates_omega_%s_%s.xlsx",
                                           model_name,
                                           names(model$ranLevels[[random_level_number]])))
            writexl::write_xlsx(values,path = file_omega)
        }
    }
    
    dev.off()
    
}



