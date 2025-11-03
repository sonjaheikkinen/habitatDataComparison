# SCRIPT FOR CALCULATING PARAMETER ESTIMATES AND VARIANCE PARTITIONING

support_level <- 0.95

file_parameter_estimates <- file.path(dir_results,"/parameter_estimates.txt")
append_to_file("Additional information regarding parameter estimates\n\n", file_parameter_estimates)


# GET MODELS WITH HIGHEST THINNING VALUE
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)

variance_partitionings <- list()
raw_variance_partitionings <- list()

associations <- list()


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
        variance_partitionings[[model_name]] <- variance_partitioning
        
        
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
        #partitioning_colors <- rainbow(nrow(variance_partitioning$vals))
        partitioning_colors <- c("gold", "yellow", "dodgerblue", "blue", "lightblue", "darkgreen", "green", "lightgreen", "forestgreen")
        partitioning_colors <- rev(partitioning_colors)
        plotVariancePartitioning(hM = model, 
                                 VP = variance_partitioning,
                                 main = sprintf("Proportion of explained variance, %s", model_name), 
                                 cex.main = 0.8, 
                                 cols = partitioning_colors, 
                                 args.leg = list(bg = "white",cex=0.7))
        
        
        # SAME FOR RAW VARIANCE
        #raw_variance_partition_order <- 1:number_of_species
        raw_variance_partition_order <- order(R2, decreasing = TRUE)

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
            raw_variance_partitionings[[model_name]] <- raw_variance_partitioning
            
            
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
                                     las = 2,
                                     cex.names = 0.5,
                                     cols = partitioning_colors, 
                                     args.leg = list(bg = "white", cex = 0.7),
                                     ylim = c(0,1))
            
            
            ranked_vals <- apply(variance_partitioning$vals, 2, rank, ties.method = "first")
            
            # Keep the same row and column names
            rownames(ranked_vals) <- rownames(variance_partitioning$vals)
            colnames(ranked_vals) <- colnames(variance_partitioning$vals)
            
        }
    }
    
    # CREATE BETA PLOT (EFFECT OF ENVIRONMENTAL VARIABLES ON SPECIES)
    # AND SAVE THE SUPPORT VALUES FOR EACH ESTIMATE
    # IF THERE ARE OTHER COVARIATES BESIDES INTERCEPT
    if (number_of_covariates > 1) {
        
        # GET POSTERIOR ESTIMATES FOR ENVIRONMENTAL VARIABLE EFFECTS
        posterior_estimates_beta <- getPostEstimate(model, parName = "Beta")
        file_beta_estimates <- file.path(dir_results, sprintf("parameter_estimates_Beta_%s.xlsx", model_name))
        # posterior mean df: covariates x species
        beta_posterior_means <- as.data.frame(t(posterior_estimates_beta$mean))
        beta_posterior_means <- cbind(model$spNames, beta_posterior_means)
        colnames(beta_posterior_means) <- c("Species", model$covNames)
        
        # GET SUPPORT VALUES FOR THE ESTIMATES
        support_for_beta_positive <- as.data.frame(t(posterior_estimates_beta$support))
        support_for_beta_positive <- cbind(model$spNames, support_for_beta_positive)
        colnames(support_for_beta_positive) <- c("Species", model$covNames)
        
        # GET NEGATIVE SUPPORT VALUES FOR THE ESTIMATES
        support_for_beta_negative <- as.data.frame(t(posterior_estimates_beta$supportNeg))
        support_for_beta_negative <- cbind(model$spNames, support_for_beta_negative)
        colnames(support_for_beta_negative) <- c("Species", model$covNames)
        
        values <- list("Posterior mean" = beta_posterior_means,
                       "Pr(x>0)" = support_for_beta_positive,
                       "Pr(x<0)" = support_for_beta_negative)
        writexl::write_xlsx(values, path = file_beta_estimates)
        
        
        # CREATE BETA PLOT
        
        old_par <- par(no.readonly = TRUE) 
        par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
        plotBeta(model, 
                 post = posterior_estimates_beta, 
                 supportLevel = support_level, 
                 param = "Sign",
                 covNamesNumbers = c(TRUE,FALSE),
                 spNamesNumbers = c(TRUE, FALSE),
                 colors = colorRampPalette(c("blue", "white", "red")),
                 cex = c(0.6,0.6,0.8))
        title(main = sprintf("BetaPlot, sign of posterior estimate \n%s", 
                             model_name), 
              line = -0.5)
        plotBeta(model, 
                 post = posterior_estimates_beta, 
                 supportLevel = support_level, 
                 param = "Mean",
                 covNamesNumbers = c(TRUE,FALSE),
                 spNamesNumbers = c(TRUE, FALSE),
                 colors = colorRampPalette(c("blue", "white", "red")),
                 cex = c(0.6,0.6,0.8))
        title(main = sprintf("BetaPlot, mean of posterior estimate \n%s", 
                             model_name), 
              line = -0.5)
        beta_posterior_means_for_heatmap <- beta_posterior_means[,3:ncol(beta_posterior_means)]
        support_for_beta_positive_for_heatmap <- support_for_beta_positive[,3:ncol(support_for_beta_positive)]
        support_for_beta_negative_for_heatmap <- support_for_beta_negative[,3:ncol(support_for_beta_negative)]
        cells_with_no_support <- (support_for_beta_positive_for_heatmap < support_level) & (support_for_beta_negative_for_heatmap < support_level)
        beta_posterior_means_for_heatmap[cells_with_no_support] <- NA
                                                                                
        beta_minimum_posterior_mean <- min(beta_posterior_means_for_heatmap, na.rm = TRUE)
        beta_maximum_posterior_mean <- max(beta_posterior_means_for_heatmap, na.rm = TRUE)
        number_of_negative_breaks <- 100
        number_of_positive_breaks <- 100
        negative_breaks <- seq(beta_minimum_posterior_mean, 0, length.out = number_of_negative_breaks)
        positive_breaks <- seq(0, beta_maximum_posterior_mean, length.out = number_of_positive_breaks)
        negative_breaks <- negative_breaks[1:length(negative_breaks) - 1]
        
        color_palette_breaks <- c(negative_breaks, positive_breaks)
        colors <- c(
            colorRampPalette(c("blue", "lightblue"))(number_of_negative_breaks - 1),  # zero was removed  
            colorRampPalette(c("pink", "red"))(number_of_positive_breaks)  
        )
        pheatmap(beta_posterior_means_for_heatmap,
                 color = colors,
                 breaks = color_palette_breaks,
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 main = "Posterior means without intercept")
        pheatmap(t(apply(beta_posterior_means[,3:ncol(beta_posterior_means)], 1, rank)),
                 main = "Ranked posterior means without intercept")
        par(old_par)
    }
    
    # PLOT CREDIBLE INTERVAL
    
    

    
    # GAMMA PLOT (EFFECTS OF TRAITS ON SPECIES RESPONSES TO COVARIATES)
    # IF TRAITS AND COVARIATES INCLUDE OTHER THINGS BESIDES INTERCEPT
    if (number_of_traits > 1 & number_of_covariates > 1) {
        
        # GET POSTERIOR ESTIMATES FOR GAMMA
        posterior_estimates_gamma <- getPostEstimate(model, parName = "Gamma")
        gamma_posterior_means <- as.data.frame(t(posterior_estimates_gamma$mean))
        gamma_posterior_means <- cbind(model$trNames, gamma_posterior_means)
        rownames(gamma_posterior_means) <- model$trNames
        colnames(gamma_posterior_means) <- c("Trait", model$covNames)
        
        # GET SUPPORT VALUES FOR THE ESTIMATES
        support_for_gamma_positive <- as.data.frame(t(posterior_estimates_gamma$support))
        support_for_gamma_positive <- cbind(model$trNames, support_for_gamma_positive)
        rownames(support_for_gamma_positive) <- model$trNames
        colnames(support_for_gamma_positive) <- c("Trait", model$covNames)
        
        # GET NEGATIVE SUPPORT VALUES FOR THE ESTIMATES
        support_for_gamma_negative <- as.data.frame(t(posterior_estimates_gamma$supportNeg))
        support_for_gamma_negative <- cbind(model$trNames, support_for_gamma_negative)
        rownames(support_for_gamma_negative) <- model$trNames
        colnames(support_for_gamma_negative) <- c("Trait", model$covNames)
        
        
        plotGamma(model, 
                  post = posterior_estimates_gamma, 
                  supportLevel = support_level, 
                  param = "Sign",
                  covNamesNumbers = c(TRUE,FALSE),
                  trNamesNumbers = c(TRUE, FALSE),
                  cex = c(0.6,0.6,0.8))
        title(main = sprintf("GammaPlot sign %s", model_name), 
              line = 2.5,
              cex.main = 0.8)

        
        gamma_posterior_means_for_heatmap <- gamma_posterior_means[,3:ncol(gamma_posterior_means)]
        support_for_gamma_positive_for_heatmap <- support_for_gamma_positive[,3:ncol(support_for_gamma_positive)]
        support_for_gamma_negative_for_heatmap <- support_for_gamma_negative[,3:ncol(support_for_gamma_negative)]
        cells_with_no_support <- (support_for_gamma_positive_for_heatmap < support_level) & (support_for_gamma_negative_for_heatmap < support_level)
        gamma_posterior_means_for_heatmap[cells_with_no_support] <- NA
        
        max_absolute_value_of_gamma <- max(abs(max(gamma_posterior_means_for_heatmap, na.rm = TRUE)),
                                          abs(min(gamma_posterior_means_for_heatmap, na.rm = TRUE)))
        color_palette_breaks <- seq(-max_absolute_value_of_gamma, 
                                    max_absolute_value_of_gamma,
                                    length.out = 101)
        pheatmap(gamma_posterior_means_for_heatmap,
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 breaks = color_palette_breaks,
                 cluster_cols = FALSE,
                 cluster_rows = FALSE,
                 main = "Gamma posterior means without intercept")
        pheatmap(t(apply(gamma_posterior_means[,3:ncol(gamma_posterior_means)], 1, rank)),
                 main = "Ranked gamma posterior means without intercept")
        
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
            
            show_species_names <- TRUE
 
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
            title <- sprintf("Association signs, %s: %s", 
                             model_name,
                             names(model$ranLevels)[[random_level_number]])
            if (model$ranLevels[[random_level_number]]$sDim > 0) {
                posterior <- convertToCodaObject(model)
                alpha_values <- unlist(poolMcmcChains(posterior$Alpha[[random_level_number]][,1]))
                title <- sprintf("%s\n, E[alpha%s] = %s, Pr[alpha%s > 0] = %f",
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
            
            # HEATMAP FOR CORRELATIONS
            species_associations_df <- random_level_mean
            species_associations_df <- omega_support_values * species_associations_df
            pheatmap(species_associations_df[plot_order, plot_order],
                     cluster_rows = FALSE,
                     cluster_cols = FALSE,
                     main = sprintf("Species associations %s\n%s",
                                    names(model$ranLevels)[[random_level_number]],
                                    model_name))
            
            # SAVE ASSOCIATIONS FOR CORRELATION TEST
            association_name <- sprintf("%s_%s", 
                                        names(model$ranLevels)[random_level_number],
                                        model_name)
            associations[[association_name]] <- species_associations_df[plot_order, plot_order]
            
            
            # WRITE ESTIMATES TO FILE
            omega_posterior_means <- as.data.frame(omega_matrices[[random_level_number]]$mean)
            omega_posterior_means <- cbind(model$spNames, omega_posterior_means)
            colnames(omega_posterior_means)[1] <- ""
            support_for_omega_positive <- as.data.frame(omega_matrices[[random_level_number]]$support)
            support_for_omega_positive <- cbind(model$spNames, support_for_omega_positive)
            colnames(support_for_omega_positive)[1] = ""
            support_for_omega_negative <- as.data.frame(1-omega_matrices[[random_level_number]]$support)
            support_for_omega_negative <- cbind(model$spNames, support_for_omega_negative)
            colnames(support_for_omega_negative)[1] <- ""
            values <- list("Posterior mean" = omega_posterior_means,
                           "Pr(x>0)" = support_for_omega_positive,
                           "Pr(x<0)" = support_for_omega_negative)
            file_omega = file.path(dir_results, 
                                   sprintf("parameter_estimates_omega_%s_%s.xlsx",
                                           model_name,
                                           names(model$ranLevels[[random_level_number]])))
            writexl::write_xlsx(values,path = file_omega)
        }
    }
    
    dev.off()
    
}


# VARIANCE PLOTS AND CORRELATIONS


plot(colSums(raw_variance_partitionings[[1]]$vals),
     colSums(raw_variance_partitionings[[2]]$vals))

cor(colSums(raw_variance_partitionings[[1]]$vals),
     colSums(raw_variance_partitionings[[2]]$vals))



corine_habitat_variables <- c("Effort:Havumetsät.kivennäismaalla",
                              "Effort:Sekametsät.kivennäismaalla",
                              "Effort:Sekametsät.turvemaalla",
                              "Effort:Lehtimetsät.kivennäismaalla",
                              "Effort:Havumetsät.kalliomaalla")

comparable_variances_corine <- data.frame(t(raw_variance_partitionings[[1]]$vals), check.names = FALSE)
comparable_variances_corine$HabitatVariables <- rowSums(comparable_variances_corine[,corine_habitat_variables])
comparable_variances_corine <- t(comparable_variances_corine[,setdiff(colnames(comparable_variances_corine),
                                                              corine_habitat_variables)])

natura_habitat_variables <- c("Effort:Luonnonmetsät",
                              "Effort:Tunturikoivikot",
                              "Effort:Lehdot",
                              "Effort:Tulvametsät")

comparable_variances_natura <- data.frame(t(raw_variance_partitionings[[2]]$vals), check.names = FALSE)
comparable_variances_natura$HabitatVariables <- rowSums(comparable_variances_natura[,natura_habitat_variables])
comparable_variances_natura <- t(comparable_variances_natura[,setdiff(colnames(comparable_variances_natura),
                                                                      natura_habitat_variables)])


colors <- rainbow(ncol(comparable_variances_corine))

plot(comparable_variances_corine[,1] ,
     comparable_variances_natura[,1],
     col = colors[1])
for (i in 2:ncol(comparable_variances_corine)) {
    points(comparable_variances_corine[,i],
           comparable_variances_natura[,i],
           col = colors[i])
}

correlations <- c()
for (species in colnames(comparable_values_corine)) {
    correlations <- c(correlations,
                      cor(comparable_variances_corine[,species],
                          comparable_variances_natura[,species]))
}

correlation_data <- data.frame(correlation = correlations,
                               species = colnames(comparable_values_corine))

ggplot(correlation_data, aes(x = reorder(species, correlation), y = correlation)) +
    geom_bar(stat = "identity", fill = "dodgerblue") +
    labs(title = "Explained variance correlation between models", x = "Species", y = "Correlation") +
    coord_flip()
    theme_minimal()



rowmeans_long_df <- data.frame(mean = c(rowMeans(comparable_variances_corine),
                                        rowMeans(comparable_variances_natura)),
                               variable = c(rownames(comparable_variances_corine),
                                            rownames(comparable_variances_natura)),
                               type = c(rep("corine", nrow(comparable_variances_corine)),
                                        rep("natura", nrow(comparable_variances_natura))))



ggplot(rowmeans_long_df, aes(x = variable, y = mean, fill = type)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = "Mean explained variance by each variable", x = "Variable", y = "Mean explained variance") +
    theme_minimal()



# SPECIES ASSOCIATIONS PLOTS AND CORRELATIONS

corine_transect_associations <- associations$Transect_3_probit_corine_forest_thin_100_samples_500_fitted
natura_transect_associations <- associations$Transect_3_probit_natura_forest_thin_100_samples_500_fitted

natura_transect_associations <- natura_transect_associations[rownames(corine_transect_associations), 
                                                             colnames(corine_transect_associations)]


# Extract upper triangle values
corine_upper_triangle_indices <- which(upper.tri(corine_transect_associations, diag = F), arr.ind = TRUE )
corine_values <- data.frame(col = dimnames(corine_transect_associations)[[2]][corine_upper_triangle_indices[,2]],
                            row = dimnames(corine_transect_associations)[[1]][corine_upper_triangle_indices[,1]],
                            value = corine_transect_associations[corine_upper_triangle_indices])
rownames(corine_values) <- paste(corine_values$row, corine_values$col, sep = " ")
natura_upper_triangle_indices <- which(upper.tri(natura_transect_associations, diag = F), arr.ind = TRUE )
natura_values <- data.frame(col = dimnames(natura_transect_associations)[[2]][natura_upper_triangle_indices[,2]],
                            row = dimnames(natura_transect_associations)[[1]][natura_upper_triangle_indices[,1]],
                            value = natura_transect_associations[natura_upper_triangle_indices])
rownames(natura_values) <- paste(natura_values$row, natura_values$col, sep = " ")






