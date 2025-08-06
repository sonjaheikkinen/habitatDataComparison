# SCRIPT FOR COMPARING MODELFITS


# FUNCTIONS | GENERAL
extract_thinning_value <- function(filename) {
    parts <- strsplit(filename, "_")[[1]] 
    thinning_value <- as.numeric(parts[length(parts) - 1])
    print(thinning_value)
    return(thinning_value)
}


# FUNCTIONS | MODEL QUALITY


create_modelfit_plot <- function(explanatory_power,
                                 predictive_power, 
                                 predictive_power_type, 
                                 type, 
                                 model_name, 
                                 thinning_value) {
    
    if (!type %in% colnames(explanatory_power)) {
        return("Type not yet calculated or does not apply for this model")
    }
    if (!is.null(explanatory_power[,type])) {
        plot(explanatory_power[,type],
             predictive_power[,type],
             xlim = c(-1,1),
             ylim = c(-1,1),
             xlab = "explanatory power",
             ylab = sprintf("predictive power %s", predictive_power_type),
             cex.main = 0.8,
             main = sprintf("%s\n thin = %s: %s. \nexp pw: [%.3f, %.3f], mean = %.3f, sd = %.3f \npred pw: [%.3f, %.3f], mean = %.3f, sd = %.3f\n",
                            model_name,
                            as.character(thinning_value),
                            type,
                            min(explanatory_power[,type]),
                            max(explanatory_power[,type]),
                            mean(explanatory_power[,type], na.rm = TRUE),
                            sd(explanatory_power[,type]),
                            min(predictive_power[,type]),
                            max(predictive_power[,type]),
                            mean(predictive_power[,type], na.rm = TRUE),
                            sd(predictive_power[,type])))
        abline(h = 0)
        abline(v = 0)
    }
    
}


plot_fits_per_species <- function(type, data) {
    values <- data[,type]
    names(values) <- data$Species
    ordered_values <- sort(values, decreasing = FALSE)
    barplot(ordered_values,
            horiz = TRUE,
            las = 1,                    
            main = type)
}

plot_fits_against_continuous <- function(xdata, 
                                         xlab, 
                                         predictive_power, 
                                         explanatory_power,
                                         predictive_power_year,
                                         waic) {
    plot(xdata, 
         explanatory_power$AUC,
         xlab = xlab,
         ylab = "AUC",
         main = "Explanatory power")
    plot(xdata, 
         explanatory_power$RMSE,
         xlab = xlab,
         ylab = "RMSE",
         main = "Explanatory power")
    plot(xdata, 
         explanatory_power$TjurR2,
         xlab = xlab,
         ylab = "TjurR2",
         main = "Explanatory power")
    
    plot(xdata, 
         predictive_power_transect$AUC,
         xlab = xlab,
         ylab = "AUC",
         main = "Predictive power transect")
    plot(xdata, 
         predictive_power_transect$RMSE,
         xlab = xlab,
         ylab = "RMSE",
         main = "Predictive power transect")
    plot(xdata, 
         predictive_power_transect$TjurR2,
         xlab = xlab,
         ylab = "TjurR2",
         main = "Predictive power transect")
    
    plot(xdata, 
         predictive_power_year$AUC,
         xlab = xlab,
         ylab = "AUC",
         main = "Predictive power year")
    plot(xdata, 
         predictive_power_year$RMSE,
         xlab = xlab,
         ylab = "RMSE",
         main = "Predictive power year")
    plot(xdata, 
         predictive_power_year$TjurR2,
         xlab = xlab,
         ylab = "TjurR2",
         main = "Predictive power year")
    
    plot(xdata, 
         waic,
         xlab = xlab,
         ylab = "WAIC",
         main = "WAIC")
}

plot_fits_against_categorical <- function(xdata, 
                                          xlab, 
                                          explanatory_power, 
                                          predictive_power_transect, 
                                          predictive_power_year,
                                          waic) {
    boxplot(explanatory_power$AUC ~ xdata,
            xlab = xlab,
            ylab = "AUC",
            main = "Explanatory power")
    boxplot(explanatory_power$RMSE ~ xdata,
            xlab = xlab,
            ylab = "RMSE",
            main = "Explanatory power")
    boxplot(explanatory_power$TjurR2 ~ xdata,
            xlab = xlab,
            ylab = "TjurR2",
            main = "Explanatory power")
    
    boxplot(predictive_power_transect$AUC ~ xdata,
            xlab = xlab,
            ylab = "AUC",
            main = "Predictive power transect")
    boxplot(predictive_power_transect$RMSE ~ xdata,
            xlab = xlab,
            ylab = "RMSE",
            main = "Predictive power transect")
    boxplot(predictive_power_transect$TjurR2 ~ xdata,
            xlab = xlab,
            ylab = "TjurR2",
            main = "Predictive power transect")
    
    boxplot(predictive_power_year$AUC ~ xdata,
            xlab = xlab,
            ylab = "AUC",
            main = "Predictive power year")
    boxplot(predictive_power_year$RMSE ~ xdata,
            xlab = xlab,
            ylab = "RMSE",
            main = "Predictive power year")
    boxplot(predictive_power_year$TjurR2 ~ xdata,
            xlab = xlab,
            ylab = "TjurR2",
            main = "Predictive power year")
    
    boxplot(waic ~ xdata,
            xlab = xlab,
            ylab = "WAIC",
            main = "WAIC")
}





# FUNCTIONS | MODEL COMPARISON

# FUNCTION FOR PLOTTING MODEL PERFORMANCE METRICS
plot_metric <- function(data) {
    metric_name <- data$Metric[1]
    ggplot(data, 
           aes(x = Value, y = reorder(Model, Value))) +
        geom_bar(stat = "identity", fill = "dodgerblue") +
        geom_text(aes(label = round(Value, 3)), hjust = -0.2, size = 4) +
        labs(x = "Value", y = "Model", title = paste("Model Comparison -", metric_name)) +
        theme_minimal()
}

# SCRIPT STARTS

# LOAD IN MODELFITS AND DATA FOR COMPARISONS
modelfit_files <- list.files(dir_modelfits, pattern="*.RData", full.names=TRUE)
load(file.path(dir_data, "species_prevalences.RData"))
load(file.path(dir_data, "trait_data.RData"))
trait_data <- as.data.frame(trait_data)
rownames(trait_data) <- trait_data$Species
load(file = file.path(dir_data, "phylogeny_data.RData"))



# QUALITY OF MODELS


for (modelfit_file_number in 1:length(modelfit_files)) {
    
    # GET MODEL INFORMATION
    load(modelfit_files[modelfit_file_number])
    model_name <- strsplit(basename(modelfit_files[modelfit_file_number]), "\\.")[[1]][1]
    thinning_value <- extract_thinning_value(model_name)
    modelfit_file <- file.path(dir_modelfits, sprintf("modelfit_%s.RData", model_name))
    
    # CREATE PDF FOR PLOTTING THE FIT VALUES
    pdf(file = file.path(dir_results, sprintf("modelfit_results_%s.pdf", model_name)))
    
    # CREATE PLOTS FOR PREDICTIVE POWER VS. EXPLANATORY POWER
    create_modelfit_plot(explanatory_power, 
                         predictive_power_transect, 
                         "Transect",
                         "TjurR2", 
                         model_name, 
                         thinning_value)
    create_modelfit_plot(explanatory_power, 
                         predictive_power_transect,
                         "Transect",
                         "AUC", 
                         model_name, 
                         thinning_value)
    create_modelfit_plot(explanatory_power,
                         predictive_power_transect,
                         "Transect",
                         "RMSE", 
                         model_name, 
                         thinning_value)
    create_modelfit_plot(explanatory_power, 
                         predictive_power_year, 
                         "Year",
                         "TjurR2", 
                         model_name, 
                         thinning_value)
    create_modelfit_plot(explanatory_power, 
                         predictive_power_year,
                         "Year",
                         "AUC", 
                         model_name, 
                         thinning_value)
    create_modelfit_plot(explanatory_power,
                         predictive_power_year,
                         "Year",
                         "RMSE", 
                         model_name, 
                         thinning_value)
    
    old_par <- par(no.readonly = TRUE) 
    par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
    waic_values <- waic_by_column
    names(waic_values) <- rownames(explanatory_power)
    ordered_values <- sort(waic_values, decreasing = FALSE)
    barplot(ordered_values,
            horiz = TRUE,
            las = 1,    
            cex.names = 0.7,
            main = sprintf("%s\n thin = %s: %s. \n[%.3f, %.3f], mean = %.3f, sd = %.3f",
                           model_name,
                           as.character(thinning_value),
                           "WAIC",
                           min(waic_values),
                           max(waic_values),
                           mean(waic_values, na.rm = TRUE),
                           sd(waic_values)))
    par(old_par)
    
    
    

    
    # LOOK AT VALUES FOR EACH SPECIES
    explanatory_power$Species <- rownames(explanatory_power)
    predictive_power_transect$Species <- rownames(predictive_power_transect)
    predictive_power_year$Species <- rownames(predictive_power_year)
    waic_df <- data.frame(WAIC = waic_by_column)
    waic_df$Species <- rownames(explanatory_power)
    old_par <- par(no.readonly = TRUE) 
    par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
    par(mfrow = c(1, 3)) 
    plot_fits_per_species("TjurR2", explanatory_power)
    plot_fits_per_species("AUC", explanatory_power)
    plot_fits_per_species("RMSE", explanatory_power)
    title(main = "Explanatory power",
          outer = TRUE,
          line = -3)
    plot_fits_per_species("TjurR2", predictive_power_transect)
    plot_fits_per_species("AUC", predictive_power_transect)
    plot_fits_per_species("RMSE", predictive_power_transect)
    title(main = "Predictive power transect",
          outer = TRUE,
          line = -3)
    plot_fits_per_species("TjurR2", predictive_power_year)
    plot_fits_per_species("AUC", predictive_power_year)
    plot_fits_per_species("RMSE", predictive_power_year)
    title(main = "Predictive power year",
          outer = TRUE,
          line = -3)
    plot_fits_per_species("WAIC", waic_df)
    par(old_par)
    
    
    # SPECIES COMPARISON
    # For metrics where smaller is better = multiply by -1 to change the direction
    # Scale all metrics to 0-1 for comparison
    species_comparison_df <- data.frame(exp_pw_TjurR2 = scale_to_range(explanatory_power$TjurR2, 0, 1),
                                        exp_pw_AUC = scale_to_range(explanatory_power$AUC, 0, 1),
                                        exp_pw_RMSE = scale_to_range(-1 * explanatory_power$RMSE, 0, 1),
                                        pred_pw_transect_TjurR2 = scale_to_range(predictive_power_transect$TjurR2, 0, 1),
                                        pred_pw_transect_AUC = scale_to_range(predictive_power_transect$AUC, 0, 1),
                                        pred_pw_transect_RMSE = scale_to_range(-1 * predictive_power_transect$RMSE, 0, 1),
                                        pred_pw_year_TjurR2 = scale_to_range(predictive_power_year$TjurR2, 0, 1),
                                        pred_pw_year_AUC = scale_to_range(predictive_power_year$AUC, 0, 1),
                                        pred_pw_year_RMSE = scale_to_range(-1 * predictive_power_year$RMSE, 0, 1),
                                        waic = scale_to_range(-1 * waic_by_column, 0, 1))
    rownames(species_comparison_df) <- rownames(explanatory_power)
    pheatmap(species_comparison_df, 
             main = sprintf("Species model fit comparison \n%s", model_name))
    
    
    
    
    # FIT AND PREDICTIVE POWER AGAINST PREVALENCE 
    plot_fits_against_continuous(species_prevalences[rownames(explanatory_power),]$prevalences,
                                 "Prevalence",
                                 explanatory_power,
                                 predictive_power_transect,
                                 predictive_power_year,
                                 waic_by_column)
    
    
    # PLOTS AGAINST TRAITS
    plot_fits_against_categorical(trait_data[rownames(explanatory_power), ]$Feeding,
                                  "Feeding",
                                  explanatory_power,
                                  predictive_power_transect,
                                  predictive_power_year,
                                  waic_by_column)
    plot_fits_against_categorical(trait_data[rownames(explanatory_power), ]$Mig,
                                  "Migration",
                                  explanatory_power,
                                  predictive_power_transect,
                                  predictive_power_year,
                                  waic_by_column)
    plot_fits_against_continuous(trait_data[rownames(explanatory_power), ]$Mass,
                                 "Mass",
                                 explanatory_power,
                                 predictive_power_transect,
                                 predictive_power_year,
                                 waic_by_column)
    plot_fits_against_continuous(trait_data[rownames(explanatory_power), ]$LogMass,
                                 "LogMass",
                                 explanatory_power,
                                 predictive_power_transect,
                                 predictive_power_year,
                                 waic_by_column)
    

    
    # CLOSE PDF
    dev.off()
    
}

# COMPARE MODEL PERFORMANCE BETWEEN THINNING VALUES
fitted_models <- list.files(dir_fitted, pattern="*.RData", full.names=TRUE)
load(fitted_models[5])
natura_thin_100 <- fitted_model
load(fitted_models[6])
natura_thin_1000 <- fitted_model
explanatory_power_natura_thin_100 <- evaluateModelFit(hM = natura_thin_100, 
                                                      predY = computePredictedValues(natura_thin_100))
explanatory_power_natura_thin_100 <- as.data.frame(explanatory_power_natura_thin_100)
explanatory_power_natura_thin_1000 <- evaluateModelFit(hM = natura_thin_1000, 
                                                      predY = computePredictedValues(natura_thin_1000))
explanatory_power_natura_thin_1000 <- as.data.frame(explanatory_power_natura_thin_1000)
cor(explanatory_power_natura_thin_100$RMSE, explanatory_power_natura_thin_1000$RMSE)
cor(explanatory_power_natura_thin_100$AUC, explanatory_power_natura_thin_1000$AUC)
cor(explanatory_power_natura_thin_100$TjurR2, explanatory_power_natura_thin_1000$TjurR2)

posterior_natura_thin_100 <- convertToCodaObject(natura_thin_100, 
                                                 spNamesNumbers = c(T,F), 
                                                 covNamesNumbers = c(T,F))
posterior_natura_thin_1000 <- convertToCodaObject(natura_thin_1000, 
                                                 spNamesNumbers = c(T,F), 
                                                 covNamesNumbers = c(T,F))

fixed_effect_parameter_names  <- c("Beta", "Gamma", "V")
random_effect_parameter_names <- c("Eta", "Lambda", "Omega", "Alpha")
randomlevel_names <- names(natura_thin_100$ranLevels) 

for (chain in 1:4) {
    print(sprintf("Chain %s", chain))
    print("")
    for (parameter_name in fixed_effect_parameter_names) {
        natura_thin_100_parameter_values <- as.vector(posterior_natura_thin_100 [[parameter_name]][[chain]])
        natura_thin_1000_parameter_values <- as.vector(posterior_natura_thin_1000[[parameter_name]][[chain]])
        print(sprintf("Natura %s correlation between 100 and 1000, chain %s: %s",
                parameter_name,
                chain,
                cor(natura_thin_100_parameter_values, natura_thin_1000_parameter_values)))
        print("")
    }
    for (parameter_name in random_effect_parameter_names) {
        for (randomlevel_number in 1:2) {
            natura_thin_100_parameter_values <- as.vector(posterior_natura_thin_100 [[parameter_name]][[randomlevel_number]][[chain]])
            natura_thin_1000_parameter_values <- as.vector(posterior_natura_thin_1000[[parameter_name]][[randomlevel_number]][[chain]])
            print(sprintf("Natura %s %s correlation between 100 and 1000, chain %s: %s",
                    parameter_name,
                    randomlevel_names[randomlevel_number],
                    chain,
                    cor(natura_thin_100_parameter_values, natura_thin_1000_parameter_values)))
        }
        print("")
    }
    print("")
}


load(fitted_models[2])
corine_thin_100 <- fitted_model
load(fitted_models[3])
corine_thin_1000 <- fitted_model
explanatory_power_corine_thin_100 <- evaluateModelFit(hM = corine_thin_100, 
                                                      predY = computePredictedValues(corine_thin_100))
explanatory_power_corine_thin_100 <- as.data.frame(explanatory_power_corine_thin_100)
explanatory_power_corine_thin_1000 <- evaluateModelFit(hM = corine_thin_1000, 
                                                       predY = computePredictedValues(corine_thin_1000))
explanatory_power_corine_thin_1000 <- as.data.frame(explanatory_power_corine_thin_1000)
cor(explanatory_power_corine_thin_100$RMSE, explanatory_power_corine_thin_1000$RMSE)
cor(explanatory_power_corine_thin_100$AUC, explanatory_power_corine_thin_1000$AUC)
cor(explanatory_power_corine_thin_100$TjurR2, explanatory_power_corine_thin_1000$TjurR2)

posterior_corine_thin_100 <- convertToCodaObject(corine_thin_100, 
                                                 spNamesNumbers = c(T,F), 
                                                 covNamesNumbers = c(T,F))
posterior_corine_thin_1000 <- convertToCodaObject(corine_thin_1000, 
                                                  spNamesNumbers = c(T,F), 
                                                  covNamesNumbers = c(T,F))

fixed_effect_parameter_names  <- c("Beta", "Gamma", "V")
random_effect_parameter_names <- c("Eta", "Lambda", "Omega", "Alpha")
randomlevel_names <- names(natura_thin_100$ranLevels) 


for (chain in 1:4) {
    print(sprintf("Chain %s", chain))
    print("")
    for (parameter_name in fixed_effect_parameter_names) {
        corine_thin_100_parameter_values <- as.vector(posterior_corine_thin_100 [[parameter_name]][[chain]])
        corine_thin_1000_parameter_values <- as.vector(posterior_corine_thin_1000[[parameter_name]][[chain]])
        print(sprintf("Corine %s correlation between 100 and 1000, chain %s: %s",
                      parameter_name,
                      chain,
                      cor(corine_thin_100_parameter_values, corine_thin_1000_parameter_values)))
        print("")
    }
    for (parameter_name in random_effect_parameter_names) {
        for (randomlevel_number in 1:2) {
            corine_thin_100_parameter_values <- as.vector(posterior_corine_thin_100 [[parameter_name]][[randomlevel_number]][[chain]])
            corine_thin_1000_parameter_values <- as.vector(posterior_corine_thin_1000[[parameter_name]][[randomlevel_number]][[chain]])
            print(sprintf("Corine %s %s correlation between 100 and 1000, chain %s: %s",
                          parameter_name,
                          randomlevel_names[randomlevel_number],
                          chain,
                          cor(corine_thin_100_parameter_values, corine_thin_1000_parameter_values)))
        }
        print("")
    }
    print("")
}





# START WITH OVERALL COMPARISONS

# EXTRACT MEAN MODELFIT VALUES FROM EACH FILE
waic_list <- list()
rmse_list <- list()
auc_list <- list()
tjurr2_list <- list()
for (file_number in 1:length(modelfit_files)) {
    
    model_name <- strsplit(basename(modelfit_files[file_number]), "\\.")[[1]][1]
    
    load(modelfit_files[file_number])
    
    waic_list[[model_name]] <- waic
    rmse_list[[model_name]] <- mean(explanatory_power$RMSE)
    auc_list[[model_name]] <- mean(explanatory_power$AUC)
    tjurr2_list[[model_name]] <- mean(explanatory_power$TjurR2)
    
    
}


# CONVERT LISTS TO DATAFRAMES FOR PLOTTING
waic_df <- data.frame(Model = names(waic_list), Value = unlist(waic_list), Metric = "WAIC")
rmse_df <- data.frame(Model = names(rmse_list), Value = unlist(rmse_list), Metric = "RMSE")
auc_df <- data.frame(Model = names(auc_list), Value = unlist(auc_list), Metric = "AUC")
tjurr2_df <- data.frame(Model = names(tjurr2_list), Value = unlist(tjurr2_list), Metric = "TjurR2")


# PLOT METRICS
plot_metric(waic_df)
plot_metric(rmse_df)
plot_metric(auc_df)
plot_metric(tjurr2_df)


rmse_df <- NULL
auc_df<- NULL
tjurr2_df <- NULL

# NEXT COMPARISONS BY SPECIES
for (file_number in 1:length(modelfit_files)) {
    
    model_name <- strsplit(basename(modelfit_files[file_number]), "\\.")[[1]][1]
    
    load(modelfit_files[file_number])
    
    initial_dataframe <- data.frame(species = rownames(explanatory_power))
    
    if (is.null(rmse_df)) {
        rmse_df <- initial_dataframe
    }
    if (is.null(auc_df)) {
        auc_df <- initial_dataframe
    }
    if(is.null(tjurr2_df)) {
        tjurr2_df <- initial_dataframe
    }
    
    rmse_df[,model_name] <- explanatory_power$RMSE
    auc_df[,model_name] <- explanatory_power$AUC
    tjurr2_df[,model_name] <- explanatory_power$TjurR2
    
}

rownames(rmse_df) <- rmse_df$species
rmse_df$species <- NULL
rownames(auc_df) <- auc_df$species
auc_df$species <- NULL
rownames(tjurr2_df) <- tjurr2_df$species
tjurr2_df$species <- NULL

pheatmap(rmse_df)
pheatmap(auc_df)
pheatmap(tjurr2_df)

# FIT LINEAR MODEL TO EXPLAIN FIT USING MODELLING APPROACH
# TO DO: Add species as random effect?
long_rmse <- stack(rmse_df)
long_auc <- stack(auc_df)
long_tjurr2 <- stack(tjurr2_df)

lm_rmse <- lm(values ~ ind, data = long_rmse)
summary(lm_rmse)
lm_auc <- lm(values ~ ind, data = long_auc)
summary(lm_auc)
lm_tjurr2 <- lm(values ~ ind, data = long_tjurr2)
summary(lm_tjurr2)


plot(auc_df$modelfit_probit_corine_thin_300_fitted, 
     auc_df$modelfit_probit_natura_thin_300_fitted)
abline(a = 0, b = 1)


