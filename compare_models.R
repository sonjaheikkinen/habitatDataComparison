# SCRIPT FOR COMPARING MODELFITS


# FUNCTIONS | GENERAL
extract_thinning_value <- function(filename) {
    parts <- strsplit(filename, "_")[[1]] 
    thinning_value <- as.numeric(parts[length(parts) - 1])
    print(thinning_value)
    return(thinning_value)
}


# FUNCTIONS | MODEL QUALITY


create_modelfit_plot <- function(explanatory_power, predictive_power, type, model_name, thinning_value) {
    
    if (!type %in% colnames(explanatory_power)) {
        return("Type not yet calculated or does not apply for this model")
    }
    if (!is.null(explanatory_power[,type])) {
        plot(explanatory_power[,type],
             predictive_power[,type],
             xlim = c(-1,1),
             ylim = c(-1,1),
             xlab = "explanatory power (MF)",
             ylab = "predictive power (MFCV)",
             main = sprintf("%s\n thin = %s: %s. \nmean(MF) = %s, mean(MFCV) = %s",
                            model_name,
                            as.character(thinning_value),
                            type,
                            as.character(mean(explanatory_power[,type], na.rm = TRUE)),
                            as.character(mean(predictive_power[,type], na.rm = TRUE))))
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

plot_fits_against_continuous <- function(xdata, xlab, predictive_power, explanatory_power) {
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
         predictive_power$AUC,
         xlab = xlab,
         ylab = "AUC",
         main = "Predictive power")
    plot(xdata, 
         predictive_power$RMSE,
         xlab = xlab,
         ylab = "RMSE",
         main = "Predictive power")
    plot(xdata, 
         predictive_power$TjurR2,
         xlab = xlab,
         ylab = "TjurR2",
         main = "Predictive power")
}

plot_fits_against_categorical <- function(xdata, xlab, explanatory_power, predictive_power) {
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
    boxplot(predictive_power$AUC ~ xdata,
            xlab = xlab,
            ylab = "AUC",
            main = "Predictive power")
    boxplot(predictive_power$RMSE ~ xdata,
            xlab = xlab,
            ylab = "RMSE",
            main = "Predictive power")
    boxplot(predictive_power$TjurR2 ~ xdata,
            xlab = xlab,
            ylab = "TjurR2",
            main = "Predictive power")
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

# LOAD IN MODELFITS
modelfit_files <- list.files(dir_modelfits, pattern="*.RData", full.names=TRUE)
load(file.path(dir_data, "species_prevalences.RData"))
load(file.path(dir_data, "trait_data.RData"))
trait_data <- as.data.frame(trait_data)
rownames(trait_data) <- trait_data$Species



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
    if (!is.null(predictive_power_transect)) {
        create_modelfit_plot(explanatory_power, 
                             predictive_power_transect, 
                             "TjurR2", 
                             model_name, 
                             thinning_value)
        create_modelfit_plot(explanatory_power, 
                             predictive_power_transect,
                             "AUC", 
                             model_name, 
                             thinning_value)
        create_modelfit_plot(explanatory_power,
                             predictive_power_transect,
                             "RMSE", 
                             model_name, 
                             thinning_value)
    }
    
    # LOOK AT VALUES FOR EACH SPECIES
    explanatory_power$Species <- rownames(explanatory_power)
    predictive_power$Species <- rownames(predictive_power)
    old_par <- par(no.readonly = TRUE) 
    par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
    par(mfrow = c(1, 3)) 
    plot_fits_per_species("TjurR2", explanatory_power)
    plot_fits_per_species("AUC", explanatory_power)
    plot_fits_per_species("RMSE", explanatory_power)
    title(main = "Explanatory power",
          outer = TRUE,
          line = -3)
    plot_fits_per_species("TjurR2", predictive_power)
    plot_fits_per_species("AUC", predictive_power)
    plot_fits_per_species("RMSE", predictive_power)
    title(main = "Predictive power",
          outer = TRUE,
          line = -3)
    par(old_par)
    
    
    # COMPARE EACH VALUE FOR EACH SPECIES
    
    explanatory_power_ordered <- explanatory_power[order(explanatory_power$AUC), ]
    explanatory_power_ordered$Species <- NULL
    long_explanatory_power <- stack(explanatory_power_ordered)
    long_explanatory_power$Species <- rep(rownames(explanatory_power_ordered), 3)
    long_explanatory_power$Species <- factor(long_explanatory_power$Species, 
                                             levels = rev(rownames(explanatory_power_ordered)))
    ggplot(long_explanatory_power, aes(x = values, y = Species, fill = ind)) +
        geom_col(position = "identity", alpha = 0.5) +
        labs(x = "Value", y = "Species", title = "Explanatory power") +
        scale_fill_manual(values = c("red", "yellow", "blue")) +
        theme_minimal()
    
    predictive_power_ordered <- predictive_power[order(predictive_power$AUC), ]
    predictive_power_ordered$Species <- NULL
    long_predictive_power <- stack(predictive_power_ordered)
    long_predictive_power$Species <- rep(rownames(predictive_power_ordered), 3)
    long_predictive_power$Species <- factor(long_predictive_power$Species, 
                                             levels = rev(rownames(predictive_power_ordered)))
    ggplot(long_predictive_power, aes(x = values, y = Species, fill = ind)) +
        geom_col(position = "identity", alpha = 0.5) +
        labs(x = "Value", y = "Species", title = "Predictive power") +
        scale_fill_manual(values = c("red", "yellow", "blue")) +
        theme_minimal()
    
    
    # FIT AND PREDICTIVE POWER AGAINST PREVALENCE 
    plot_fits_against_continuous(species_prevalences[rownames(explanatory_power),]$prevalences,
                                 "Prevalence",
                                 explanatory_power,
                                 predictive_power)
    
    
    # PLOTS AGAINST TRAITS
    plot_fits_against_categorical(trait_data[rownames(explanatory_power), ]$Feeding,
                                  "Feeding",
                                  explanatory_power,
                                  predictive_power)
    plot_fits_against_categorical(trait_data[rownames(explanatory_power), ]$Mig,
                                  "Migration",
                                  explanatory_power,
                                  predictive_power)
    plot_fits_against_continuous(trait_data[rownames(explanatory_power), ]$Mass,
                                 "Mass",
                                 explanatory_power,
                                 predictive_power)
    plot_fits_against_continuous(trait_data[rownames(explanatory_power), ]$LogMass,
                                 "LogMass",
                                 explanatory_power,
                                 predictive_power)
    
    

    
    # CLOSE PDF
    dev.off()
    
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


