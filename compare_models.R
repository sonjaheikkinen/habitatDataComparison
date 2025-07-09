# SCRIPT FOR COMPARING MODELFITS


# FUNCTIONS | GENERAL
extract_thinning_value <- function(filename) {
    parts <- strsplit(filename, "_")[[1]] 
    thinning_value <- as.numeric(parts[length(parts) - 1])
    print(thinning_value)
    return(thinning_value)
}


# FUNCTIONS | MODEL QUALITY


create_modelfit_plot_old <- function(explanatory_power, predictive_power, type, model_name, thinning_value) {
    
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


create_modelfit_plot <- function(explanatory_power, type, model_name, thinning_value, validation_type) {
    
    if (!type %in% colnames(explanatory_power)) {
        return("Type not yet calculated or does not apply for this model")
    }
    
    values <- explanatory_power[,type]
    
    if (!is.null(values)) {
        
        title <- sprintf("%s\n thin = %s | %s | %s \n mean = %s",
                         model_name,
                         as.character(thinning_value),
                         type,
                         validation_type,
                         as.character(mean(explanatory_power[,type], na.rm = TRUE)))
        
        plot(1:length(values),
             values,
             ylim = c(min(values), max(values)),
             xlab = "Species",
             ylab = sprintf("%s", type),
             main = title)
        
        hist(values, xlab = sprintf("%s", type), main = title)
    }
    
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


for (modelfit_file_number in 1:length(modelfit_files)) {
    
    # GET MODEL INFORMATION
    load(modelfit_files[modelfit_file_number])
    model_name <- strsplit(basename(modelfit_files[modelfit_file_number]), "\\.")[[1]][1]
    thinning_value <- extract_thinning_value(model_name)
    modelfit_file <- file.path(dir_modelfits, sprintf("modelfit_%s.RData", model_name))
    
    # CREATE PDF FOR PLOTTING THE FIT VALUES
    pdf(file = file.path(dir_results, sprintf("modelfit_results_%s.pdf", model_name)))
    
    # CREATE PLOTS OF FIT VALUES
    create_modelfit_plot(explanatory_power, "TjurR2", model_name, thinning_value, "Explanatory power")
    create_modelfit_plot(explanatory_power, "R2", model_name, thinning_value, "Explanatory power")
    create_modelfit_plot(explanatory_power, "AUC", model_name, thinning_value, "Explanatory power")
    create_modelfit_plot(explanatory_power, "SR2", model_name, thinning_value, "Explanatory power")
    create_modelfit_plot(explanatory_power, "RMSE", model_name, thinning_value, "Explanatory power")
    create_modelfit_plot(predictive_power_transect, "TjurR2", model_name, thinning_value, "Predictive power")
    create_modelfit_plot(predictive_power_transect, "R2", model_name, thinning_value, "Predictive power")
    create_modelfit_plot(predictive_power_transect, "AUC", model_name, thinning_value, "Predictive power")
    create_modelfit_plot(predictive_power_transect, "SR2", model_name, thinning_value, "Predictive power")
    create_modelfit_plot(predictive_power_transect, "RMSE", model_name, thinning_value, "Predictive power")
    
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


