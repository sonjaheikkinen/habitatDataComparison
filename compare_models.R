# SCRIPT FOR COMPARING MODELFITS

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

# EXTRACT MODELFIT VALUES FROM EACH FILE
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


