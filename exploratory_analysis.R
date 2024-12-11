# HELPER FUNCTIONS

create_clusters(data, cluster_amount) {
    
}


# EXPLORATORY ANALYSIS FUNCTIONS

explore_bird_data <- function(occurrence, abundance, spatiotemporal_context, type) {
    
    pdf(file.path(dir_results, sprintf("exploratory_analysis_observations_%s.pdf", type)))
    
    
    
    # Total abundance in samples
    sample_abundances <- rowSums(abundance, na.rm = TRUE)
    hist(sample_abundances,
         main = "Histogram of total abundance in samples",
         xlab = sprintf("Abundance, min %s, max %s", min(sample_abundances), max(sample_abundances)),
         ylab = sprintf("Frequency, %s total samples", nrow(abundance)))
    
    # Total species richness in samples
    sample_occurrences <- rowSums(occurrence, na.rm = TRUE)
    hist(sample_occurrences,
         main = "Histogram of species richness in samples",
         xlab = sprintf("Occurrence, min %s, max %s", min(sample_occurrences), max(sample_occurrences)),
         ylab = sprintf("Frequency, total samples: %s, total richness. %s", 
                        nrow(occurrence),
                        ncol(occurrence)))
    
    
    # Adjust margin to fit species names to the left
    old_par <- par(no.readonly = TRUE) 
    par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
    
    # How many samples was each species found in?
    species_sample_occurrences <- sort(colSums(occurrence))
    barplot(species_sample_occurrences, 
            horiz = TRUE, 
            las = 1, 
            main = "How many samples was each species observed in",
            xlab = "Samples",
            cex.names = 0.7)
    
    # Most common species by occurrence
    common_species_occurrences <- tail(species_sample_occurrences, 10)
    barplot(common_species_occurrences, 
            horiz = TRUE, 
            las = 1, 
            main = "How many samples were common species observed in",
            xlab = "Samples",
            cex.names = 0.7)
    
    # Most rare species by occurrence
    rare_species_occurrences <- species_sample_occurrences[species_sample_occurrences < 5]
    barplot(rare_species_occurrences, 
            horiz = TRUE, 
            las = 1, 
            main = "How many samples were rare species observed in",
            xlab = "Samples",
            cex.names = 0.7)
    
    # Use old margins again
    par(old_par)
    
    
    # Spatial and temporal patterns
    
    # Calculate sample level mean abundances across species
    sample_total_abundances <- matrix(0,
                                      nrow = length(unique(spatiotemporal_context$Year)),
                                      ncol = length(unique(spatiotemporal_context$Transect)))
    rownames(sample_total_abundances) <- rev(sort(unique(spatiotemporal_context$Year)))
    colnames(sample_total_abundances) <- unique(spatiotemporal_context$Transect)
    for (row in 1:nrow(spatiotemporal_context)) {
        year <- spatiotemporal_context[row, ]$Year
        transect <- spatiotemporal_context[row, ]$Transect
        sample <- sprintf("%s%s", year, transect)
        year_index <- match(year, rownames(sample_total_abundances))
        transect_index <- match(transect, colnames(sample_total_abundances))
        row_abundances <- abundance[row, ]
        row_abundances[is.na(row_abundances)] <- 0 
        sample_total_abundances[year_index, transect_index] <- sum(row_abundances)
    }
    pheatmap(sample_total_abundances,
             cluster_rows = FALSE)
    
    
    # Spatial variation in transect total abundance 
    # Transect total abundances: sum of species abundances in each sample
    # Transect mean total abundances: mean of transect total abundances for each transect
    transect_mean_total_abundances <- apply(sample_total_abundances, 
                                            2, 
                                            function(col) {
                                                non_zero_values <- col[col != 0]
                                                if (length(non_zero_values) > 0) {
                                                    mean(non_zero_values)
                                                } else {
                                                    NA  # Return NA if there are no non-zero values
                                                }
                                            })
    spatial_plot_data <- data.frame(Transect = names(transect_mean_total_abundances),
                                    Abundance = transect_mean_total_abundances)
    x_coords <- c()
    y_coords <- c()
    for (transect in spatial_plot_data$Transect) {
        x <- spatiotemporal_context[spatiotemporal_context$Transect == transect,]$x[1]
        x_coords <- c(x_coords, x)
        y <- spatiotemporal_context[spatiotemporal_context$Transect == transect,]$y[1]
        y_coords <- c(y_coords, y)
    }
    spatial_plot_data$x <- x_coords
    spatial_plot_data$y <- y_coords
    plot(spatial_plot_data$x, spatial_plot_data$y, 
         cex = spatial_plot_data$Abundance / max(spatial_plot_data$Abundance) * 1.5, 
         pch = 21, 
         col = "black", # Dot appearance
         xlab = "X Coordinate", ylab = "Y Coordinate", 
         main = "Transect mean total abundances")
    
    
    # Temporal variation in total abundance
    # Scaled by how many samples were taken each year
    year_total_abundances <- apply(sample_total_abundances, 
                                   1, 
                                   function(row) {
                                       non_zero_values <- row[row != 0]
                                       if (length(non_zero_values) > 0) {
                                           mean(non_zero_values)
                                       } else {
                                           NA  # Return NA if there are no non-zero values
                                       }
                                   })
    names(year_total_abundances) <- rownames(sample_total_abundances)
    barplot(year_total_abundances,
            horiz = TRUE,
            las = 1)
    
    
    
    # Spatial variation in species composition
    transect_occurrence_data_list <- list()
    for (transect in unique(spatiotemporal_context$Transect)) {
        transect_indices <- match(transect, spatiotemporal_context$Transect)
        transect_occurrences <- occurrence[transect_indices, , drop = FALSE]
        species_occurrences_in_transect <- apply(transect_occurrences,
                                                 2,
                                                 function(col) {
                                                     any(col == 1) * 1
                                                 })
        transect_occurrence_data_list[transect] <- list(species_occurrences_in_transect)
    }
    transect_occurrence_data <- do.call(rbind, transect_occurrence_data_list)
    
    
    # Plot species clusters over transects
    
    
    
    
    dev.off()
    
}




# SCRIPT STARTS 


# LOAD RAW DATA
load(file = file.path(dir_data, "occurrence_raw.RData")) 
load(file = file.path(dir_data, "abundance_raw.RData")) 
load(file = file.path(dir_data, "spatiotemporal_context_raw.RData"))


explore_bird_data(occurrence, abundance, spatiotemporal_context, "raw")

