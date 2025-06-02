# LOAD RAW DATA
observations <- read.csv(file.path(dir_data, "observations_not_year_limited.csv"),
                         sep = ";")

# LOAD SELECTED DATA
load(file = file.path(dir_data, "occurrence.RData")) 
load(file = file.path(dir_data, "abundance.RData")) 
load(file = file.path(dir_data, "spatiotemporal_context.RData"))
load(file = file.path(dir_data, "transect_coordinates.RData"))
load(file = file.path(dir_data, "env_data_natura.RData")) 
load(file = file.path(dir_data, "env_data_corine.RData"))
load(file = file.path(dir_data, "fractions_natura.RData")) 
load(file = file.path(dir_data, "fractions_corine.RData"))
load(file = file.path(dir_data, "pca_results_natura.RData"))
load(file = file.path(dir_data, "pca_results_corine.RData"))
load(file = file.path(dir_data, "phylogeny_data.RData"))
load(file = file.path(dir_data, "trait_data.RData"))
transects_shp <- vect(file.path(dir_data, "transects_selected.shp"))
natura_classification <- read_excel(file.path(dir_data, "Ylalappi_luokitus.xls"))
corine_classification <- read_excel(file.path(dir_data, "CorineMaanpeite2018Luokat.xls"))


# NAME MAPPING FOR HABITAT
names_natura <- data.frame(value = natura_classification$Value, name = natura_classification$NaturaTyyppi)
names_corine <- data.frame(value = corine_classification$Value, name = corine_classification$Level4Suo)

# TRANSECT NAMES FOR COORDINATE DATA
transect_coordinates$Transect <- rownames(transect_coordinates)




pdf(file.path(dir_results, sprintf("exploratory_analysis_final.pdf")))


# SPATIOTEMPORAL_CONTEXT OVERVIEW
years <- sort(unique(observations$Year))
transects <- unique(observations$SiteID)
samples <- matrix(0, nrow = length(transects), ncol = length(years))
colnames(samples) <- years
rownames(samples) <- transects
for (sample in 1:nrow(observations)) {
    year <- as.character(observations[sample, ]$Year)
    transect <- observations[sample, ]$SiteID
    samples[transect, year] <- 1
}
pheatmap(samples,
         color = c("lightblue", "navy"),
         cluster_cols = FALSE,
         legend_breaks = c(0, 1),
         main = "Transect visits")



# HABITAT RASTERS OVERVIEW
par(mfrow = c(1, 2))
plot(natura_raster,
     main = "Natura",
     legend = FALSE)
lines(transects_shp, col = "red")
corine_cropped <- crop(corine_raster, natura_raster)
plot(corine_cropped, 
     main = "Corine",
     legend = FALSE)
lines(transects_shp, col = "red")
title(main = "Natura and Corine types in study area, with included transects",
      outer = TRUE,
      line = -3)




# HABITAT TYPES OVERVIEW

# Frequency of values in the entire raster
natura_total_frequencies <- freq(natura_raster)
natura_total_frequencies$fraction <- natura_total_frequencies$count / sum(natura_total_frequencies$count)
natura_total_frequencies$region <- "Total"
corine_cropped <- crop(corine_raster, natura_raster)
corine_total_frequencies <- freq(corine_cropped)
corine_total_frequencies$fraction <- corine_total_frequencies$count / sum(corine_total_frequencies$count)
corine_total_frequencies$region <- "Total"

# Read lines and buffer
buffer_area <- buffer(transects_shp, width = buffer_width)

# Mask raster with buffer
natura_masked <- mask(natura_raster, buffer_area)
corine_masked <- mask(corine_cropped, buffer_area)

# Frequency within buffer
natura_buffer_frequencies <- freq(natura_masked)
natura_buffer_frequencies$fraction <- natura_buffer_frequencies$count / sum(natura_buffer_frequencies$count)
natura_buffer_frequencies$region <- "Buffer"
corine_buffer_frequencies <- freq(corine_masked)
corine_buffer_frequencies$fraction <- corine_buffer_frequencies$count / sum(corine_buffer_frequencies$count)
corine_buffer_frequencies$region <- "Buffer"

# Combine both data frames
natura_frequencies_combined <- rbind(natura_total_frequencies, natura_buffer_frequencies)
natura_frequencies_combined$value <- as.factor(natura_frequencies_combined$value)
corine_frequencies_combined <- rbind(corine_total_frequencies, corine_buffer_frequencies)
corine_frequencies_combined$value <- as.factor(corine_frequencies_combined$value)

# Add habitat type names
corine_frequencies_combined <- merge(
    corine_frequencies_combined, 
    names_corine[, c("value", "name")], 
    by.x = "value", 
    by.y = "value",
    all.x = TRUE
)

# Add habitat type names
natura_frequencies_combined <- merge(
    natura_frequencies_combined, 
    names_natura[, c("value", "name")], 
    by.x = "value", 
    by.y = "value",
    all.x = TRUE
)


natura_barplot <- ggplot(natura_frequencies_combined, aes(x = name, y = fraction, fill = region)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() + 
    labs(
        title = "1A Natura total fractions",
        x = "Habitat type",
        y = "Fraction",
        fill = "Region"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Total" = "dodgerblue", "Buffer" = "darkblue"))

corine_barplot <- ggplot(corine_frequencies_combined, aes(x = name, y = fraction, fill = region)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(
        title = "1B Corine total fractions",
        x = "Land cover type",
        y = "Fraction",
        fill = "Region"
    ) +
    theme_minimal() +
    scale_fill_manual(values = c("Total" = "dodgerblue", "Buffer" = "darkblue"))





# Convert transect habitat fractions to long format
natura_fractions_long <- data.frame(
    variable = rep(colnames(fractions_natura), each = nrow(fractions_natura)),
    fraction = as.vector(as.matrix(fractions_natura))
)
corine_fractions_long <- data.frame(
    variable = rep(colnames(fractions_corine), each = nrow(fractions_corine)),
    fraction = as.vector(as.matrix(fractions_corine))
)


natura_boxplot <- ggplot(natura_fractions_long, aes(x = variable, y = fraction)) +
    geom_boxplot() +
    coord_flip() +
    labs(title = "2A Natura transect fractions", 
         x = "Habitat type", 
         y = "Fraction") +
    theme_minimal()

corine_boxplot <- ggplot(corine_fractions_long, aes(x = variable, y = fraction)) +
    geom_boxplot() +
    coord_flip() +
    labs(title = "2B Corine transect fractions", 
         x = "Land cover type", 
         y = "Fraction") +
    theme_minimal()


# Arrange all plots to grid
grid.arrange(
    natura_barplot,
    corine_barplot,
    natura_boxplot,
    corine_boxplot,
    ncol = 2
)



par(mfrow = c(3, 4))

for (habitat_type in colnames(fractions_natura)) {
    plot_histogram(env_data_natura[,habitat_type],
                   title = sprintf("%s", habitat_type))
}


par(mfrow = c(5, 4))

for (habitat_type in colnames(fractions_corine)) {
    plot_histogram(env_data_corine[,habitat_type],
                   title = sprintf("%s", habitat_type))
}

par(mfrow = c(3, 2))
plot_histogram(env_data_natura$Temperature,
               title = "Temperature (April and May average, celcius)")
plot_histogram(env_data_natura$Rainfall,
               title = "Rainfall (April and May average, mm)") 
plot_histogram(env_data_natura$PatchDensity,
               title = "Natura patch density within buffer areas")
plot_histogram(env_data_corine$PatchDensity,
               title = "Corine patch density within buffer areas")
plot_histogram(env_data_natura$Effort,
               title = "Transect length (meters)")




 

# SPATIAL AND TEMPORAL PLOTS OF ENVIRONMENTAL DATA


# Average over transect
natura_transect_averages <- aggregate(env_data_natura, 
                             by = list(Transect = spatiotemporal_context$Transect),
                             FUN = mean, na.rm = TRUE)
natura_transect_averages <- merge(natura_transect_averages, 
                                  transect_coordinates, 
                                  by.x = "Transect", 
                                  by.y = "Transect",
                                  all.x = TRUE)
corine_transect_averages <- aggregate(env_data_corine, 
                                      by = list(Transect = spatiotemporal_context$Transect),
                                      FUN = mean, na.rm = TRUE)
corine_transect_averages <- merge(corine_transect_averages, 
                                  transect_coordinates, 
                                  by.x = "Transect", 
                                  by.y = "Transect",
                                  all.x = TRUE)

# Average over year
natura_year_averages <- aggregate(env_data_natura, 
                                  by = list(Year = spatiotemporal_context$Year),
                                  FUN = mean, na.rm = TRUE)
corine_year_averages <- aggregate(env_data_corine, 
                                  by = list(Year = spatiotemporal_context$Year),
                                  FUN = mean, na.rm = TRUE)



plot_spatially <- function(data, variable, title, text = FALSE) {
    color_palette <- colorRampPalette(c("blue", "white", "red"))
    number_of_colors <- 100 
    scaled_values <- scale_to_range(data[, variable], 1, number_of_colors)
    colors <- color_palette(number_of_colors)[round(scaled_values)]
    plot(data[, c("x")], data[, c("y")], 
         col = colors, pch = 19,
         xlab = "X Coordinate", 
         ylab = "Y Coordinate",
         main = title)
    if (text) {
        text(x = data[, c("x")],
             y = data[, c("y")],
             labels = round(data[, variable], 1),  
             pos = 3, 
             col = "black", 
             cex = 0.5)
    }
}

par(mfrow = c(2, 4))
plot_spatially(natura_transect_averages, "Cluster", title = "Natura cluster", text=TRUE)
plot_spatially(corine_transect_averages, "Cluster", title = "Corine cluster", text=TRUE)
plot_spatially(natura_transect_averages, "Temperature", title = "Average temperature", text=TRUE)
plot_spatially(natura_transect_averages, "Rainfall", title = "Average rainfall", text=TRUE)
plot_spatially(natura_transect_averages, "PatchDensity", title = "Natura Patch Density", text=TRUE)
plot_spatially(corine_transect_averages, "PatchDensity", title = "Corine Patch Density", text=TRUE)
plot_spatially(natura_transect_averages, "Effort", title = "Transect length", text=TRUE)


plot_temporally <- function(data, variable, title) {
    barplot(data[, variable], 
            names.arg = data[, c("Year")], 
            xlab = "Year", 
            ylab = sprintf("%s", variable),
            main = title)
}
par(mfrow = c(2, 1))
plot_temporally(natura_year_averages, "Temperature", title = "Average temperature")
plot_temporally(natura_year_averages, "Rainfall", title = "Average rainfall")




# WITHIN DATASET CORRELATIONS
natura_correlation_columns <- c(colnames(fractions_natura),
                                "Temperature",
                                "Rainfall",
                                "PatchDensity")
corine_correlation_columns <- c(colnames(fractions_corine),
                                "Temperature",
                                "Rainfall",
                                "PatchDensity")
pheatmap(cor(env_data_natura[,natura_correlation_columns]), 
         display_numbers = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
pheatmap(cor(env_data_corine[,corine_correlation_columns]), 
         display_numbers = TRUE,
         cluster_rows = FALSE,
         cluster_cols = FALSE)



# ABUNDANCE DISTRIBUTIONS
abundance_df <- as.data.frame(abundance)
abundance_long <- stack(abundance_df)
colnames(abundance_long) <- c("Abundance", "Species")

ggplot(abundance_long, aes(x = Abundance, color = Species)) +
    geom_freqpoly(binwidth = 5, size = 1) +
    scale_colour_manual(values = rainbow(ncol(abundance))) +
    theme_minimal() +
    labs(title = "Species abundance histograms as lines",
         x = "Abundance", y = "Count") +
    theme(legend.position = "none")

#ggplot(abundance_long, aes(x = Abundance, fill = Species)) +
#    geom_histogram(position = "identity", alpha = 0.2, bins = 30, color = "black") +
#    scale_colour_manual(values = rainbow(ncol(abundance))) +
#    theme_minimal() +
#    labs(title = "Overlapping Histograms of Species Abundance",
#         x = "Abundance", y = "Count") +
#    theme(legend.position = "none")

#ggplot(abundance_long, aes(x = Abundance, color = Species)) +
#    stat_bin(geom = "step", binwidth = 5, position = "identity", size = 1) +
#    scale_colour_manual(values = rainbow(ncol(abundance))) +
#    theme_minimal() +
#    labs(title = "Species Abundance Distributions (Step Histogram)",
#         x = "Abundance", y = "Count") +
#    theme(legend.position = "none")





# CLOSE PDF
dev.off()
