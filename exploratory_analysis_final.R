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












# CLOSE PDF
dev.off()
