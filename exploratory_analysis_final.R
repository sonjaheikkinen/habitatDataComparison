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
load()


pdf(file.path(dir_results, sprintf("exploratory_analysis_final.pdf")))


# SPATIOTEMPORAL_CONTEXT OVERVIEW
years <- sort(unique(spatiotemporal_context$Year))
transects <- unique(spatiotemporal_context$Transect)
samples <- matrix(0, nrow = length(transects), ncol = length(years))
colnames(samples) <- years
rownames(samples) <- transects
for (sample in 1:nrow(spatiotemporal_context)) {
    year <- as.character(spatiotemporal_context[sample, ]$Year)
    transect <- spatiotemporal_context[sample, ]$Transect
    samples[transect, year] <- 1
}
pheatmap(samples,
         cluster_cols = FALSE)



# HABITAT TYPES OVERVIEW
plot(natura_raster)
corine_cropped <- crop(corine_raster, natura_raster)
plot(corine_cropped)
natura_frequencies <- freq(natura_raster)
natura_frequencies$fraction <- frequencies_to_fractions(natura_frequencies$count)
corine_frequencies <- freq(corine_cropped)
corine_frequencies$fraction <- frequencies_to_fractions(corine_frequencies$count)
old_par <- par(no.readonly = TRUE) 
par(mar = c(old_par$mar[1], 10, old_par$mar[3], old_par$mar[4]))
barplot(natura_frequencies$fraction,
        horiz = TRUE,
        las = 2,
        names.arg = natura_frequencies$value)
par(old_par)



# CLOSE PDF
dev.off()
