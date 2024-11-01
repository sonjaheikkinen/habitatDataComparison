# EXPLORE RAW DATA

# Check that landscape is suitable for the metrics
#check_landscape(habitats) # Habitats ok!
#print(list_lsm(level = "landscape"), n = 66) # list all available metrics

# Open pdf for plots exploratory analysis
pdf(file.path(dir_results, "exploratory_analysis.pdf"))

# set filename for exploratory analysis text file
file <- file.path(dir_results, "exploratory_analysis.txt")


# LOAD RAW DATA
load(file = file.path(dir_data, "occurrence_raw.RData")) # occurrence
load(file = file.path(dir_data, "natura_raster.RData")) # natura_raster
load(file = file.path(dir_data, "fractions_natura.RData")) # fractions_natura



# BASIC NUMBERS ON RAW DATA


# Occurrence
number_of_species <- length(colnames(occurrence))
append_to_file(sprintf("Number of species: %s\n", number_of_species), file)
# Natura habitats
number_of_habitats_natura <- length(unique(natura_raster)[,1])
append_to_file(sprintf("Number of natura habitats: %s\n", number_of_habitats_natura), file)
number_of_known_habitats_within_buffer_natura <- length(colnames(fractions_natura))
append_to_file(sprintf("Number of known natura habitats within buffer: %s\n", 
                       number_of_known_habitats_within_buffer_natura), file)



# BASIC PLOTS FROM RAW DATA

# On how many years was each transect visited 
transect_visits <- aggregate(SampleEvent ~ Transect, data = observations,
                             FUN = function(x) length(unique(x)))
plot(as.factor(transect_visits$Transect), 
     transect_visits$SampleEvent, 
     type="l",
     main="Number of visits for each transect")
hist(transect_visits$SampleEvent, main = "Histogram for number of visits per transect")

# Plot amount of transects that was visited each year
year_transect_amounts <- aggregate(Transect ~ Year, data = observations,
                                   FUN = function(x) length(unique(x)))
plot(year_transect_amounts$Year, year_transect_amounts$Transect, type="h")
hist(year_transect_amounts$Transect, main = "Histogram for number of transects visited each year")

# Plot transect lengths
transect_names <- unique(observations
                         $Transect)
transect_lengths <- c()
for (transect in transect_names) {
    transect_lengths <- c(transect_lengths, mean(observations
                                                 [observations
                                                     $Transect == transect, "Effort"]))
}
transect_efforts <- data.frame(Transect = transect_names, Effort = transect_lengths)
hist(transect_efforts$Effort, main = "Histogram of transect lengths")






# Produce data
transect_natura_data_list <- list()
for (habitat_data_type in habitat_data_types) {
    natura_data_on_transects <- get_transect_habitat_data(buffer_width, 
                                                          natura_raster, 
                                                          transects_shp,
                                                          names_natura,
                                                          habitat_data_type)
    transect_natura_data_list[[habitat_data_type]] <- natura_data_on_transects
}
transect_species_data_list <- list()
for (species_data_type in species_data_types) {
    species_data_on_transects <- get_transect_species_data(species_traits,
                                                           observations,
                                                           species_data_type)
    transect_species_data_list[[species_data_type]] <- species_data_on_transects
}




# Close pdf
dev.off()
