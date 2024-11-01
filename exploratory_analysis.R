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




# Close pdf
dev.off()
