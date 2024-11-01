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

# BASIC NUMBERS ON RAW DATA
# Occurrence
number_of_species <- length(colnames(occurrence))
append_to_file(sprintf("Number of species: %s\n", number_of_species), file)
# Natura raster
number_of_habitats_natura <- length(unique(habitat_types_raster)[,1])
append_to_file(sprintf("Number of natura habitats: %s\n", number_of_habitats_natura), file)




# Close pdf
dev.off()
