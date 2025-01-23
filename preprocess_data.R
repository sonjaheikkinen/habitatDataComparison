# THIS SCRIPT PREPROCESSES DATA: 
# The data is limited to only the study area
# The data is cleaned to remove problematic observations


# FUNCTIONS FOR PREPROCESSING

get_transects_that_overlap_raster <- function(transects, raster) {
    raster_values_for_transects <- extract(raster, transects)
    transects_overlapping_raster <- c()
    for (transect_id in unique(raster_values_for_transects$ID)) {
        values_for_transect <- raster_values_for_transects[raster_values_for_transects$ID == transect_id,]$natura_2393
        if (!all(is.na(values_for_transect))) {
            transects_overlapping_raster <- c(transects_overlapping_raster, transect_id)
        }
    }
    filtered_transects <- transects[transects$ID %in% transects_overlapping_raster]
    return(filtered_transects)
}

get_observations_from_given_transects <- function(transects, observations) {
    transects_to_keep <- transects$Numero
    filtered_observations <- observations[observations$Transect %in% transects_to_keep,]
    
    return(filtered_observations)
}

# SCRIPT STARTS
###################################################################################################################

# READ DATA
transects_shp <- vect(file.path(dir_data, "bird_transects.shp"))
observations <- read.csv(file.path(dir_data, "birdtransects_toRecbase.csv"))
id_mapping <- read.csv(file.path(dir_data, "All_GIS_Routes_LinjaRouteID_Matched.csv"), sep = ";")


# LIMIT DATA TO STUDY AREA

# Limit transects to study area
transects_shp$ID <- 1:nrow(transects_shp) # Adding numerical id, because existing chr identifier caused issues
transects_that_overlap_natura_raster <- get_transects_that_overlap_raster(transects_shp, natura_raster)
transects <- transects_that_overlap_natura_raster

# Limit observations to study area
# Transects are numbered differently in observations and transects shp
# First the transect numbers from shp need to be added for observations
# This is done using the id mapping file which contains the matching numbers from both numbering systems
transect_numbers_in_observations_data <- observations$SiteID
transect_number_list_observation <- id_mapping$RouteID_REC
transect_number_list_transect_shp <- id_mapping$GIS_route
transect_number_indices_in_id_mapping <- match(transect_numbers_in_observations_data, 
                                               transect_number_list_observation)
observations$Transect <- transect_number_list_transect_shp[transect_number_indices_in_id_mapping]
# Remove observations that are from transects that do not overlap natura data
observations_from_study_area <- get_observations_from_given_transects(transects, observations)
observations <- observations_from_study_area


# LIMIT DATA TO STUDY YEARS
observations_2006_onwards <- observations$Year >= 2006
observations <- observations[observations_2006_onwards, ]



# CLEAN DATA

# Remove species that are not birds
not_birds <- c("Myodes rufocanus")
observations <- observations[!observations$Species %in% not_birds,]
#Remove species that do not have full species names TO DO: could these somehow be incorporated?
observations <- observations[sapply(strsplit(observations$Species, " "), length) == 2, ]
# Filter out transects that do not have any observations
transect_numbers_that_have_observations <- unique(observations$Transect)
transects <- subset(transects,
                    transects$Numero %in% transect_numbers_that_have_observations)


# SAVE PREPROCESSED DATA
# Save preprocessed data
writeVector(transects, 
            file.path(dir_data, "transects_preprocessed.shp"),
            overwrite = TRUE)
write.table(observations, 
            file.path(dir_data, "observations_preprocessed.csv"), 
            sep = ";",
            row.names = FALSE)