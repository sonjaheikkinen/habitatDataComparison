# THIS SCRIPT PREPROCESSES DATA: 
# The data is limited to only the study area
# The data is cleaned to remove problematic observations


# FUNCTIONS FOR PREPROCESSING

# SCRIPT STARTS
###################################################################################################################

# READ DATA
transects_shp <- vect(file.path(dir_data, "bird_transects.shp"))
observations <- read.csv(file.path(dir_data, "birdtransects_toRecbase.csv"))

# LIMIT DATA TO STUDY AREA

# CLEAN DATA

# SAVE PREPROCESSED DATA