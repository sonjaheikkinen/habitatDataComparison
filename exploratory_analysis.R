# EXPLORE RAW DATA

# Check that landscape is suitable for the metrics
#check_landscape(habitats) # Habitats ok!
#print(list_lsm(level = "landscape"), n = 66) # list all available metrics

# Open pdf for exploratory analysis
pdf(file.path(dir_results, "exploratory_analysis.pdf"))

dev.off()
