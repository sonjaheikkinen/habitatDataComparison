# habitatDataComparison
A script for my thesis to see how well EU habitats directive habitat types predict bird species distribution in comparison to Corine land cover data.

Script was built with the help of the example pipeline available at [https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/past-hmsc-workshops](https://www.helsinki.fi/en/researchgroups/statistical-ecology/software/past-hmsc-workshops) 


## Running the scripts

1. Select your root directory in filepath_config.txt. Under this directory, create folders: "Data", "Models", "Results". Under models create folder "Fitted". Under results create folder "Modelfits". 
2. Select settings in main.R
    1. Select which scripts are run. Note: Automatic running of scripts not tested lately. Some of them might have broken unused parts, and need to be run manually to bypass those rows. 
        - run_preprocess_data: Cleans and filters data. 
        - run_format_data: Extracts data from files and formats it so that it can be utilized by Hmsc library
        - run_exploratory_analysis: Creates some graphs from raw and selected data. 
        - run_select_data: Filters out rare species and environmental variables that are left out from final model
        - run_define_models: Creates definitions for Hmsc models
        - run_fit_models: Fits the Hmsc models. Note: Can be very slow
        - run_check_model_convergence: Calculates convergence values for the fitted models
        - run_check_parameter_effects: Extracts and plots parameter values
        - run_check_model_fits: Runs model validation. Note: Even slower than fitting. 
    2. Select buffer width (how far away from the transect line are the habitat values extracted)
    3. Select number of samples per markov chain
    4. Select thinning value for markov chain
    5. Select how many groups the samples are divided in in cross-validation. Note: More groups --> slower validation
    6. Update your data filenames to match the ones you are using under comment "Read in rasters". 
3. Install necessary packages. Most libraries are listed in main, but there might be some used in analysis that are yet to be moved here.
4. Run main script to define, fit and validate models. Exploratory analysis might be better to to run manually
5. Run manually the other analysis scripts not configured for running trough main.R. Note: These are a work in progress. Might be easier to write your own analysis. 
    - compare_models: Creates some plots to compare model validation results
    - make_predictions: Calculates predictions, but their analysis and comparison is also currently here. Better to run manually, and skip the prediction creation rows, if you do not want to create new predictions every time, which takes some time. 

