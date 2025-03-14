#########################################################################
# Replication material for
# "Violence, Displacement, and Support for Internally Displaced Persons:
# Evidence from Syria"
# Authors: Alexandra Hartman, Benjamin Morse, Sigrid Weber
# Last update: 15/02/2021 
#########################################################################

Structure of the folder:
==========================
- data.csv = anonymized survey data from Syria (respondent ID and village are replaced to ensure anonymity)
- syria_shape = folder containing shapefiles for Syrian sub-districts
- codebook.xlsx = description of all variables in the data (+ those constructed through the analysis scripts)
- 1_data_cleaning.R = R script to clean the available survey data and prepare for analysis
- 2_data_analysis.R = R script to use the cleaned survey data for conjoint and observational analysis
- sensitivity_analysis = folder containing all necessary files to run the sensitivity analysis presented in the appendix of the paper


Steps necessary to replicate the findings:
==========================================
1. First run the 1_data_cleaning.R script which uses the data.csv to read the survey data, clean variables, construct indices, and reformat the data for the conjoint analysis. The script produces two datasets as output: "conjoint_data" (for conjoint analysis),"obs_data" (for observational analysis).

2. After cleaning the data, run the 2_data_analysis.R script. The script produces all models (conjoint & analysis) and replicates the main maps, figures, and tables in the paper in order of appearance. The script also replicates all findings from the online appendix. 

(3. The sensitivitiy analysis (following Oster) has been done with STATA. To replicate this analysis, access the sensitivity_analysis folder. The survey data is saved there in dta format (stata_file.dta). Run the Sensitivity_analysis.do file which requires the psacalc.do file to replicate the findings.)

 