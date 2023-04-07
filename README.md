# Instruction
This folder contains the R-files that conduct the simulations of the thesis.
In the following, it is listed which files generate which results.

## Results with generated data
Moderate sample size (Section 9.2.1): `RDD-with-Generated-Data.R`\
Large sample size (Section 9.2.1): `RDD-with-Generated-Data.R`\
Redundant covariates (Section 9.3): `RDD-with-Generated-Data_Redundant-Covariates.R`

## Results with data set
Data set on Austrian labor market (Section 9.4): `RDD-with-Dataset.R`

## Requirements
The requirements such that the above R-files run are:
* R must be installed
* All required packages must be installed (`utils`, `haven`, `RDRobust`, `RDHonest`, `mvtnorm`, `parallel`)
* The working directory must be set to this folder
* This folder must contain another folder `data` containing the file `Card_analysis_dataset.dta`
* This folder must contain another folder `R` containing the files `functions.R` and `RDD_functions.R`

Those requirements are given in each R-file in more details as a comment at the beginning.
In fact, the installation of the required packages as well as the configuration of the working directory
should be done automatically within these R-files. Only in the case of unexpected incompatibilities, there might
be a need to do this manually.
Moreover, the folder `R` contains a documentation of all functions included in it.