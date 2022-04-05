# MonoTunedTiming

These scripts were used to run the analysis for 

Hendrikx, E., Paul, J.M., van Ackooij, M., van der Stoep, N., & Harvey, B.M. (2021) Visual timing-tuned responses in human association cortices and response dynamics in early visual cortex. Nature Communications (In revision). 

The scripts were written in Matlab 2017a.

## Running the pipeline
### creating the dataset
monoTuned_run_all_anonymous basically calls all other scripts you need. You still need to provide paths where the subjects' information is stored, as well as the path where you want your results to be saved. It initially creates datasets with the required parameters for all the analyses (monoTuned_load_timing_data_anonymous_incl_vfminfo and monoTuned_load_time_series_anonymous).
Note that the datasets created here for our data are available at the following DOIs: visual field map parameters (doi.org/10.6084/m9.figshare.19146131), timing map parameters (doi.org/10.6084/m9.figshare.17122706), visual field map time series (doi.org/10.6084/m9.figshare.19146092) & timing map time series (doi.org/10.6084/m9.figshare.17122718) 
It then runs the model comparisons (monoTuned_model_comparisons_anonymous), model parameter comparisons (monoTuned_model_parameters_anonymous), and eccentricity comparisons (monoTuned_ecc_comparison_anonymous) using these datasets.

### Creating meshes & voxel data
As the meshes require 3d information about the voxel coordinates, this information cannot be extracted for the previously created datasets. They are created separately (monoTuned_mesh_comparisons_anonymous)
Furthermore, the time series and model predictions for a single voxel cannot be extracted from summary statistics (like variance explained of a voxel), which is why this script (fitting_figures_anonymous) is also independent from the previously created datasets.

### Model validation
In order to validate the models, data with a known ground-truth was simulated (MakeValidationPredictions_anonymous) and models were run to predict values for this simulated dataset (crossval2component and dependencies available in other github repositories: vistasoft (github.com/vistalab/vistasoft); vistasoftAddOns (github.com/benharvey/vistasoftAddOns); fMRI_preproc (github.com/MvaOosterhuis/fMRI_preproc )). The results were then analyzed (AnalysisValidationPredictions_anonymous).
Note that the simulated data and the predictions of the models are available at the following DOIs: parameters used during validation (doi.org/10.6084/m9.figshare.17122727) & validation time series (doi.org/10.6084/m9.figshare.17122748)

## Licence
Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a
[Creative Commons Attribution 4.0 International License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg

## Contact
For any questions or remarks, please send an email to b.m.harvey@uu.nl
