%% run all: contains basic info all other scripts need. You can load the data in and run all analyses from here
%% general: data paths & subjects
paths{1}= PATH 1; %TODO fill in and adapt new_subjNames
paths{2}= PATH 2;
paths{3}= PATH 3;
paths{4}= PATH 4;
paths{5}= PATH 5;
paths{6}= PATH 6;
paths{7}= PATH 7;
paths{8}= PATH 8;

new_subjNames = {'S1','S2','S3','S4','S5','S6','S7','S8'}; % {'S3','S2','S6','S5','S4','S7','S8','S1'};

%% load in data ==> needs rmGet (and all the files that is dependent on); initHiddenGray; intersectCols
save_path = SAVE_PATH DATA; %TODO fill in

timing_maps = 1; % set however you like TM = 1, VFM = 0

if timing_maps == 1
    maps = 'timing_maps';
elseif timing_maps == 0
    maps = 'visual_field_maps';
end
if exist(strcat(save_path, '/ve_data_incl_vfm_',maps,'.mat'), 'file') ~= 2 % Only needs to be done if you don't already have ve_data
    ve_data = monoTuned_load_timing_data_anonymous_incl_vfminfo(new_subjNames, paths, timing_maps,save_path);
end
if exist(strcat(save_path, '/time_series_',maps,'.mat'), 'file') ~= 2 % Only needs to be done if you don't already have time_series
    time_series = monoTuned_load_time_series_anonymous(new_subjNames, paths, timing_maps,save_path); % never used, but nice to have
end

%% run analyses
minDuration=0.06;maxDuration=0.99;
minPeriod=0.06;maxPeriod=0.99;

save_path = SAVE_PATH RESULTS; %TODO fill in

    %% model comparisons (figure 3 and 7; suppl. figure 6 and 7)
medianPerSubj=1; % use 1, we take the medians per subject for model comparisons
minVE = 0.2;% use 0.2 for model comparisons
model_comparison = 1;

timing_maps = 1; % set however you like TM = 1, VFM = 0
monoTuned_model_comparisons_anonymous(save_path, minVE,timing_maps,medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod);
monoTuned_export_stats_anonymous(timing_maps, model_comparison, save_path)

    %% model parameter comparisons (figure 4)
medianPerSubj=1; % use 1, we take the medians per subject for model comparisons
minVE = 0.2;% use 0.2 for model parameter comparisons
timing_maps = 0; % set to 0 (VFM), because timing maps have already been done

monoTuned_model_parameters_anonymous(save_path, minVE,timing_maps,medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod);

    %% eccentricity comparisons (figure 5 and 6; suppl. figure 8)
medianPerSubj = 0; % use 0, we take the means per subject for ecc comparisons
minVE = 0;% use 0 for ecc comparisons
model_comparison = 0;
timing_maps = 0; % use 0: eccentricity is only compared in the VFM, because the timing maps do not include the entire range

monoTuned_ecc_comparison_anonymous(save_path, minVE, timing_maps,medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod);
monoTuned_export_stats_anonymous(timing_maps, model_comparison, save_path)


    %% mesh comparisons (figure 2; suppl. figure 5)
minVE = 0.1;% 0.2 for model comparisons, 0.1 for supplementary figure
monoTuned_mesh_comparisons_anonymous(new_subjNames, paths, save_path, minVE,minDuration,maxDuration,minPeriod,maxPeriod); % mrVista doesn't actually love this. If you get an error about  "datastructure"-something: run the mesh comparison file regularly (not as a function)

    %% data of separate voxels (figure 1; suppl. figure 1)
 % !!!NOTE!!! doesn't actually work as a function, because you need to debug another function within this
 % run the fitting figures file piece by piece from that file (and not as a function)
 % also, figures were saved manually

% Fig 1: whichVoxel = 19265; datapoint = 129; subj_id= 4; roi_name = 'VFMafni_leftsPCS2'; cv_run = 'TimingSweeps' 
% Suppl. Fig 1: whichVoxel = 360460; datapoint = 693; subj_id = 3; roi_name = 'VFMafni_rightV1'; cv_run = 'OddScans' 
%       and whichVoxel = 30126; datapoint = 354; subj_id = 3; roi_name = 'VFMafni_leftIPS2'; cv_run = 'OddScans' 

fitting_figures_anonymous(paths, subj_id, whichVoxel, datapoint, roi_name, cv_run);

    %% model validation (suppl. figure 2)
subj_id = 3; % we chose 3, you can chose anyone who has more voxels than you have output

% !!!NOTE!!! doesn't actually work as a function, because you need to debug another function within this
% run the Make Validation Predictions file piece by piece from that file (and not as a function)
save_path = SAVE_PATH DATA; %TODO fill in
if exist(strcat(save_path, '/test_models.mat'), 'file') ~= 2 % Only needs to be done if you don't already have test_models 
    MakeValidationPredictions_anonymous(paths, subj_id);
    
    % Now run crossValidate2Component to do the training of the models and crossvalidation
end

save_path = SAVE_PATH RESULTS; %TODO fill in
minVE = 0.2;% because this value used for model (parameter) comparisons
AnalysisValidationPredictions_anonymous(paths, subj_id, save_path,minVE, minDuration,maxDuration,minPeriod,maxPeriod);

