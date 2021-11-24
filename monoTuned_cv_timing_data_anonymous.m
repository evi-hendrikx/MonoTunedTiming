function cv_data = monoTuned_cv_timing_data_anonymous(ve_data, use_models)
%% cv timing data: builds new dataset including only cross validated datapoints
% basically kicks out TimingSweepsL1 of ve_data and rearranges order of selected models if necessary
%   ve_data: struct with voxel information (e.g., roi, variance explained, etc.)
%   use_models: allows selection of 2 of the included response models in ve_data. Our main analyses were done on MonoOcc (the best-performing monotonic model) and TunedLin2d


modelFieldNames = use_models;
subjNames = fieldnames(ve_data.(modelFieldNames{1}));

for subj = 1:length(subjNames)
    condNames = fieldnames(ve_data.(modelFieldNames{1}).(subjNames{subj}));
    
    % VERY IMPORTANT: makes sure the tuned model (if it's included) is always
    % second, which is what the next scripts assume
    if find(strcmp(modelFieldNames,'TunedLin2d')) == 1
        model_order = length(modelFieldNames):-1:1;
    else
        model_order = 1:length(modelFieldNames);
    end
    
    for models = 1:length(modelFieldNames)
        cv_data.(modelFieldNames{model_order(models)}).(subjNames{subj}).(condNames{2}) = ve_data.(modelFieldNames{model_order(models)}).(subjNames{subj}).(condNames{2});
        cv_data.(modelFieldNames{model_order(models)}).(subjNames{subj}).(condNames{3}) = ve_data.(modelFieldNames{model_order(models)}).(subjNames{subj}).(condNames{3});
    end
end
