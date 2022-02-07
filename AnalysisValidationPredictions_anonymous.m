function AnalysisValidationPredictions_anonymous(paths, subj_id, save_path,minVE, minDuration,maxDuration,minPeriod,maxPeriod)
%% Analyzes model outputs on ground truth dataset
% paths: where the information about each participant is stored
% subj_id: index of subject you want (where in paths)
% save_path: general folder where you want results to be stored
% minVE: threshold for inclusion of the voxels
% minDuration,maxDuration,minPeriod,maxPeriod: ranges outside which the tuned model is considered invalid

save_path_stuct = strrep(save_path, [RESULTS],[DATA]); %TODO fill in
if exist(strcat(save_path_struct, '/parameters_validation.mat'), 'file') == 2
    load(strcat(save_path_struct, '/parameters_validation.mat'))
else
    exclude_outside_range=0;
    [actualValues, evalModels] = monoTuned_load_test_model_data_anonymous(paths, subj_id,save_path_stuct, minDuration,maxDuration,minPeriod,maxPeriod, exclude_outside_range);
    monoTuned_load_time_series_test_model_anonymous(paths, subj_id, save_path);
end
%% get all info

% Tuned dataset
X_tuned_tuned=[evalModels.Dtuned.Mtuned.X.Odd evalModels.Dtuned.Mtuned.X.Even];
Y_tuned_tuned=[evalModels.Dtuned.Mtuned.Y.Odd evalModels.Dtuned.Mtuned.Y.Even];
includedTuned=X_tuned_tuned > minDuration & X_tuned_tuned < maxDuration & Y_tuned_tuned > minPeriod & Y_tuned_tuned < maxPeriod;

% for proportion
VE_tuned_tuned=[evalModels.Dtuned.Mtuned.VE.Odd evalModels.Dtuned.Mtuned.VE.Even];
VE_tuned_tuned_incl=VE_tuned_tuned(includedTuned);
VE_tuned_tuned_excl=VE_tuned_tuned(~includedTuned);
VE_tuned_mono=[evalModels.Dtuned.Mmono.VE.Odd evalModels.Dtuned.Mmono.VE.Even];
VE_tuned_mono_incl=VE_tuned_mono(includedTuned);
VE_tuned_mono_excl=VE_tuned_mono(~includedTuned);

best_correct_incl_tuned=VE_tuned_tuned_incl>VE_tuned_mono_incl;
best_correct_excl_tuned=VE_tuned_tuned_excl>VE_tuned_mono_excl;

% for threshold
nonCV_VE_tuned_tuned=[evalModels.Dtuned.Mtuned.not_crossvalidated_VE.Odd evalModels.Dtuned.Mtuned.not_crossvalidated_VE.Even];
nonCV_VE_tuned_tuned_incl=nonCV_VE_tuned_tuned(includedTuned);
nonCV_VE_tuned_tuned_excl=nonCV_VE_tuned_tuned(~includedTuned);
nonCV_VE_tuned_mono=[evalModels.Dtuned.Mmono.not_crossvalidated_VE.Odd evalModels.Dtuned.Mmono.not_crossvalidated_VE.Even];
nonCV_VE_tuned_mono_incl=nonCV_VE_tuned_mono(includedTuned);
nonCV_VE_tuned_mono_excl=nonCV_VE_tuned_mono(~includedTuned);

nonCV_VE_tuned_best_incl=max([nonCV_VE_tuned_tuned_incl; nonCV_VE_tuned_mono_incl]);
nonCV_VE_tuned_best_excl=max([nonCV_VE_tuned_tuned_excl; nonCV_VE_tuned_mono_excl]);
threshold_id_incl_tuned = nonCV_VE_tuned_best_incl > -0.2;
threshold_id_excl_tuned = nonCV_VE_tuned_best_excl > -0.2;

% monotonic dataset
X_mono_tuned=[evalModels.Dmono.Mtuned.X.Odd evalModels.Dmono.Mtuned.X.Even];
Y_mono_tuned=[evalModels.Dmono.Mtuned.Y.Odd evalModels.Dmono.Mtuned.Y.Even];
includedMono=X_mono_tuned > minDuration & X_mono_tuned < maxDuration & Y_mono_tuned > minPeriod & Y_mono_tuned < maxPeriod;

% for proportion
VE_mono_mono=[evalModels.Dmono.Mmono.VE.Odd evalModels.Dmono.Mmono.VE.Even];
VE_mono_mono_incl=VE_mono_mono(includedMono);
VE_mono_mono_excl=VE_mono_mono(~includedMono);
VE_mono_tuned=[evalModels.Dmono.Mtuned.VE.Odd evalModels.Dmono.Mtuned.VE.Even];
VE_mono_tuned_incl=VE_mono_tuned(includedMono);
VE_mono_tuned_excl=VE_mono_tuned(~includedMono);

best_correct_incl_mono=VE_mono_mono_incl>VE_mono_tuned_incl;
best_correct_excl_mono=VE_mono_mono_excl>VE_mono_tuned_excl;

% for threshold
nonCV_VE_mono_mono=[evalModels.Dmono.Mmono.not_crossvalidated_VE.Odd evalModels.Dmono.Mmono.not_crossvalidated_VE.Even];
nonCV_VE_mono_mono_incl=nonCV_VE_mono_mono(includedMono);
nonCV_VE_mono_mono_excl=nonCV_VE_mono_mono(~includedMono);
nonCV_VE_mono_tuned=[evalModels.Dmono.Mtuned.not_crossvalidated_VE.Odd evalModels.Dmono.Mtuned.not_crossvalidated_VE.Even];
nonCV_VE_mono_tuned_incl=nonCV_VE_mono_tuned(includedMono);
nonCV_VE_mono_tuned_excl=nonCV_VE_mono_tuned(~includedMono);

nonCV_VE_mono_best_incl=max([nonCV_VE_mono_mono_incl; nonCV_VE_mono_tuned_incl]);
nonCV_VE_mono_best_excl=max([nonCV_VE_mono_mono_excl; nonCV_VE_mono_tuned_excl]);
threshold_id_incl_mono = nonCV_VE_mono_best_incl > -0.2;
threshold_id_excl_mono = nonCV_VE_mono_best_excl > -0.2;

%% make plots per noise bin
cd(strcat(save_path,'/test_models'))

noise_bins_tuned = unique(actualValues.Dtuned.Noise);
actual_noise_tuned = [actualValues.Dtuned.Noise' actualValues.Dtuned.Noise'];
noise_bins_mono = unique(actualValues.Dmono.Noise);
actual_noise_mono = [actualValues.Dmono.Noise' actualValues.Dmono.Noise'];


for bin = 1:length(noise_bins_tuned)
    noise_id_incl = actual_noise_tuned(includedTuned) == noise_bins_tuned(bin);
    combined_id_incl = noise_id_incl & threshold_id_incl_tuned;
    noise_id_excl = actual_noise_tuned(~includedTuned) == noise_bins_tuned(bin);
    combined_id_excl = noise_id_excl & threshold_id_excl_tuned;
    
    prop_incl_tuned(bin) = sum(best_correct_incl_tuned(combined_id_incl))/sum(combined_id_incl);
    prop_excl_tuned(bin) = sum(best_correct_excl_tuned(combined_id_excl))/sum(combined_id_excl);
    
end

figure; plot(noise_bins_tuned, prop_incl_tuned,'Color', [1,.56,.56]);
hold on; plot(noise_bins_tuned, prop_excl_tuned,'Color', [1,0,0]);


for bin = 1:length(noise_bins_mono)
    noise_id_incl = actual_noise_mono(includedMono) == noise_bins_mono(bin);
    combined_id_incl = noise_id_incl & threshold_id_incl_mono;
    noise_id_excl = actual_noise_mono(~includedMono) == noise_bins_mono(bin);
    combined_id_excl = noise_id_excl & threshold_id_excl_mono;
    
    prop_incl_mono(bin) = sum(best_correct_incl_mono(combined_id_incl))/sum(combined_id_incl);
    prop_excl_mono(bin) = sum(best_correct_excl_mono(combined_id_excl))/sum(combined_id_excl);
    
end

hold on; plot(noise_bins_mono, prop_incl_mono,'Color', [.56,.56,1]);
hold on; plot(noise_bins_mono, prop_excl_mono,'Color', [0,0,1]);

title(['Proportions correct dataset'])
xlabel('noise SDs')
ylabel('Proportions')
axis([0,6,0,1]);
axis square
legend([{'Included Tuned'}, {'Excluded Tuned'}, {'Included Mono'}, {'Excluded Mono'}])
export_fig('noise_bins.eps','-eps','-r600','-painters');
close all

%% Bin by best VE, separating by included and excluded
binsize = 0.05;

binVE=(binsize:binsize:1)-binsize/2;
for b=1:length(binVE)
    bii=nonCV_VE_tuned_best_incl>=binVE(b)-binsize/2 & nonCV_VE_tuned_best_incl<binVE(b)+binsize/2;
    proportionCorrectTuned(b)=sum(best_correct_incl_tuned(bii))/sum(bii);
    
    bii=nonCV_VE_mono_best_incl>=binVE(b)-binsize/2 & nonCV_VE_mono_best_incl<binVE(b)+binsize/2;
    proportionCorrectMono(b)=sum(best_correct_incl_mono(bii))/sum(bii);
    
    bii=nonCV_VE_tuned_best_excl>=binVE(b)-binsize/2 & nonCV_VE_tuned_best_excl<binVE(b)+binsize/2;
    proportionCorrectTunedExcluded(b)=sum(best_correct_excl_tuned(bii))/sum(bii);
    
    bii=nonCV_VE_mono_best_excl>=binVE(b)-binsize/2 & nonCV_VE_mono_best_excl<binVE(b)+binsize/2;
    proportionCorrectMonoExcluded(b)=sum(best_correct_excl_mono(bii))/sum(bii);
end
figure; plot(binVE, proportionCorrectMono,'Color', [.56,.56,1]);
hold on; plot(binVE, proportionCorrectTuned,'Color', [1,.56,.56]);
hold on; plot(binVE, proportionCorrectMonoExcluded,'Color', [0,0,1]);
hold on; plot(binVE, proportionCorrectTunedExcluded,'Color', [1,0,0]);
title(['Proportions correct dataset'])
xlabel('VE')
ylabel('Proportions')
axis([0,1,0,1]);
axis square
legend([{'Included Mono'}, {'Included Tuned'}, {'Excluded Mono'}, {'Excluded Tuned'}])
export_fig('VE_bins.eps','-eps','-r600','-painters');
close all

%% binning by actual variable
variable = 'sMaj'; %theta is more difficult, use round(..., 5)
variable = char(variable);
proportionCorrectIncluded = [];
proportionCorrectExcluded = [];

%1D plots for bins by variable
var_tuned=[actualValues.Dtuned.(variable)' actualValues.Dtuned.(variable)'];
var_tuned_incl = var_tuned(includedTuned);
var_tuned_excl = var_tuned(~includedTuned);
bins_incl=unique(round(var_tuned_incl,5));
bins_excl=unique(round(var_tuned_excl,5));
noise_id_incl = nonCV_VE_tuned_best_incl>minVE;
noise_id_excl = nonCV_VE_tuned_best_excl>minVE;

for b=1:length(bins_incl)
    bii = round(var_tuned_incl,5)==bins_incl(b) & noise_id_incl;
    proportionCorrectIncluded(b)=sum(best_correct_incl_tuned(bii))/sum(bii);
end

for b=1:length(bins_excl)
    bii = round(var_tuned_excl,5)==bins_excl(b) & noise_id_excl;
    proportionCorrectExcluded(b)=sum(best_correct_excl_tuned(bii))/sum(bii);
end

if string(variable) == 'Theta'
    figure; plot([rad2deg(bins_incl) 180], [proportionCorrectIncluded proportionCorrectIncluded(1)],'Color', [1,.56,.56]);
    hold on; plot([rad2deg(bins_excl) 180], [proportionCorrectExcluded proportionCorrectExcluded(1)],'Color', [1,0,0]);
    axis([0,180,0,1]);
else
    figure; plot(bins_incl, proportionCorrectIncluded,'Color', [1,.56,.56]);
    hold on; plot(bins_excl, proportionCorrectExcluded,'Color', [1,0,0]);
    axis([0,max(bins_excl),0,1]);
end

title(['Proportions correct per ' variable ' (threshold 0.2)'])
xlabel(variable)
ylabel('Proportion correct')
legend([{'Included Tuned'}, {'Excluded Tuned'}])
axis square

save_name = strcat(variable, '_bins_minVE=', str(minVE), '.eps');
export_fig(save_name,'-eps','-r600','-painters');
close all

end