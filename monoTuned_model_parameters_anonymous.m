function monoTuned_model_parameters_anonymous(save_path, minVE, timing_maps,medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod)
%% model_parameters: compares the parameters of a response models between various ROIs
% Compares the summary statistic of parameters 
% Runs ANOVA (subj, roi), post-hocs (Dunn), and gives graphs for figure 4
%
% save_path: general folder where you want results to be stored
% minVE: threshold for inclusion of the voxels
% timing maps: 1 if analyze timing maps; 0 if analyze visual field maps
% medianPerSubj: determines the summary statistic: 1 = means; 2 = medians; 3 = 75th percentile; 4 = 90th percentile
% minDuration,maxDuration,minPeriod,maxPeriod: ranges outside which the tuned model is considered invalid


%% load and prepare data
if timing_maps == 1
    load('ve_data_timing_maps.mat')
elseif timing_maps == 0
    load('ve_data_visual_field_maps.mat')
end

save_path_model_properties = strcat(save_path, 'model_property_comparisons/');

%choose which models you want to use
modelFieldNames = fieldnames(ve_data);
use_models = {modelFieldNames{1},modelFieldNames{2}}; % fill in models you want

cv_data = monoTuned_cv_timing_data_anonymous(ve_data, use_models);
modelFieldNames = fieldnames(cv_data);
subjNames = fieldnames(cv_data.(modelFieldNames{1}));
condNames = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{1}));
ROILabels = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{1}).(condNames{1}).crossValidated); %possible, because S3 has all maps
ROIs = unique(erase(ROILabels, ["Right","Left","right","left"]),'stable');

%% prepare data for ANOVA
hemispheres = {'left','right'};
meanSub = [];
ANOVAprop = [];

% get mean values per roi per run, put in struct with runs and hemispheres
% in same column
for subj = 1:length(subjNames)
    for roi = 1:length(ROILabels)
        
        if roi <= length(ROIs) % left
            roi_id = roi;
            hemisphere = 1;
        else
            roi_id = roi - length(ROIs); % right
            hemisphere = 2;
        end
        
        roi_names_subj = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{subj}).(condNames{1}).crossValidated);
        for condition = 1:length(condNames)
            
            if any(strcmp(roi_names_subj, ROILabels{roi}))
                
                %NOTE: if tuned model is involved, it's always the second model (due to build cv data)      
                % variance explained on not cross-validated % necessary for voxel selection
                ve_mod1_Current_pre_cv = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).used_for_crossValidation.(ROILabels{roi}).varianceExplained';
                ve_mod2_Current_pre_cv = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).used_for_crossValidation.(ROILabels{roi}).varianceExplained';
                
                % compressive exponents (and if second model is tuned preferred duration/period)
                xs_mod1_Current = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).xs';
                ys_mod1_Current = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).ys';
                xs_mod2_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).xs';
                ys_mod2_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).ys';
                
                % beta mono
                beta_xs_mod1_Current = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).used_for_crossValidation.(ROILabels{roi}).beta1;
                beta_ys_mod1_Current = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).used_for_crossValidation.(ROILabels{roi}).beta2;
                
                % if parameters out of range for the tuned model: set VE to 0
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    sigma_major_tuned_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).sigma_major';
                    sigma_minor_tuned_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).sigma_minor';
                    theta_tuned_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).theta';
                    exp_tuned_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).exponent';
                    
                    % NOTE: these voxels are later excluded, because 0< minVE (IF YOU CHANGE MINVE CHECK THAT THIS STILL HAPPENS)
                    outsideRange = find(sum([xs_mod2_Current<=minDuration, xs_mod2_Current>=maxDuration,ys_mod2_Current<=minPeriod,ys_mod2_Current>=maxPeriod],2)>0);
                    ve_mod2_Current_pre_cv(outsideRange)=0;
                    
                else % if second model is also monotonic
                    beta_xs_mod2_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).beta1;
                    beta_ys_mod2_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).beta2;
                end
                
                
                % NOTE: different model selection than in the other analyses. Here we
                % want voxels selected for model 1 that do well for model 1 (>minVE),
                % and select voxels for model 2 that do well for model 2
                % (>minVE) separately
                % voxel selection at subject level, remove if < minVE for uncrossvalidated data
                removeVE_mod1 = find(ve_mod1_Current_pre_cv<=minVE);
                xs_mod1_Current(removeVE_mod1) = [];
                ys_mod1_Current(removeVE_mod1) = [];
                
                % these (the "beta" parameters) have a slightly different voxel selection and
                % therefore sometimes iCrds is slightly not the same size
                % as the "beta" values
                rmCrds1 = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).iCoords(removeVE_mod1);
                betaCrds = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).beta_Coords;
                [~,~,rmBetaCrds1] = intersect(rmCrds1,betaCrds);
                
                % remove voxel from analysis if ratio cannot properly be
                % calculated (i.e., at least one of the betas = 0)
                beta_xs_mod1_Current(rmBetaCrds1,:) = []; % so also works for S3 and S2
                beta_ys_mod1_Current(rmBetaCrds1,:) = [];% so also works for S3 and S2
                
                beta_ratio_mod1_Current = beta_xs_mod1_Current./beta_ys_mod1_Current;
                beta_ratio_mod1_Current(beta_ratio_mod1_Current ==0)=NaN;
                beta_ratio_mod1_Current(beta_ratio_mod1_Current ==Inf)=NaN;
                tmp = nanmean(beta_ratio_mod1_Current,2); % so also works for S3 and S2
                beta_ratio_mod1_Current = tmp(~isnan(tmp));
                
                removeVE_mod2 = find(ve_mod2_Current_pre_cv<= minVE);
                xs_mod2_Current(removeVE_mod2) = [];
                ys_mod2_Current(removeVE_mod2) = [];
                
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    sigma_major_tuned_Current(removeVE_mod2)=[];
                    sigma_minor_tuned_Current(removeVE_mod2)=[];
                    theta_tuned_Current(removeVE_mod2)=[];
                    exp_tuned_Current(removeVE_mod2)=[];
                else
                    rmCrds2 = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).iCoords(removeVE_mod2);
                    betaCrds = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).beta_Coords;
                    [~,~,rmBetaCrds2] = intersect(rmCrds2,betaCrds);
                    
                    beta_xs_mod2_Current(rmBetaCrds2,:) = [];
                    beta_ys_mod2_Current(rmBetaCrds2,:) = [];
                    
                    beta_ratio_mod2_Current = beta_xs_mod2_Current./beta_ys_mod2_Current;
                    beta_ratio_mod2_Current(beta_ratio_mod2_Current ==0)=NaN;
                    beta_ratio_mod2_Current(beta_ratio_mod2_Current ==Inf)=NaN;
                    tmp = nanmean(beta_ratio_mod2_Current,2); % so also works for S3 and S2
                    beta_ratio_mod2_Current = tmp(~isnan(tmp));
                end
                
            else
                xs_mod1_Current = [];
                ys_mod1_Current = [];
                xs_mod2_Current = [];
                ys_mod2_Current = [];
                beta_xs_mod1_Current = [];
                beta_ys_mod1_Current = [];
                beta_ratio_mod1_Current = [];
                
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    sigma_major_tuned_Current=[];
                    sigma_minor_tuned_Current=[];
                    theta_tuned_Current=[];
                    exp_tuned_Current=[];
                else
                    beta_xs_mod2_Current = [];
                    beta_ys_mod2_Current = [];
                    beta_ratio_mod2_Current = [];
                end
                
            end
            
            %% calculate one value per hemisphere            
            if medianPerSubj == 0  
                meanSub.xs.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmean(xs_mod1_Current);
                meanSub.xs.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(xs_mod2_Current);
                meanSub.ys.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmean(ys_mod1_Current);
                meanSub.ys.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(ys_mod2_Current);
                
                meanSub.beta_xs.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmean(mean(beta_xs_mod1_Current,2));% so also works for S3 and S2
                meanSub.beta_ys.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmean(mean(beta_ys_mod1_Current,2));% so also works for S3 and S2
                meanSub.beta_ratio.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmean(beta_ratio_mod1_Current);
                
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    meanSub.major.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(sigma_major_tuned_Current);
                    meanSub.minor.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(sigma_minor_tuned_Current);
                    meanSub.theta_rad.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(theta_tuned_Current);
                    meanSub.theta.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(rad2deg(theta_tuned_Current));
                    meanSub.exp.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(exp_tuned_Current);
                    meanSub.ratio.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(sigma_major_tuned_Current./sigma_minor_tuned_Current);
                else
                    meanSub.beta_xs.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(mean(beta_xs_mod2_Current,2));% so also works for S3 and S2
                    meanSub.beta_ys.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(mean(beta_ys_mod2_Current,2));% so also works for S3 and S2
                    meanSub.beta_ratio.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(beta_ratio_mod2_Current);
                end
            elseif medianPerSubj == 1
                meanSub.xs.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmedian(xs_mod1_Current);
                meanSub.xs.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(xs_mod2_Current);
                meanSub.ys.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmedian(ys_mod1_Current);
                meanSub.ys.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(ys_mod2_Current);
                
                meanSub.beta_xs.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmedian(mean(beta_xs_mod1_Current,2));% so also works for S3 and S2
                meanSub.beta_ys.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmedian(mean(beta_ys_mod1_Current,2));% so also works for S3 and S2
                meanSub.beta_ratio.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmedian(beta_ratio_mod1_Current);
                
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    meanSub.major.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(sigma_major_tuned_Current);
                    meanSub.minor.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(sigma_minor_tuned_Current);
                    meanSub.theta_rad.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(theta_tuned_Current);
                    meanSub.theta.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(rad2deg(theta_tuned_Current));
                    meanSub.exp.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(exp_tuned_Current);
                    meanSub.ratio.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(sigma_major_tuned_Current./sigma_minor_tuned_Current);
                else
                    meanSub.beta_xs.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(mean(beta_xs_mod2_Current,2));% so also works for S3 and S2
                    meanSub.beta_ys.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(mean(beta_ys_mod2_Current,2));% so also works for S3 and S2
                    meanSub.beta_ratio.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(beta_ratio_mod2_Current);
                end
            elseif medianPerSubj == 2  
                meanSub.xs.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(xs_mod1_Current,75);
                meanSub.xs.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(xs_mod2_Current,75);
                meanSub.ys.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(ys_mod1_Current,75);
                meanSub.ys.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(ys_mod2_Current,75);
                
                meanSub.beta_xs.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(mean(beta_xs_mod1_Current,2),75);% so also works for S3 and S2
                meanSub.beta_ys.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(mean(beta_ys_mod1_Current,2),75);% so also works for S3 and S2
                meanSub.beta_ratio.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(beta_ratio_mod1_Current,75);
                
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    meanSub.major.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(sigma_major_tuned_Current,75);
                    meanSub.minor.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(sigma_minor_tuned_Current,75);
                    meanSub.theta_rad.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(theta_tuned_Current,75);
                    meanSub.theta.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(rad2deg(theta_tuned_Current),75);
                    meanSub.exp.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(exp_tuned_Current,75);
                    meanSub.ratio.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(sigma_major_tuned_Current./sigma_minor_tuned_Current,75);
                else
                    meanSub.beta_xs.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(mean(beta_xs_mod2_Current,2),75);% so also works for S3 and S2
                    meanSub.beta_ys.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(mean(beta_ys_mod2_Current,2),75);% so also works for S3 and S2
                    meanSub.beta_ratio.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(beta_ratio_mod2_Current,75);
                end
            elseif medianPerSubj == 3            
                meanSub.xs.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(xs_mod1_Current,90);
                meanSub.xs.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(xs_mod2_Current,90);
                meanSub.ys.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(ys_mod1_Current,90);
                meanSub.ys.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(ys_mod2_Current,90);
                
                meanSub.beta_xs.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(mean(beta_xs_mod1_Current,2),90);% so also works for S3 and S2
                meanSub.beta_ys.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(mean(beta_ys_mod1_Current,2),90);% so also works for S3 and S2
                meanSub.beta_ratio.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(beta_ratio_mod1_Current,90);
                
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    meanSub.major.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(sigma_major_tuned_Current,90);
                    meanSub.minor.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(sigma_minor_tuned_Current,90);
                    meanSub.theta_rad.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(theta_tuned_Current,90);
                    meanSub.theta.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(rad2deg(theta_tuned_Current),90);
                    meanSub.exp.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(exp_tuned_Current,90);
                    meanSub.ratio.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(sigma_major_tuned_Current./sigma_minor_tuned_Current,90);
                else
                    meanSub.beta_xs.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(mean(beta_xs_mod2_Current,2),90);% so also works for S3 and S2
                    meanSub.beta_ys.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(mean(beta_ys_mod2_Current,2),90);% so also works for S3 and S2
                    meanSub.beta_ratio.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(beta_ratio_mod2_Current,90);
                end
            end
            
            meanSubSubj(subj,hemisphere,condition,roi_id) = subjNames(subj);
            meanSubRoi(subj,hemisphere,condition,roi_id) = ROIs(roi_id);
            meanSubHemi(subj,hemisphere,condition,roi_id) = hemispheres(hemisphere);
            meanSubCond(subj,hemisphere,condition,roi_id) = condNames(condition);
            
        end
    end
end

%% get info into new struct
ANOVAprop.roi.(modelFieldNames{1}) = meanSubRoi(:);
ANOVAprop.subj.(modelFieldNames{1}) = meanSubSubj(:);
ANOVAprop.roi.(modelFieldNames{2}) = meanSubRoi(:);
ANOVAprop.subj.(modelFieldNames{2}) = meanSubSubj(:);

ANOVAprop.xs.(modelFieldNames{1}) = meanSub.xs.(modelFieldNames{1})(:);
ANOVAprop.ys.(modelFieldNames{1}) = meanSub.ys.(modelFieldNames{1})(:);
ANOVAprop.xs.(modelFieldNames{2}) = meanSub.xs.(modelFieldNames{2})(:);
ANOVAprop.ys.(modelFieldNames{2}) = meanSub.ys.(modelFieldNames{2})(:);

ANOVAprop.beta_ys.(modelFieldNames{1}) = meanSub.beta_ys.(modelFieldNames{1})(:);
ANOVAprop.beta_xs.(modelFieldNames{1}) = meanSub.beta_xs.(modelFieldNames{1})(:);
ANOVAprop.beta_ratio.(modelFieldNames{1}) = meanSub.beta_ratio.(modelFieldNames{1})(:);

% different removal voxels betas
ANOVAprop.roi_beta.(modelFieldNames{1}) = meanSubRoi(:); 
ANOVAprop.subj_beta.(modelFieldNames{1}) = meanSubSubj(:);

if find(strcmp(modelFieldNames,'TunedLin2d'))
    ANOVAprop.major.(modelFieldNames{2})= meanSub.major.(modelFieldNames{2})(:);
    ANOVAprop.minor.(modelFieldNames{2})= meanSub.minor.(modelFieldNames{2})(:);
    ANOVAprop.theta_rad.(modelFieldNames{2})= meanSub.theta_rad.(modelFieldNames{2})(:);
    ANOVAprop.theta.(modelFieldNames{2})= meanSub.theta.(modelFieldNames{2})(:);
    ANOVAprop.exp.(modelFieldNames{2})= meanSub.exp.(modelFieldNames{2})(:);
    ANOVAprop.ratio.(modelFieldNames{2})= meanSub.ratio.(modelFieldNames{2})(:);
else
    ANOVAprop.roi_beta.(modelFieldNames{2}) = meanSubRoi(:);% different removal voxels
    ANOVAprop.subj_beta.(modelFieldNames{2}) = meanSubSubj(:);% different removal voxels
    
    ANOVAprop.beta_ys.(modelFieldNames{2}) = meanSub.beta_ys.(modelFieldNames{2})(:);
    ANOVAprop.beta_xs.(modelFieldNames{2}) = meanSub.beta_xs.(modelFieldNames{2})(:);
    ANOVAprop.beta_ratio.(modelFieldNames{2}) = meanSub.beta_ratio.(modelFieldNames{2})(:);
end

%% change format and remove nan values from new struct and 
% if voxel is missing for one model in one parameter, it's missing for the
% others too (due to exclusion with minVE threshold of this voxel)
    % exception: beta values, which include slightly different voxels and
    % exclude invalid ratio-voxels
ANOVAprop.roi.(modelFieldNames{1})(isnan(meanSub.xs.(modelFieldNames{1})(:))) = []; 
ANOVAprop.subj.(modelFieldNames{1})(isnan(meanSub.xs.(modelFieldNames{1})(:))) = []; 
ANOVAprop.roi.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];
ANOVAprop.subj.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];

ANOVAprop.xs.(modelFieldNames{1})(isnan(meanSub.xs.(modelFieldNames{1})(:))) = []; 
ANOVAprop.ys.(modelFieldNames{1})(isnan(meanSub.xs.(modelFieldNames{1})(:))) = [];
ANOVAprop.xs.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];
ANOVAprop.ys.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];

ANOVAprop.beta_ratio.(modelFieldNames{1})(isnan(meanSub.beta_ratio.(modelFieldNames{1})(:))) = [];
ANOVAprop.roi_beta.(modelFieldNames{1})(isnan(meanSub.beta_ratio.(modelFieldNames{1})(:))) = [];
ANOVAprop.subj_beta.(modelFieldNames{1})(isnan(meanSub.beta_ratio.(modelFieldNames{1})(:))) = [];

if find(strcmp(modelFieldNames,'TunedLin2d'))
    ANOVAprop.major.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];
    ANOVAprop.minor.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];
    ANOVAprop.theta_rad.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];
    ANOVAprop.theta.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];
    ANOVAprop.exp.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];
    ANOVAprop.ratio.(modelFieldNames{2})(isnan(meanSub.xs.(modelFieldNames{2})(:))) = [];
else
    ANOVAprop.beta_ratio.(modelFieldNames{2})(isnan(meanSub.beta_ratio.(modelFieldNames{2})(:))) = [];
    ANOVAprop.roi_beta.(modelFieldNames{2})(isnan(meanSub.beta_ratio.(modelFieldNames{2})(:))) = [];
    ANOVAprop.subj_beta.(modelFieldNames{2})(isnan(meanSub.beta_ratio.(modelFieldNames{2})(:))) = [];
end


savename = 'ANOVAprop_model_params_';
if medianPerSubj == 1
    savename = strcat(savename, 'medianPerSubj');
elseif medianPerSubj == 0
    savename = strcat(savename, 'meanPerSubj');
elseif medianPerSubj == 2
    savename = strcat(savename, '75thPerSubj');
elseif medianPerSubj == 3
    savename = strcat(savename, '90thPerSubj');
end

if timing_maps == 1
    savename = strcat(savename, '_timing_maps');
elseif timing_maps == 0
    savename = strcat(savename, '_visual_field_maps');
end
savename = strcat(savename, '_minVE=', string(minVE), '.mat');
save(savename, 'ANOVAprop')


%% reshape meansub because later scripts expect all data points per roi under each other
meanSub.xs.(modelFieldNames{1}) = reshape(meanSub.xs.(modelFieldNames{1}),length(subjNames)*length(condNames)*2,length(ROIs));
meanSub.ys.(modelFieldNames{1}) = reshape(meanSub.ys.(modelFieldNames{1}),length(subjNames)*length(condNames)*2,length(ROIs));
meanSub.xs.(modelFieldNames{2}) = reshape(meanSub.xs.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
meanSub.ys.(modelFieldNames{2}) = reshape(meanSub.ys.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));

meanSub.beta_xs.(modelFieldNames{1}) = reshape(meanSub.beta_xs.(modelFieldNames{1}),length(subjNames)*length(condNames)*2,length(ROIs));
meanSub.beta_ys.(modelFieldNames{1}) = reshape(meanSub.beta_ys.(modelFieldNames{1}),length(subjNames)*length(condNames)*2,length(ROIs));
meanSub.beta_ratio.(modelFieldNames{1}) = reshape(meanSub.beta_ratio.(modelFieldNames{1}),length(subjNames)*length(condNames)*2,length(ROIs));
if find(strcmp(modelFieldNames,'TunedLin2d'))
    meanSub.major.(modelFieldNames{2}) = reshape(meanSub.major.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.minor.(modelFieldNames{2}) = reshape(meanSub.minor.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.theta_rad.(modelFieldNames{2}) = reshape(meanSub.theta_rad.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.theta.(modelFieldNames{2}) = reshape(meanSub.theta.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.exp.(modelFieldNames{2}) = reshape(meanSub.exp.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.ratio.(modelFieldNames{2}) = reshape(meanSub.ratio.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
else
    meanSub.beta_xs.(modelFieldNames{2}) = reshape(meanSub.beta_xs.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.beta_ys.(modelFieldNames{2}) = reshape(meanSub.beta_ys.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.beta_ratio.(modelFieldNames{2}) = reshape(meanSub.beta_ratio.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
end

meanSub.subj = reshape(meanSubSubj,length(subjNames)*length(condNames)*2,length(ROIs));
meanSub.roi = reshape(meanSubRoi,length(subjNames)*length(condNames)*2,length(ROIs));
meanSub.hemi = reshape(meanSubHemi,length(subjNames)*length(condNames)*2,length(ROIs));
meanSub.cond = reshape(meanSubCond,length(subjNames)*length(condNames)*2,length(ROIs));

cd(save_path_model_properties)

savename = 'meanSub_model_params_';
if medianPerSubj == 1
    savename = strcat(savename, 'medianPerSubj');
elseif medianPerSubj == 0
    savename = strcat(savename, 'meanPerSubj');
elseif medianPerSubj == 2
    savename = strcat(savename, '75thPerSubj');
elseif medianPerSubj == 3
    savename = strcat(savename, '90thPerSubj');
end

if timing_maps == 1
    savename = strcat(savename, '_timing_maps');
elseif timing_maps == 0
    savename = strcat(savename, '_visual_field_maps');
end
savename = strcat(savename, '_minVE=', string(minVE), '.mat');
save(savename, 'meanSub')

%% do the ANOVA and post-hocs
terms(1,1) = 1; % main subj
terms(2,2) = 1; % main map
stat = [];

if find(strcmp(modelFieldNames,'TunedLin2d'))
    parameters = {'xs','ys', 'beta_ratio','major','minor','theta','exp','ratio'};
else
    parameters = {'xs','ys','beta_ratio'};
end

counter = 25; % random number, higher than figures already on the screen
for parameter = 1:length(parameters)
    
    if parameter == 1 || parameter == 2
        modelNames = {modelFieldNames{1}, modelFieldNames{2}};
    elseif parameter == 3
        modelNames = {modelFieldNames{1}};
    elseif parameter > 3 && find(strcmp(modelFieldNames,'TunedLin2d'))
        % some analyses can only be done for tuned
        modelNames = {modelFieldNames{2}};
    end
    
    for model = 1:length(modelNames)
        
        if parameter == 3
            [stat.(parameters{parameter}).(modelNames{model}).p,stat.(parameters{parameter}).(modelNames{model}).tbl,stat.(parameters{parameter}).(modelNames{model}).stats,stat.(parameters{parameter}).(modelNames{model}).terms] = anovan(ANOVAprop.(parameters{parameter}).(modelNames{model}),{ANOVAprop.subj_beta.(modelNames{model}),ANOVAprop.roi_beta.(modelNames{model})},...
                'model',terms, 'varnames', {'subject', 'map'});
        else
            [stat.(parameters{parameter}).(modelNames{model}).p,stat.(parameters{parameter}).(modelNames{model}).tbl,stat.(parameters{parameter}).(modelNames{model}).stats,stat.(parameters{parameter}).(modelNames{model}).terms] = anovan(ANOVAprop.(parameters{parameter}).(modelNames{model}),{ANOVAprop.subj.(modelNames{model}),ANOVAprop.roi.(modelNames{model})},...
                'model',terms, 'varnames', {'subject', 'map'});
        end
        
        % for paired comparisons: load data without NaNs
        data_rois = [];
        nr_rois = [];
        for roi = 1:length(ROIs)
            stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi}) = meanSub.(parameters{parameter}).(modelNames{model})(~isnan(meanSub.(parameters{parameter}).(modelNames{model})(:,roi)),roi);
            data_rois = [data_rois meanSub.(parameters{parameter}).(modelNames{model})(~isnan(meanSub.(parameters{parameter}).(modelNames{model})(:,roi)),roi)'];
            stat.(parameters{parameter}).(modelNames{model}).N.(ROIs{roi}) = length(stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi}));
            nr_rois = [nr_rois repelem(roi, stat.(parameters{parameter}).(modelNames{model}).N.(ROIs{roi}))];
        end
        
        pnorm_jb = [];
        pnorm_sw = [];
        for roi = 1:length(ROIs)
            [stat.(parameters{parameter}).(modelNames{model}).normality_jb.(ROIs{roi}).h, stat.(parameters{parameter}).(modelNames{model}).normality_jb.(ROIs{roi}).p] = jbtest(stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi}));
            pnorm_jb = [pnorm_jb stat.(parameters{parameter}).(modelNames{model}).normality_jb.(ROIs{roi}).p];
            
            
            [stat.(parameters{parameter}).(modelNames{model}).normality_sw.(ROIs{roi}).h, stat.(parameters{parameter}).(modelNames{model}).normality_sw.(ROIs{roi}).p,stat.(parameters{parameter}).(modelNames{model}).normality_sw.(ROIs{roi}).w] = swtest(stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi}));
            pnorm_sw = [pnorm_sw stat.(parameters{parameter}).(modelNames{model}).normality_jb.(ROIs{roi}).p];
        end
        
        [~, ~, ~, adj_pnorm_jb]=fdr_bh(pnorm_jb);
        [~, ~, ~, adj_pnorm_sw]=fdr_bh(pnorm_sw);
        for roi = 1:length(ROIs)
            stat.(parameters{parameter}).(modelNames{model}).normality_jb.(ROIs{roi}).adj_pnorm = adj_pnorm_jb(roi);
            stat.(parameters{parameter}).(modelNames{model}).normality_sw.(ROIs{roi}).adj_pnorm = adj_pnorm_sw(roi);
        end
        
        % NOTE: Chose to go with only non-parametric tests, since not all parameters were normally distributed and we wanted to keep them similar
        pairw_comp = dunn(data_rois,nr_rois);
        stat.(parameters{parameter}).(modelNames{model}).dunncompare.tbl = pairw_comp;
        id_sign = find(string(pairw_comp(:,6))=="Reject H0");
        stat.(parameters{parameter}).(modelNames{model}).dunncompare.pairs = ROIs(str2double(split(pairw_comp(id_sign,1),'-')));
        stacked_plots_median = 1; % for plots
        
        %% plot altogether - stacked
        % mean & SE values per roi
        if stacked_plots_median > 0
            
            % use NaN removed data
            for roi = 1:length(ROIs)
                midpoints(roi) =median(stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi}));
                barpoints(:,roi) = bootci(1000, {@median, stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi})},'alpha',0.05);
                stat.(parameters{parameter}).(modelNames{model}).CI.(ROIs{roi})= barpoints(:,roi);
            end
            
            if string(parameters{parameter}) == "beta_ratio"
                barpoints = log10(barpoints);
                midpoints = log10(midpoints);
                stat.(parameters{parameter}).(modelNames{model}).CI.(ROIs{roi}).log = barpoints;
            end
        
            barpoints_low =  midpoints - barpoints(1,:);
            barpoints_high = barpoints(2,:) -  midpoints;
        else
            
            
            for roi = 1:length(ROIs)
                midpoints(roi) =mean(stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi}));
                std_ve = nanstd(stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi}));
                se_ve(roi) = (std_ve/sqrt(stat.(parameters{parameter}).(modelNames{model}).N.(ROIs{roi})));
            end
            barpoints_low = se_ve;
            barpoints_high = se_ve;
            
        end
        
     
        
        %% plot stacked models
        % set axis and orders
        if timing_maps ==1
            scatter_order=[4 10 9 8 2 6 5 7 1 3];
        elseif timing_maps == 0
            scatter_order=[14:16,10:13,17,1:9];
        end
        
        x_LR=1:length(ROIs);
        if string(modelNames{model}) == 'TunedLin2d'
            bar=[1,.56,.56];
        else
            bar=[.56,.56,1];
        end
        
        
        figure(counter)
        hold on
        title(parameters{parameter})
        
        if string(parameters{parameter}) == 'theta'
            plot_max = rad2deg(3/8*pi);
            plot_min = rad2deg(1/8*pi);
        elseif string(parameters{parameter}) == 'ratio'
            plot_max = 7;
            plot_min = 0;
        elseif string(parameters{parameter}) == 'beta_ratio'
            plot_max = 1;
            plot_min = -1;
        elseif string(parameters{parameter}) == 'major'
            plot_max = 1.1;
            plot_min = 0;
        else
            plot_max = 1;
            plot_min = 0;
        end
        
        axis([0,length(ROIs)+1,plot_min,plot_max])
        
        set(gcf,'units','centimeters','position',[0.1 0.1 20 12]);
        xticks(1:1:length(ROIs)+1);
        xticklabels([ROIs(scatter_order(1:end))]);
        xtickangle(90);box off
        
        if string(parameters{parameter}) == 'theta'
            yticks(plot_min:rad2deg(1/8*pi):plot_max);
            yticklabels(plot_min:rad2deg(1/8*pi):plot_max);
            
        elseif string(parameters{parameter}) == 'ratio'
            yticks(plot_min:1:plot_max);
            yticklabels(plot_min:1:plot_max);
        elseif string(parameters{parameter}) == 'beta_ratio'
            yticks(plot_min:0.5:plot_max);
            yticklabels(plot_min:0.5:plot_max);
            
        else
            yticks(plot_min:0.05:plot_max);
            yticklabels(plot_min:0.05:plot_max);
        end
        
        % create bars     
        err_bar = errorbar(x_LR,midpoints(1,scatter_order),barpoints_low(1,(scatter_order)),barpoints_high(1,(scatter_order)),'vertical', 'LineStyle', 'none', 'Color',bar,'LineWidth',3, 'Marker','o','MarkerSize',5.3,'MarkerEdgeColor','k','MarkerFaceColor',bar,'Displayname', 'TODO');
        drawnow;
        err_bar.MarkerHandle.LineWidth = 0.05;

        % redraw markers on top (otherwised covered by error bars
        scatter(x_LR,midpoints(1,scatter_order),24,'o','MarkerFaceColor',bar,'MarkerEdgeColor','k', 'HandleVisibility', 'off');
        
        
        set(gcf,'Color','w');
        if timing_maps ==0
            xlabel('Visual Field Maps','FontWeight','bold');
        elseif timing_maps == 1
            xlabel('Timing Maps','FontWeight','bold');
        end
        
        
        savename = 'model_parameters_bootciDefault_';
        if medianPerSubj == 1
            savename = strcat(savename, 'medianPerSubj');
        elseif medianPerSubj == 0
            savename = strcat(savename, 'meanPerSubj');
        elseif medianPerSubj == 2
            savename = strcat(savename, '75thPerSubj');
        elseif medianPerSubj == 3
            savename = strcat(savename, '90thPerSubj');
        end
        if sum(stacked_plots_median) > 0
            savename = strcat(savename, '_median');
        else
            savename = strcat(savename, '_mean');
        end
        savename = strcat(savename,'_',parameters{parameter},'_',modelNames{model});
        if timing_maps == 1
            savename = strcat(savename, '_timing_maps');
        elseif timing_maps == 0
            savename = strcat(savename, '_visual_field_maps');
        end
        savename = strcat(savename, '_minVE=', string(minVE));
        savename = strcat(savename, '.eps');
        export_fig(savename,'-eps','-r600','-painters');
        disp(savename)
        %close all
        hold off
        counter = counter + 1;
        
        close all
    end
end


savename = 'stat_model_params_';
if medianPerSubj == 1
    savename = strcat(savename, 'medianPerSubj');
elseif medianPerSubj == 0
    savename = strcat(savename, 'meanPerSubj');
elseif medianPerSubj == 2
    savename = strcat(savename, '75thPerSubj');
elseif medianPerSubj == 3
    savename = strcat(savename, '90thPerSubj');
end

if timing_maps == 1
    savename = strcat(savename, '_timing_maps');
elseif timing_maps == 0
    savename = strcat(savename, '_visual_field_maps');
end
savename = strcat(savename, '_minVE=', string(minVE), '.mat');
save(savename, 'stat')

%% last part: make figures with median values for various areas
monoTuned_parameter_examples_anonymous(save_path, stat, timing_maps, medianPerSubj, stacked_plots_median, minVE)
end
