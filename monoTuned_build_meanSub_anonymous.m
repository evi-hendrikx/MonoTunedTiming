function [meanSub, ANOVAstruc] = monoTuned_build_meanSub_anonymous(cv_data, minVE, ecc_analysis, near_ecc, far_ecc, max_ecc,medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod)
%% build meanSub: builds structures (meanSub and ANOVAstruc) on which the analyses will be done on
% creates a structure called meanSub, that contains the summary statistic of the variance explained (or difference in VE) 
% of voxels with a variance explained above a threshold (minVE) for at least one of the models in cv_data
% summary statistics are created on the level of hemisphere and cross-validation run. i.e., 8 participants x 2 hemispheres x 2 runs = 32 datapoints per ROI
% 
% also creates ANOVAstruc that contains the same data and will be used for the ANOVAs (needs sligthly different order/ removal of NaN values)
%
%   cv_data: structure containing coordinates, variance explained, rss, rawrss, and model parameters for evenscans L1 and oddscans L1 for both models in all areas for each participant
%   minVE: threshold for inclusion of the voxels
%   ecc_analysis: 1 if following analyses are on eccentricity ranges, 0 if not (i.e., in model comparisons)
   %    near_ecc: threshold, voxel with an eccentricity below this value are included
   %    far_ecc: threshold, voxel with an eccentricity above this value are included
   %    max_ecc: threshold, voxel with an eccentricity above this value are excluded
%   medianPerSubj: determines the summary statistic: 1 = means; 2 = medians; 3 = 75th percentile; 4 = 90th percentile
%   minDuration,maxDuration,minPeriod,maxPeriod: ranges outside which the tuned model is considered invalid

%% get variable names
modelFieldNames = fieldnames(cv_data);
subjNames = fieldnames(cv_data.(modelFieldNames{1}));
condNames = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{1})); % even & odd
ROILabels = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{1}).(condNames{1}).crossValidated); %possible, because S3 has all maps
ROIs = unique(erase(ROILabels, ["Right","Left","right","left"]),'stable');
hemispheres = {'left','right'};

%% prepare data for ANOVA

for subj = 1:length(subjNames)
    for roi = 1:length(ROILabels)  
        
        if roi <= length(ROIs) % left
            roi_id = roi;
            hemisphere = 1;
        else % right
            roi_id = roi - length(ROIs); 
            hemisphere = 2;
        end
        
        roi_names_subj = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{subj}).(condNames{1}).crossValidated); % because not all subjects have all ROIs
        for condition = 1:length(condNames)
            
            if any(strcmp(roi_names_subj, ROILabels{roi}))

                %NOTE: if tuned model is involved, it's always the second model (due to build cv data)
                ve_mod1_Current = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).varianceExplained';
                ve_mod2_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).varianceExplained';
                
                % necessary for voxel selection
                ve_mod1_Current_pre_cv = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).used_for_crossValidation.(ROILabels{roi}).varianceExplained';
                ve_mod2_Current_pre_cv = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).used_for_crossValidation.(ROILabels{roi}).varianceExplained';
                
                % if parameters out of range for the tuned model: set VE to 0 (these are later excluded because <minVE)
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    xs_tuned_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).xs';
                    ys_tuned_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).ys';

                    outsideRange = find(sum([xs_tuned_Current<=minDuration, xs_tuned_Current>=maxDuration,ys_tuned_Current<=minPeriod,ys_tuned_Current>=maxPeriod],2)>0);
                    ve_mod2_Current(outsideRange)=0;
                    ve_mod2_Current_pre_cv(outsideRange)=0;
                end
                 
                % calculate difference between the models for each voxel
                % NOTE: depends on order modelFieldNames EXCEPT if there is a tuned model: always tuned - mono
                ve_diff_Current = ve_mod2_Current - ve_mod1_Current;

                % voxel selection, remove if < minVE for both models in uncrossvalidated data
                % It is allowed to be below treshold for one of them
                removeVE = find(sum([ve_mod1_Current_pre_cv<=minVE,ve_mod2_Current_pre_cv<= minVE],2)>1); 
                ve_mod1_Current(removeVE)=[];
                ve_mod2_Current(removeVE)=[];
                ve_diff_Current(removeVE)=[];
                
                
                % if eccentricity data, remove irrelevant eccentricities 
                if ecc_analysis == 1
                    eccClose = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{1})).crossValidated.(ROILabels{roi}).ecc';  %condition and model don't matter: the same ecc everywhere
                    eccClose(removeVE) = [];
                    
                    if near_ecc ~= far_ecc
                        mid_eccs = find(sum([eccClose >= near_ecc, eccClose <= far_ecc],2)>1);% I want a voxel removed if it falls between the low and high eccentricity ranges
                        ve_mod1_Current(mid_eccs) = [];
                        ve_mod2_Current(mid_eccs) = [];
                        ve_diff_Current(mid_eccs) = [];
                        eccClose(mid_eccs)=[];
                    end
                    
                    ve_mod1_Current(eccClose>=max_ecc) = []; % I want a voxel removed if it falls above max ecc
                    ve_mod2_Current(eccClose>=max_ecc)=[];
                    ve_diff_Current(eccClose>=max_ecc)=[];
                    eccClose(eccClose>=max_ecc) = [];
                    
                    % split the data in near and far eccentricities
                    ve_mod1_Current_near = ve_mod1_Current(eccClose<near_ecc);
                    ve_mod1_Current_far = ve_mod1_Current(eccClose>far_ecc);
                    ve_mod2_Current_near = ve_mod2_Current(eccClose<near_ecc);
                    ve_mod2_Current_far = ve_mod2_Current(eccClose>far_ecc);
                    ve_diff_Current_near = ve_diff_Current(eccClose<near_ecc);
                    ve_diff_Current_far = ve_diff_Current(eccClose>far_ecc);
                end
                
            % make empty if ROI doesn't exist (otherwise filled with 0 instead of NaN)
            else
                ve_mod1_Current = [];
                ve_mod2_Current = [];
                ve_diff_Current = [];
                if ecc_analysis == 1
                    ve_mod1_Current_near = [];
                    ve_mod1_Current_far = [];
                    ve_mod2_Current_near = [];
                    ve_mod2_Current_far = [];
                    ve_diff_Current_near = [];
                    ve_diff_Current_far = [];
                end
            end
            
            if ecc_analysis == 1
                if medianPerSubj == 0
                    meanSub.(modelFieldNames{1}).near(subj,hemisphere,condition,roi_id)=nanmean(ve_mod1_Current_near);
                    meanSub.(modelFieldNames{1}).far(subj,hemisphere,condition,roi_id)=nanmean(ve_mod1_Current_far);
                    meanSub.(modelFieldNames{2}).near(subj,hemisphere,condition,roi_id)=nanmean(ve_mod2_Current_near);
                    meanSub.(modelFieldNames{2}).far(subj,hemisphere,condition,roi_id)=nanmean(ve_mod2_Current_far);
                    meanSub.diff.near(subj,hemisphere,condition,roi_id)=nanmean(ve_diff_Current_near);
                    meanSub.diff.far(subj,hemisphere,condition,roi_id)=nanmean(ve_diff_Current_far);
                elseif medianPerSubj == 1
                    meanSub.(modelFieldNames{1}).near(subj,hemisphere,condition,roi_id)=nanmedian(ve_mod1_Current_near);
                    meanSub.(modelFieldNames{1}).far(subj,hemisphere,condition,roi_id)=nanmedian(ve_mod1_Current_far);
                    meanSub.(modelFieldNames{2}).near(subj,hemisphere,condition,roi_id)=nanmedian(ve_mod2_Current_near);
                    meanSub.(modelFieldNames{2}).far(subj,hemisphere,condition,roi_id)=nanmedian(ve_mod2_Current_far);
                    meanSub.diff.near(subj,hemisphere,condition,roi_id)=nanmedian(ve_diff_Current_near);
                    meanSub.diff.far(subj,hemisphere,condition,roi_id)=nanmedian(ve_diff_Current_far);
                elseif medianPerSubj == 2
                    meanSub.(modelFieldNames{1}).near(subj,hemisphere,condition,roi_id)=prctile(ve_mod1_Current_near,75);
                    meanSub.(modelFieldNames{1}).far(subj,hemisphere,condition,roi_id)=prctile(ve_mod1_Current_far,75);
                    meanSub.(modelFieldNames{2}).near(subj,hemisphere,condition,roi_id)=prctile(ve_mod2_Current_near,75);
                    meanSub.(modelFieldNames{2}).far(subj,hemisphere,condition,roi_id)=prctile(ve_mod2_Current_far,75);
                    meanSub.diff.near(subj,hemisphere,condition,roi_id)=prctile(ve_diff_Current_near,75);
                    meanSub.diff.far(subj,hemisphere,condition,roi_id)=prctile(ve_diff_Current_far,75);
                elseif medianPerSubj == 3
                    meanSub.(modelFieldNames{1}).near(subj,hemisphere,condition,roi_id)=prctile(ve_mod1_Current_near,90);
                    meanSub.(modelFieldNames{1}).far(subj,hemisphere,condition,roi_id)=prctile(ve_mod1_Current_far,90);
                    meanSub.(modelFieldNames{2}).near(subj,hemisphere,condition,roi_id)=prctile(ve_mod2_Current_near,90);
                    meanSub.(modelFieldNames{2}).far(subj,hemisphere,condition,roi_id)=prctile(ve_mod2_Current_far,90);
                    meanSub.diff.near(subj,hemisphere,condition,roi_id)=prctile(ve_diff_Current_near,90);
                    meanSub.diff.far(subj,hemisphere,condition,roi_id)=prctile(ve_diff_Current_far,90);
                end
                
                
            else
                if medianPerSubj == 0
                    meanSub.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmean(ve_mod1_Current);
                    meanSub.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmean(ve_mod2_Current);
                    meanSub.diff(subj,hemisphere,condition,roi_id)=nanmean(ve_diff_Current);
                    
                elseif medianPerSubj == 1
                    meanSub.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=nanmedian(ve_mod1_Current);
                    meanSub.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=nanmedian(ve_mod2_Current);
                    meanSub.diff(subj,hemisphere,condition,roi_id)=nanmedian(ve_diff_Current);
                    
                elseif medianPerSubj == 2
                    meanSub.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(ve_mod1_Current,75);
                    meanSub.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(ve_mod2_Current,75);
                    meanSub.diff(subj,hemisphere,condition,roi_id)=prctile(ve_diff_Current,75);
                    
                elseif medianPerSubj == 3
                    meanSub.(modelFieldNames{1})(subj,hemisphere,condition,roi_id)=prctile(ve_mod1_Current,90);
                    meanSub.(modelFieldNames{2})(subj,hemisphere,condition,roi_id)=prctile(ve_mod2_Current,90);
                    meanSub.diff(subj,hemisphere,condition,roi_id)=prctile(ve_diff_Current,90);
                end
                
            end
            
            meanSubSubj(subj,hemisphere,condition,roi_id) = subjNames(subj);
            meanSubRoi(subj,hemisphere,condition,roi_id) = ROIs(roi_id);
            meanSubHemi(subj,hemisphere,condition,roi_id) = hemispheres(hemisphere);
            meanSubCond(subj,hemisphere,condition,roi_id) = condNames(condition);
                   
        end
    end
end

ANOVAstruc.roi = [meanSubRoi(:);meanSubRoi(:)];
ANOVAstruc.subj = [meanSubSubj(:);meanSubSubj(:)];

if ecc_analysis ==  0 
    ANOVAstruc.VE = [meanSub.(modelFieldNames{1})(:);meanSub.(modelFieldNames{2})(:)];
    ANOVAstruc.model = repelem(modelFieldNames,length(meanSub.(modelFieldNames{1})(:)));
    
    % since we're working with a threshold on variance explained,
    % sometimes average value doesn't exist.
    ANOVAstruc.roi(isnan(ANOVAstruc.VE))=[];
    ANOVAstruc.subj(isnan(ANOVAstruc.VE))=[];
    ANOVAstruc.model(isnan(ANOVAstruc.VE))=[];
    ANOVAstruc.VE(isnan(ANOVAstruc.VE))=[];

elseif ecc_analysis == 1 
  
    ANOVAstruc.VE.(modelFieldNames{1}) = [meanSub.(modelFieldNames{1}).near(:);meanSub.(modelFieldNames{1}).far(:)];
    ANOVAstruc.VE.(modelFieldNames{2}) = [meanSub.(modelFieldNames{2}).near(:);meanSub.(modelFieldNames{2}).far(:)];
    ANOVAstruc.VE.diff = [meanSub.diff.near(:);meanSub.diff.far(:)];
    ANOVAstruc.ecc = repelem({'near','far'}, length(meanSub.(modelFieldNames{1}).near(:)));
    
    % When an entire eccentricity range does not have VE values >0 for any of the models
    % same voxels for both models, because the voxel selection takes place on their combination. So
    % model 2 can be used everywhere
    ANOVAstruc.roi(isnan(ANOVAstruc.VE.(modelFieldNames{2})))=[];
    ANOVAstruc.subj(isnan(ANOVAstruc.VE.(modelFieldNames{2})))=[];
    ANOVAstruc.ecc(isnan(ANOVAstruc.VE.(modelFieldNames{2})))=[];
    ANOVAstruc.VE.(modelFieldNames{1})(isnan(ANOVAstruc.VE.(modelFieldNames{2})))=[];
    ANOVAstruc.VE.diff(isnan(ANOVAstruc.VE.(modelFieldNames{2})))=[];
    ANOVAstruc.VE.(modelFieldNames{2})(isnan(ANOVAstruc.VE.(modelFieldNames{2})))=[];
end


%reshape meansub because later scripts expect all data points per roi under
%each other
if ecc_analysis == 1
   meanSub.(modelFieldNames{1}).near = reshape(meanSub.(modelFieldNames{1}).near,length(subjNames)*length(condNames)*2,length(ROIs));
   meanSub.(modelFieldNames{1}).far = reshape(meanSub.(modelFieldNames{1}).far,length(subjNames)*length(condNames)*2,length(ROIs));
   meanSub.(modelFieldNames{2}).near = reshape(meanSub.(modelFieldNames{2}).near,length(subjNames)*length(condNames)*2,length(ROIs));
   meanSub.(modelFieldNames{2}).far = reshape(meanSub.(modelFieldNames{2}).far,length(subjNames)*length(condNames)*2,length(ROIs));
   meanSub.diff.near = reshape(meanSub.diff.near,length(subjNames)*length(condNames)*2,length(ROIs));
   meanSub.diff.far = reshape(meanSub.diff.far,length(subjNames)*length(condNames)*2,length(ROIs));
elseif ecc_analysis == 0
    meanSub.(modelFieldNames{1}) = reshape(meanSub.(modelFieldNames{1}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.(modelFieldNames{2}) = reshape(meanSub.(modelFieldNames{2}),length(subjNames)*length(condNames)*2,length(ROIs));
    meanSub.diff = reshape(meanSub.diff,length(subjNames)*length(condNames)*2,length(ROIs));
end

  % these provide an extra check to see which datapoint belongs to which
  % participant in which roi, hemi and condition
  meanSub.subj = reshape(meanSubSubj,length(subjNames)*length(condNames)*2,length(ROIs));
  meanSub.roi = reshape(meanSubRoi,length(subjNames)*length(condNames)*2,length(ROIs));
  meanSub.hemi = reshape(meanSubHemi,length(subjNames)*length(condNames)*2,length(ROIs));
  meanSub.cond = reshape(meanSubCond,length(subjNames)*length(condNames)*2,length(ROIs));
end