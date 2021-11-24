function ve_data = monoTuned_load_timing_data_anonymous_incl_vfminfo(new_subjNames, paths, timing_maps,save_path)
%% load timing data: creates a structure called ve_data, that contains coordinates, variance explained, rss, rawrss,and model parameters for timingsweeps L1, evenscans L1 and oddscans L1 for both models in all areas for each participant
% new_subjNames: subject names you want in the structure
% paths: where the information about each participant is stored (should be as long as and corresponding with new_subjNames)
% timing_maps: 1 for timing map ROI, 0 for visual field map ROIs
% save_path: where you want the created structure to be stored

%% general info
modelFieldNames = {'MonoOcc','TunedLin2d','MonoLinDurCompFreq','MonoLinDurLinFreq'};
modelNames = {'*202109*Lin-compressiveXcompressiveYNoNormOccupancy-DurFreq-20-gFit*.mat',...
    '*Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit-fFit.mat',...
    '*202109*linearXcompressiveY*.mat','*202109*linearXlinearY*.mat'};
cv_or_not = {'crossValidated','used_for_crossValidation'};
file_VF = 'rmImported_copy_retModel-Decimate-Freexy-FreeHRF_NoDC_OneG-fFit-fFit-fFit';%'combined';


%% ROI names
if timing_maps == 1
    ROIs = {'TFI', 'TLS', 'TFS','TLO', 'TPCM','TPCI', 'TPCS', 'TPO','TTOA','TTOP'};
    ROILabels = [strcat('Left', ROIs), strcat('Right', ROIs)];
    ROINames = strcat(ROILabels,'Map.mat');
elseif timing_maps == 0
    ROIs =  {'IPS0', 'IPS1', 'IPS2','IPS3', 'IPS4','IPS5', 'sPCS1', 'sPCS2','iPCS','LO1','LO2', 'TO1', 'TO2','V1', 'V2','V3', 'V3AB'};
    ROILabels = [strcat('left', ROIs), strcat('right', ROIs)];
    ROINames = strcat('VFMafni_', ROILabels, '.mat');
end


%% get data
for subj = 1:length(new_subjNames)
    %get voxels' eccentricities
    eccentricity_path = strcat(paths{subj},'/Gray/TimingSweeps (L1)/', file_VF);
    load(fullfile(eccentricity_path))
    eccentricities = rmGet(model{1}, 'eccentricity');
    vfm_sigma = rmGet(model{1}, 'sigma'); % sigma major and minor are the same, rmGet takes the average value of them (and so theta is always 0)
    vfm_ve = rmGet(model{1}, 've');
    vfm_xs = model{1}.x0;
    vfm_ys = model{1}.y0;
    
    % starts mrVista 3
    cd(paths{subj})
    VOLUME{1} = initHiddenGray;
    stimuli = [5:7]; % 5 = TimingSweeps (L1), 6 = OddScans (L1), 7 = EvenSans (L1)
    stimNames_doc = {'TimingSweeps (L1)','OddScans (L1)','EvenScans (L1)'};
    
    for cv = 1:length(cv_or_not)
        
        %get model info voxels
        for stim = 1:length(stimuli)
            %Get correct folder for stimulus configuration
            stimName = regexprep(stimNames_doc{stim}, '[^\w'']','');
            
            for models = 1:length(modelFieldNames)
                if strcmp(modelFieldNames{models},'TunedLin2d')~=1  % monotonic model(s)
                    stimFolder = fullfile(paths{subj},'/Gray',stimNames_doc(stim), 'Monotonic');

                    if stim>1 && cv == 1 %for odd/even scans % use cross-validated data
                        stimFolder = fullfile(paths{subj}, '/Gray',stimNames_doc(stim),'Monotonic/xvalRefit');
                    elseif stim>1 && cv == 2 %for odd/even scans % use uncross-validated data
                        stimFolder = fullfile(paths{subj}, '/Gray',stimNames_doc(stim),'Monotonic/xval');
                    end
                    
                else % tuned model
                    stimFolder = fullfile(paths{subj},'/Gray',stimNames_doc(stim),'SearchFitFreeExponent');

                    if stim>1 && cv == 1 %for odd/even scans % use cross-validated data
                        stimFolder = fullfile(paths{subj},'/Gray',stimNames_doc(stim),'SearchFitFreeExponent/xvalRefit');
                    elseif stim>1 && cv == 2 %for odd/even scans % use uncross-validated data
                        stimFolder = fullfile(paths{subj},'/Gray',stimNames_doc(stim),'SearchFitFreeExponent/xval');
                    end
                end
                
                cd(char(stimFolder));
                load(strcat(ls(modelNames{models})));
                
                % get info for specific roi
                for roi = 1:length(ROINames)
                    cd(paths{subj})
                    if timing_maps == 1
                        ROI_path = strcat(paths{subj},'/Gray/ROIs/',ROINames{roi});
                    elseif timing_maps == 0
                        ROI_path = strcat(paths{subj}, '/ROIs/',ROINames{roi});
                    end
                    
                    % need to check if file exists, because not all
                    % subjects have all maps
                    if exist(ROI_path, 'file') == 2
                        load(fullfile(ROI_path))                         %Load ROIs
                        [tmp, iCrds] = intersectCols(VOLUME{1}.coords, ROI.coords);         %Isolate ROI locations
                        varianceExplained = rmGet(model{1},'ve');                          %Get variance explained (VE) Calculated with rmGet
                        varianceExplained = varianceExplained(iCrds);                      %Get VE in ROI
                        xs = model{1}.x0(iCrds);
                        ys = model{1}.y0(iCrds);
                        rss = model{1}.rss(iCrds);
                        rawrss = model{1}.rawrss(iCrds);
                        
                        grayROI=load(fullfile(strcat(paths{subj},'/Gray/ROIs/gray-Layer1.mat')));
                        grayCoords=grayROI.ROI.coords;
                        [tmp, betaCrds,beta_iCrds] = intersectCols(grayCoords, VOLUME{1}.coords(:,iCrds));
                        betas = squeeze(model{1}.beta(1,betaCrds,:));
                        if isfield(params, 'seperateRunBetas') && params.seperateRunBetas
                            beta1=betas(:,[1 2 3 4]);
                            beta2=betas(:,[5 6 7 8]);
                        else
                            beta1=betas(:,1);
                            beta2=betas(:,2);
                        end
                        
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).iCoords = iCrds;
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).xs = xs; % best duration (tuned) or compressive exponent on duration (mono)
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).varianceExplained = varianceExplained; %Save VE
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).ys = ys; % best period (tuned) or compressive exponent on period (mono)
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).rss = rss;
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).rawrss = rawrss;
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).beta1 = beta1;
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).beta2 = beta2;
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).beta_Coords = iCrds(beta_iCrds);
                        
                        
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).ecc = eccentricities(iCrds);
                       
                        if strcmp(modelFieldNames{models},'TunedLin2d')==1
                            ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).sigma_major = model{1}.sigma.major(iCrds);
                            ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).sigma_minor = model{1}.sigma.minor(iCrds);
                            ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).theta = model{1}.sigma.theta(iCrds);
                            ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).exponent = model{1}.exp(iCrds);
                        end
                        
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).VFM_sigma = vfm_sigma(iCrds);
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).VFM_ve = vfm_ve(iCrds);
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).VFM_xs = vfm_xs(iCrds);
                        ve_data.(modelFieldNames{models}).(new_subjNames{subj}).(stimName).(cv_or_not{cv}).(ROILabels{roi}).VFM_ys = vfm_ys(iCrds);
                        
                        
                    else
                        fprintf('Map does not exist: %s for %s \n',ROINames{roi}, new_subjNames{subj})
                    end
                end
            end
        end
    end
    mrvCleanWorkspace
end

cd(save_path);
savename = 've_data_';
if timing_maps == 1
    savename = strcat(savename, 'timing');
elseif timing_maps == 0
    savename = strcat(savename, 'visual_field');
end
savename = strcat(savename, '_maps.mat');
save(savename, 've_data');
end



