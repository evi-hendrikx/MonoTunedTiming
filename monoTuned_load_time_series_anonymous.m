function time_series = monoTuned_load_time_series_anonymous(new_subjNames, paths, vfm_paths, timing_maps,save_path)
%% load time series: creates a structure called time_series, that contains coordinates, and their time series timingsweeps L1, evenscans L1 and oddscans L1 in all areas for each participant
% new_subjNames: subject names you want in the structure
% paths: where the information about each participant is stored (should be as long as and corresponding with new_subjNames)
% timing_maps: 1 for timing map ROI, 0 for visual field map ROIs
% save_path: where you want the created structure to be stored


%% general info
scan_nr_timing = 4;
scan_nr_VFM = 1;

if timing_maps == 1
    ROIs = {'TFI', 'TLS', 'TFS','TLO', 'TPCM','TPCI', 'TPCS', 'TPO','TTOA','TTOP'};
    ROILabels = [strcat('Left', ROIs), strcat('Right', ROIs)];
    ROINames = strcat(ROILabels,'Map.mat');
elseif timing_maps == 0
    ROIs =  {'IPS0', 'IPS1', 'IPS2','IPS3', 'IPS4','IPS5', 'sPCS1', 'sPCS2','iPCS','LO1','LO2', 'TO1', 'TO2','V1', 'V2','V3', 'V3AB','V4','VO1','VO2','PHC'};
    ROILabels = [strcat('left', ROIs), strcat('right', ROIs)];
    ROINames = strcat('VFMafni_', ROILabels, '.mat');
end


%% get data
for subj = 1:length(new_subjNames)
    
    % starts mrVista 3
    stimuli = [5:7]; % 5 = TimingSweeps (L1), 6 = OddScans (L1), 7 = EvenSans (L1)
    stimNames_doc = {'TimingSweeps (L1)','OddScans (L1)','EvenScans (L1)'};
        
    cd(paths{subj})
    VOLUME{1} = initHiddenGray;

    cd(vfm_paths{subj})
    VOLUME{2} = initHiddenGray;
    
    %get model info voxels
    for stim = 1:length(stimuli)
        stimName = regexprep(stimNames_doc{stim}, '[^\w'']','');
        
        % get info for specific roi
        for roi = 1:length(ROINames)
            if timing_maps == 1
                ROI_path = strcat(paths{subj},'/Gray/ROIs/',ROINames{roi});
            elseif timing_maps == 0
                switch ROINames{roi}
                    case {'VFMafni_leftLO1.mat', 'VFMafni_leftLO2.mat','VFMafni_leftVO1.mat','VFMafni_leftVO2.mat', ...
                            'VFMafni_leftV4.mat','VFMafni_leftPHC.mat','VFMafni_rightLO1.mat', 'VFMafni_rightLO2.mat', ...
                            'VFMafni_rightVO1.mat','VFMafni_rightVO2.mat','VFMafni_rightV4.mat','VFMafni_rightPHC.mat'}
                        % strrep part is a trick to get original subject numbers
                        ROI_path = strcat(paths{subj},'/ROIs/', char(strrep(extractBetween(paths{subj},'c/','C'),'_','')),'NewVFM/',ROINames{roi});
                    otherwise
                        ROI_path = strcat(paths{subj}, '/ROIs/',ROINames{roi});
                end
            end
            
            % need to check if file exists, because not all
            % subjects have all maps
            if exist(ROI_path, 'file') == 2
                load(fullfile(ROI_path))                         %Load ROIs
                [tmp, iCrds] = intersectCols(VOLUME{1}.coords, ROI.coords);         %Isolate ROI locations
                
                time_series.(new_subjNames{subj}).(stimName).(ROILabels{roi}).iCoords = iCrds;
                
                for scan = 1:scan_nr_timing
                    scan_name = ['Scan', char(num2str(scan))];
                    cd(strcat(paths{subj},'/Gray/', stimNames_doc{stim}, '/TSeries/', scan_name));
                    load('tSeries1');
                    time_series.(new_subjNames{subj}).(stimName).(ROILabels{roi}).(scan_name) = tSeries(:,iCrds);
                end
                close all
                
                %Isolate ROI locations
                if stim == length(stimuli) % only needs to be done once per roi
                    [tmp, iCrds_vfm] = intersectCols(VOLUME{2}.coords, ROI.coords);
                    
                    time_series.(new_subjNames{subj}).VFM.(ROILabels{roi}).iCoords = iCrds;
                    time_series.(new_subjNames{subj}).VFM.(ROILabels{roi}).iCoords_vfm = iCrds_vfm;
                    for scan = 1:scan_nr_VFM
                        scan_name = ['Scan', char(num2str(scan))];
                        cd(strcat(vfm_paths{subj},'/Gray/VFM (L1)/TSeries/', scan_name));
                        load('tSeries1');
                        time_series.(new_subjNames{subj}).VFM.(ROILabels{roi}).(scan_name) = tSeries(:,iCrds_vfm);
                    end
                    close all
                end

                
            else
                fprintf('Map does not exist: %s for %s \n',ROILabels{roi}, new_subjNames{subj})
            end
        end
    end
    mrvCleanWorkspace
end

cd(save_path);
savename = 'time_series_';
if timing_maps == 1
    savename = strcat(savename, 'timing');
elseif timing_maps == 0
    savename = strcat(savename, 'visual_field');
end
savename = strcat(savename, '_maps.mat');
save(savename, 'time_series');



