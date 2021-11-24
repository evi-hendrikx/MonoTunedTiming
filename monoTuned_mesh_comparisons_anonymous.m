function monoTuned_mesh_comparisons_anonymous(subjNames, paths, save_path, minVE,minDuration,maxDuration,minPeriod,maxPeriod)
%% mesh comparisons: creates meshes where the VE of the best-performing (monotonic or tuned) model is portrayed for each voxel
% gets x0, y0, rss, and rawrss for each voxel for the tuned and monotonic
% model
% calculates average variance explained (VE) for each model in each voxel over cross-validated runs 
% (NOTE: if the x0 (preferred duration) or y0 (preferred frequency) of BOTH runs in the tuned model are outside of the range, this voxel's VE is 
% set to 0. During model comparisons this is done per run, since both runs can there be included as separate datapoints.
% So the results slightly seem different in some of the areas, but they are the same data)
%
% makes mesh with voxels where monotonic fits best
% makes mesh with voxels where tuned fits best
% overlays these meshes
%
% subjNames: labels of the participants
% paths: where participant's files are stored
% save_path: general folder where you want results to be stored
% minVE: threshold for inclusion of the voxels
% minDuration,maxDuration,minPeriod,maxPeriod: ranges outside which the tuned model is considered invalid

% Note to self: models are the other way around compared to the other scripts
modelFieldNames = {'TunedLin2d', 'MonoOcc'};
save_mesh = strcat(save_path, '/mesh_comparisons');


% load in data
for subj = 1:length(subjNames)
    cd(paths{subj})
    
    VOLUME{1} = initHiddenGray;
    stimuli = [5:7]; % 5 = TimingSweeps (L1), 6 = OddScans (L1), 7 = EvenSans (L1)
    
    for stim = 1:length(stimuli)
        VOLUME{1}.curDataType = stimuli(stim);
        
        for models = 1:length(modelFieldNames)
            cd(paths{subj})
            modelNames = {'*Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit-fFit.mat', '*202109*Lin-compressiveXcompressiveYNoNormOccupancy-DurFreq-20-gFit*.mat'};
            
            %Get correct folder for stimulus configuration
            if stim == 1 %timingsweeps
                if models==1 %tuned
                    stimFolder = fullfile('Gray',dataTYPES(1,stimuli(stim)).name,'SearchFitFreeExponent');
                elseif models==2 % monotonic
                    stimFolder = fullfile('Gray',dataTYPES(1,stimuli(stim)).name,'Monotonic');
                end
                
                cd(stimFolder);
                load(strcat(ls(modelNames{models})));
                
                if models == 1 %tuned
                    tuneAll = rmGet(model{1},'ve');
                    tuneAll_xs = model{1}.x0;
                    tuneAll_ys = model{1}.y0;
                    tuneAll_rss = model{1}.rss;
                    tuneAll_rawrss = model{1}.rawrss;
                end
                
                if models == 2 %monotonic
                    monoAll = rmGet(model{1},'ve');
                    monoAll_xs = model{1}.x0;
                    monoAll_ys = model{1}.y0;
                    monoAll_rss = model{1}.rss;
                    monoAll_rawrss = model{1}.rawrss;
                end
                
            elseif stim == 2 || stim == 3
                
                if models == 1 %tuned
                    stimFolder = fullfile('Gray',dataTYPES(1,stimuli(stim)).name,'SearchFitFreeExponent/xvalRefit');
                    
                    cd(stimFolder);
                    load(strcat(ls(modelNames{models})));
                    if stim == 2 %odd
                        tuneOdd = rmGet(model{1},'ve');
                        tuneOdd_xs = model{1}.x0;
                        tuneOdd_ys = model{1}.y0;
                        tuneOdd_rss = model{1}.rss;
                        tuneOdd_rawrss = model{1}.rawrss;
                    elseif stim == 3 %even
                        tuneEven = rmGet(model{1},'ve');
                        tuneEven_xs = model{1}.x0;
                        tuneEven_ys = model{1}.y0;
                        tuneEven_rss = model{1}.rss;
                        tuneEven_rawrss = model{1}.rawrss;
                    end
                elseif models == 2 % monotonic
                    stimFolder = fullfile('Gray',dataTYPES(1,stimuli(stim)).name,'Monotonic/xvalRefit');
                    
                    cd(stimFolder);
                    load(strcat(ls(modelNames{models})));
                    if stim == 2 % odd
                        monoOdd = rmGet(model{1},'ve');
                        monoOdd_xs = model{1}.x0;
                        monoOdd_ys = model{1}.y0;
                        monoOdd_rss = model{1}.rss;
                        monoOdd_rawrss = model{1}.rawrss;
                    elseif stim == 3 %even
                        monoEven = rmGet(model{1},'ve');
                        monoEven_xs = model{1}.x0;
                        monoEven_ys = model{1}.y0;
                        monoEven_rss = model{1}.rss;
                        monoEven_rawrss = model{1}.rawrss;
                    end
                end
            end
            
        end
        
    end
    
    % calculate average VEs (meshes can only portray 1 value per voxel)
    cd(paths{subj})
    
    mono_rss = monoOdd_rss + monoEven_rss;
    mono_rawrss = monoOdd_rawrss + monoEven_rawrss;
    
    mono = 1 - mono_rss ./ mono_rawrss;
    mono(~isfinite(mono)) = 0;
    mono = max(mono, 0);
    mono = min(mono, 1);
    
    tune_rss = tuneOdd_rss + tuneEven_rss;
    tune_rawrss = tuneOdd_rawrss + tuneEven_rawrss;
    
    tune = 1 - tune_rss ./ tune_rawrss;
    tune(~isfinite(tune)) = 0;
    tune = max(tune, 0);
    tune = min(tune, 1);
    
    tune(tuneOdd_xs<=minDuration & tuneEven_xs<=minDuration) = 0;
    tune(tuneOdd_xs>=maxDuration & tuneEven_xs>=maxDuration) = 0;
    tune(tuneOdd_ys<=minPeriod & tuneEven_ys<=minPeriod) = 0;
    tune(tuneOdd_ys>=maxPeriod & tuneEven_ys>=maxPeriod) = 0;
    
    % prepare mesh settings
    meshNames = {'meshBothInflated','meshLeftInflated','meshRightInflated'};
    meshSettings = load('MeshSettings.mat');
    subjSettings = length(meshSettings.settings);
    currentSettings = [subjSettings-2,subjSettings-1,subjSettings]; %back,left,right
    
    close all
    mrvCleanWorkspace
    
    saveNames = {'monoTuned_back(cv)','monoTuned_left(cv)','monoTuned_right(cv)'};
    saveNames = strcat(saveNames, '_minVE=', string(minVE));
    
    for meshes=1:3 % back left right
        
        mrVista 3
        H = open3DWindow;
        
        % load mesh
        meshPath = strcat([pwd,'/',meshNames{meshes}]);
        [VOLUME{1},OK] = meshLoad(VOLUME{1},meshPath,1);
        msh = viewGet(VOLUME{1}, 'currentmesh'); % refresh happens there
        if ~isempty(msh)
            VOLUME{1} = viewSet(VOLUME{1},'currentmesh',msh);
        end
        mrmSet(msh, 'cursoroff');
        
        MSH = viewGet(VOLUME{1}, 'Mesh');
        vertexGrayMap = mrmMapVerticesToGray( meshGet(MSH, 'initialvertices'), viewGet(VOLUME{1}, 'nodes'),...
            viewGet(VOLUME{1}, 'mmPerVox'), viewGet(VOLUME{1}, 'edges') );
        MSH = meshSet(MSH, 'vertexgraymap', vertexGrayMap); VOLUME{1} = viewSet(VOLUME{1}, 'Mesh', MSH);
        clear MSH vertexGrayMap
        meshRetrieveSettings(viewGet(VOLUME{1}, 'CurMesh'), currentSettings(meshes));
        
        
        %% check values per voxel
        model_tuned.x0=[];model_tuned.rss=[];model_tuned.rawrss=[];
        model_mono.x0=[];model_mono.rss=[];model_mono.rawrss=[];
        
        for idx=1:size(mono,2)
            if (tune(idx) > mono(idx) && (tune(idx) > minVE)) %tuned response fits better than monotonic
                model_tuned.x0(idx) = tune(idx);
                model_tuned.rss(idx) = tune_rss(idx);
                model_tuned.rawrss(idx) = tune_rawrss(idx);
                
                model_mono.x0(idx) = -1;
                model_mono.rss(idx) = -1;
                model_mono.rawrss(idx) = -1;
                
            elseif (tune(idx) < mono(idx) && (mono(idx) > minVE)) %monotonic response fits better than tuned
                model_tuned.x0(idx) = -1;
                model_tuned.rss(idx) = -1;
                model_tuned.rawrss(idx) = -1;
                
                model_mono.x0(idx) = mono(idx);
                model_mono.rss(idx) = mono_rss(idx);
                model_mono.rawrss(idx) = mono_rawrss(idx);
                
            elseif (tune(idx) <= minVE && mono(idx) <= minVE) %neither model fits well 
                model_tuned.x0(idx) = -1;
                model_tuned.rss(idx) = -1;
                model_tuned.rawrss(idx) = -1;
                
                model_mono.x0(idx) = -1;
                model_mono.rss(idx) = -1;
                model_mono.rawrss(idx) = -1;
                
            elseif (tune(idx) == mono(idx) && (tune(idx) > minVE)) %tuned and monotonic response fits equally well, coin-flip
                monotuneFlip = randi(2);
                if monotuneFlip == 1 %choose mono
                    model_tuned.x0(idx) = -1;
                    model_tuned.rss(idx) = -1;
                    model_tuned.rawrss(idx) = -1;
                    
                    model_mono.x0(idx) = mono(idx);
                    model_mono.rss(idx) = mono_rss(idx);
                    model_mono.rawrss(idx) = mono_rawrss(idx);
                    
                elseif monotuneFlip == 2 %choose tuned
                    model_tuned.x0(idx) = tune(idx);
                    model_tuned.rss(idx) = tune_rss(idx);
                    model_tuned.rawrss(idx) = tune_rawrss(idx);
                    
                    model_mono.x0(idx) = -1;
                    model_mono.rss(idx) = -1;
                    model_mono.rawrss(idx) = -1;
                end
                
            end
        end
        
        %% save tuned best
        model{1}.x0 = model_tuned.x0;
        model{1}.rss = model_tuned.rss;
        model{1}.rawrss = model_tuned.rawrss;
        
        save_name = 'monotonicTuned_tuned(cv)_';
        cd(save_mesh),save(strcat(save_name,'voxel_selection=2_minVE=', string(minVE), '_', subjNames{subj},'.mat'),'model','params');
        cd(paths{subj})
        
        rmName = char(strcat(save_mesh,'/', save_name, 'voxel_selection=2_minVE=', string(minVE), '_', subjNames{subj},'.mat'));
        VOLUME{1} = rmSelect(VOLUME{1}, 1, rmName); VOLUME{1} = rmLoadDefault(VOLUME{1});
        VOLUME{1} = setDisplayMode(VOLUME{1},'co'); VOLUME{1} = refreshScreen(VOLUME{1});
        modeInfo = VOLUME{1}.ui.(['co' 'Mode']);
        VOLUME{1} = setClipMode(VOLUME{1}, 'co', [0 1]);
        
        VOLUME{1}.ui.mapMode.numGrays=128;
        VOLUME{1}.ui.mapMode.numColors=128;
        VOLUME{1}.ui.mapMode.clipMode = [0 1];
        cmap_tuned=zeros(256,3);
        cmap_tuned(1:128,:)=[linspace(0,1,128)',linspace(0,1,128)',linspace(0,1,128)'];
        
        load vik.mat
        vik128 = vik(round(linspace(129,256,128)),:);
        cmap_tuned(129:256,:) = vik128;
        
        modeInfo.cmap = cmap_tuned;
        modeInfo.name = 'cmap_tuned';
        VOLUME{1}.ui.(['co' 'Mode']) = modeInfo;
        VOLUME{1} = setCothresh(VOLUME{1}, minVE);VOLUME{1} = refreshScreen(VOLUME{1}, 1);
        VOLUME{1} = setClipMode(VOLUME{1}, 'co', [0 1]);
        [~, ~, ~, tunedColors] = meshColorOverlay(VOLUME{1},0);
        VOLUME{1} = meshColorOverlay(VOLUME{1});
        
        %% save mono best
        model{1}.x0 = model_mono.x0;
        model{1}.rss = model_mono.rss;
        model{1}.rawrss = model_mono.rawrss;
        
        save_name = 'monotonicTuned_mono(cv)_';
        cd(save_mesh),save(strcat(save_name, 'voxel_selection=2_minVE=', string(minVE), '_', subjNames{subj},'.mat'),'model','params');
        cd(paths{subj})
        
        rmName = char(strcat(save_mesh,'/', save_name,'voxel_selection=2_minVE=', string(minVE),'_', subjNames{subj}, '.mat'));
        VOLUME{1} = rmSelect(VOLUME{1}, 1, rmName);
        VOLUME{1} = rmLoadDefault(VOLUME{1});
        VOLUME{1} = setDisplayMode(VOLUME{1},'co'); VOLUME{1} = refreshScreen(VOLUME{1});
        modeInfo = VOLUME{1}.ui.(['co' 'Mode']);
        VOLUME{1} = setClipMode(VOLUME{1}, 'co', [0 1]);
        
        VOLUME{1}.ui.mapMode.numGrays=128;
        VOLUME{1}.ui.mapMode.numColors=128;
        VOLUME{1}.ui.mapMode.clipMode = [0 1];
        cmap_monoIncrease=zeros(256,3);
        cmap_monoincrease(1:128,:)=[linspace(0,1,128)',linspace(0,1,128)',linspace(0,1,128)'];
        
        load vik.mat
        vik128 = vik(round(linspace(1,128,128)),:);
        cmap_monoIncrease(256:-1:129,:) = vik128;
        
        modeInfo.cmap = cmap_monoIncrease;
        modeInfo.name = 'cmap_monoIncrease';
        VOLUME{1}.ui.(['co' 'Mode']) = modeInfo;
        VOLUME{1} = setCothresh(VOLUME{1}, minVE);VOLUME{1} = refreshScreen(VOLUME{1}, 1);
        VOLUME{1} = setClipMode(VOLUME{1}, 'co', [0 1]);
        [~, ~, ~, increaseColors] = meshColorOverlay(VOLUME{1},0);
        VOLUME{1} = meshColorOverlay(VOLUME{1});
        
        
        %% overlay mesh colors - slow loop, not vectorised but works
        sulciColor=[160;160;160;153];gyriColor=[96;96;96;153];
        tunedIncrease_colors = nan(4,size(tunedColors,2));
        
        tuned_sulcus = sum(tunedColors == sulciColor)==4;
        mono_sulcus = sum(increaseColors == sulciColor)==4;
        tuned_gyrus = sum(tunedColors == gyriColor)==4;
        mono_gyrus = sum(increaseColors == gyriColor)==4;
        sulci = find(tuned_sulcus & mono_sulcus);
        gyri = find(tuned_gyrus & mono_gyrus);
        tuned_voxel = find((~tuned_sulcus & mono_sulcus) | (~tuned_gyrus & mono_gyrus));
        mono_voxel = find((tuned_sulcus & ~mono_sulcus) | (tuned_gyrus & ~mono_gyrus));
        
        for idx =1:length(sulci), tunedIncrease_colors(:,sulci(idx)) = sulciColor; end
        for idx = 1:length(gyri), tunedIncrease_colors(:,gyri(idx)) = gyriColor; end
        tunedIncrease_colors(:,tuned_voxel) = tunedColors(:,tuned_voxel);
        tunedIncrease_colors(:,mono_voxel) = increaseColors(:,mono_voxel);
        
        VOLUME{1}.mesh{1}.colors = tunedIncrease_colors;
        mrmSet(VOLUME{1}.mesh{1},'colors',tunedIncrease_colors');
        
        
        save_name = 'monotonicTuned_colors(cv)_';
        cd(save_mesh), save(strcat(save_name,'voxel_selection=2_minVE=', string(minVE), '_',subjNames{subj},'.mat'),'tunedIncrease_colors');
        cd(paths{subj})

        % store back, left and right mesh views
        % make sure to hide cursor
        % % % meshStoreSettings( viewGet(VOLUME{1}, 'Mesh') );
        
        msh = viewGet(VOLUME{1},'currentmesh');
        lights = meshGet(msh,'lights');
        lights{1}.diffuse=[0 0 0];
        lights{1}.ambient=[0.5 0.5 0.5];
        host = meshGet(msh, 'host');
        windowID = meshGet(msh, 'windowID');
        distance=10;
        for n = 1:length(lights)
            L.actor = lights{n}.actor;
            L.origin = distance .* lights{n}.origin;
            lights{n} = mergeStructures(lights{n}, L);
            mrMesh(host, windowID, 'set', L);
        end
        msh = meshSet(msh, 'lights', lights);
        VOLUME{1} = viewSet(VOLUME{1}, 'mesh', msh);
        meshRetrieveSettings(viewGet(VOLUME{1}, 'CurMesh'), currentSettings(meshes));
        wSize=[1056,1855]; %@1920x1080
        mrmSet(msh,'windowSize',wSize(1),wSize(2));
        
        %save screenshot
        fname = fullfile(save_mesh,[saveNames{meshes}, '_', subjNames{subj},'.png']);
        udata.rgb = mrmGet(msh,'screenshot')/255;
        imwrite(udata.rgb, fname);
        set(gcf,'userdata',udata);
        
        allMeshes = viewGet(VOLUME{1}, 'allmeshes');
        mrmSet(allMeshes, 'closeall');
        h = findobj('Name', '3DWindow (mrMesh)');
        close(h);
        VOLUME{1}=rmfield(VOLUME{1},'mesh');
        
        close all
        mrvCleanWorkspace
    end
    
end
end

