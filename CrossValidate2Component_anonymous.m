function CrossValidate2Component_anonymous(paths, whichModelMono, combinedDTMono, validation, whichModelTuned, combinedDTTuned)
%% CrossValidate2Component: run models, copy the files to the other cross-validation run, and run cross-validation
% paths: where the information about each participant is stored
% whichModel: which model you want to run
% % whichModel = [13]; % compressive duration compressive frequency
% % whichModel = [15]; % linear duration compressive frequency
% % whichModel = [16]; % linear duration linear frequency
% combinedDTMono: which cross-validation runs
% % combinedDTMono=5:7; or combinedDTMono = 8:9; (validation monotonic)
% validation: run on self-made ground-truth dataset (1) or not (0)
% whichModelTuned: whichModel for the tuned model, we used 4
% combinedDTTuned: combinedDTTuned=10:13; or combinedDTTuned=14:17; (tuned) for validation


whichSubs=1:length(paths);
if validation == 0 % currently only does monotonic, as we already had the tuned files
    %Data types for all, odd and even runs
    
    %Fitting
    for thisSub=1:length(whichSubs)
        cd(paths{whichSubs(thisSub)})
        mrVista 3;

        %The code to run for each
        if whichSubs(thisSub)<=6
            load('TimingModelParamsMix.mat')
            setAllRetParams(paramsDurationPeriod, combinedDTMono);
            rmRunDurFreq2d(VOLUME{1},combinedDTMono, 'gray-Layer1',whichModelMono,{'1g'},[],[],[],0);
        else
            load('TimingModelParamsSessions.mat')
            setAllRetParams(paramsDurationPeriod, combinedDTMono);
            rmRunDurFreq2d(VOLUME{1},combinedDTMono, 'gray-Layer1',whichModelMono,{'1g'},[],[],[],1);
        end
        close(1); mrvCleanWorkspace;
    end
    
    
    % move around files
    for thisSub=whichSubs
        cd(paths{thisSub})
        mrVista 3;
        
        %The code to run for each
        if thisSub<=6
            load('TimingModelParamsMix.mat')
        else
            load('TimingModelParamsSessions.mat')
        end
        allXvalDTs=combinedDTMono(2:3);
        
        %Monotonic component models
        for n=1:length(allXvalDTs)
            files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/', '*YNoNormOccupancy*.mat']);
            thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/'];
            otherPath=['Gray/' dataTYPES(allXvalDTs(3-n)).name, '/Monotonic/'];
            %        eval(['!mkdir ',  '"',otherPath, 'xval"']);
            %        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
            
            for whichFile=1:length(files)
                eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
            end
        end
        close(1); mrvCleanWorkspace;
        
    end
    
    
    
    %Cross validating the resulting model.
    for thisSub=whichSubs
        cd(paths{thisSub})
        mrVista 3;
        
        %The code to run for each
        if thisSub<=6
            load('TimingModelParamsMix.mat')
        else
            load('TimingModelParamsSessions.mat')
        end
        allXvalDTs=combinedDTMono(2:3);
        
        %
        for whichDT=1:length(allXvalDTs)
            folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/Monotonic/xval'];
            modelFiles=dir([folderName '/*202109*YNoNormOccupancy*']);
            for whichModel=1:length(modelFiles)
                rmMainPostSearch([1 allXvalDTs(whichDT)],'gray-Layer1',4, [folderName '/' modelFiles(whichModel).name]);
            end
        end
        delete(gcp('nocreate'))
        
        close(1); mrvCleanWorkspace;
        
    end

else
    
    %Both models on monotonic predictions
    for thisSub=1:length(whichSubs)
        cd(paths{whichSubs(thisSub)})
        mrVista 3;
        
        %The code to run for each
        load('TimingModelParamsMix.mat')
        paramsDurationPeriod(1).nDCT=0;
        paramsDurationPeriod(2).nDCT=0;
        paramsDurationPeriod(3).nDCT=0;
        paramsDurationPeriod(4).nDCT=0;
        setAllRetParams(paramsDurationPeriod, combinedDTMono);
        rmRunDurFreq2d(VOLUME{1},combinedDTMono, 'ValidationMonotonic',whichModelMono,{'1g'},[],[],[],0);
        rmRunDuration2dOval(VOLUME{1},combinedDTMono, 'ValidationMonotonic',whichModelTuned,{'1g'},[],[],[],0,'free');
        close(1); mrvCleanWorkspace;
        
    end
    
    % Both models on tuned predictions
    for thisSub=1:length(whichSubs)
        cd(paths{whichSubs(thisSub)})
        mrVista 3;
        
        %The code to run for each
        load('TimingModelParamsMix.mat')
        paramsDurationPeriod(1).nDCT=0;
        paramsDurationPeriod(2).nDCT=0;
        paramsDurationPeriod(3).nDCT=0;
        paramsDurationPeriod(4).nDCT=0;
        setAllRetParams(paramsDurationPeriod, combinedDTTuned);
        rmRunDurFreq2d(VOLUME{1},combinedDTTuned, 'ValidationTuned',whichModelMono,{'1g'},[],[],[],0);
        rmRunDuration2dOval(VOLUME{1},combinedDTTuned, 'ValidationTuned',whichModelTuned,{'1g'},[],[],[],0,'free');
        close(1); mrvCleanWorkspace;
    end
    
    
    %Now cross-validate
    % move around files
    for thisSub=whichSubs
        cd(paths{thisSub})
        mrVista 3;
        combinedDT=[combinedDTMono combinedDTTuned];
        
        %The code to run for each
        load('TimingModelParamsMix.mat')
        
        allXvalDTs=combinedDT;
        
        %Monotonic component models
        for n=1:length(allXvalDTs)
            files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/', '*compressiveXcompressiveYNoNormOccupancy*.mat']);
            thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/Monotonic/'];
            if mod(n,2)==1
                otherDTind=n+1;
            else
                otherDTind=n-1;
            end
            
            otherPath=['Gray/' dataTYPES(allXvalDTs(otherDTind)).name, '/Monotonic/'];
            %        eval(['!mkdir ',  '"',otherPath, 'xval"']);
            %        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
            
            for whichFile=1:length(files)
                eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
            end
            
            
            files=dir(['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFitFreeExponent/', '*Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit-fFit.mat']);
            thisPath=['Gray/' dataTYPES(allXvalDTs(n)).name, '/SearchFitFreeExponent/'];
            otherPath=['Gray/' dataTYPES(allXvalDTs(otherDTind)).name, '/SearchFitFreeExponent/'];
            %        eval(['!mkdir ',  '"',otherPath, 'xval"']);
            %        eval(['!mkdir ',  '"',otherPath, 'xvalRefit"']);
            
            for whichFile=1:length(files)
                eval(['!cp ', '"', thisPath, files(whichFile).name, '" "', otherPath, 'xval/xval-', files(whichFile).name, '"']);
            end
            
        end
        close(1); mrvCleanWorkspace;
        
    end
    
    
    % %
    % %     %Cross validating the resulting model. For models with detrending
    % %     for thisSub=whichSubs
    % %         cd(paths{thisSub})
    % %         mrVista 3;
    % %
    % %         %The code to run for each
    % %         load('TimingModelParamsMix.mat')
    % %         allXvalDTs=combinedDT;
    % %
    % %         %  Monotonic models
    % %         for whichDT=1:length(allXvalDTs)
    % %
    % %             folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/Monotonic/xval'];
    % %             modelFiles=dir([folderName '/*compressiveXcompressiveYNoNormOccupancy-DurFreq-DT0.5-20*']);
    % %             for whichModel=1:length(modelFiles)
    % %                 if whichDT<=2
    % %                     rmMainPostSearch([1 allXvalDTs(whichDT)],'ValidationMonotonic',4, [folderName '/' modelFiles(whichModel).name]);
    % %                 else
    % %                     rmMainPostSearch([1 allXvalDTs(whichDT)],'ValidationTuned',4, [folderName '/' modelFiles(whichModel).name]);
    % %                 end
    % %             end
    % %
    % %             %Tuned models
    % %             folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFitFreeExponent/xval'];
    % %             modelFiles=dir([folderName '/*Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit-fFit.mat']);
    % %             for whichModel=1:length(modelFiles)
    % %                 if whichDT<=2
    % %                     rmMainPostSearch([1 allXvalDTs(whichDT)],'ValidationMonotonic',4, [folderName '/' modelFiles(whichModel).name]);
    % %                 else
    % %                     rmMainPostSearch([1 allXvalDTs(whichDT)],'ValidationTuned',4, [folderName '/' modelFiles(whichModel).name]);
    % %                 end
    % %             end
    % %
    % %         end
    % %         delete(gcp('nocreate'))
    % %
    % %         close(1); mrvCleanWorkspace;
    % %     end
    % %
    
    %Cross validating the resulting model. For models without detrending
    for thisSub=whichSubs
        cd(paths{thisSub})
        mrVista 3;
        
        %The code to run for each
        load('TimingModelParamsMix.mat')
        allXvalDTs=combinedDT;
        
        %  Monotonic models
        for whichDT=1:length(allXvalDTs)
            
            folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/Monotonic/xval'];
            modelFiles=dir([folderName '/*compressiveXcompressiveYNoNormOccupancy-DurFreq-DT0-20*']);
            for whichModel=1:length(modelFiles)
                if whichDT<=2
                    rmMainPostSearch([1 allXvalDTs(whichDT)],'ValidationMonotonic',4, [folderName '/' modelFiles(whichModel).name]);
                else
                    rmMainPostSearch([1 allXvalDTs(whichDT)],'ValidationTuned',4, [folderName '/' modelFiles(whichModel).name]);
                end
            end
            
            %Tuned models
            folderName=[pwd '/Gray/' dataTYPES(allXvalDTs(whichDT)).name, '/SearchFitFreeExponent/xval'];
            modelFiles=dir([folderName '/*Lin-2dOvalGaussian-DurationPeriod-DT0-maxValue-2-expIntensity-free-fFit-fFit.mat']);
            for whichModel=1:length(modelFiles)
                if whichDT<=2
                    rmMainPostSearch([1 allXvalDTs(whichDT)],'ValidationMonotonic',4, [folderName '/' modelFiles(whichModel).name]);
                else
                    rmMainPostSearch([1 allXvalDTs(whichDT)],'ValidationTuned',4, [folderName '/' modelFiles(whichModel).name]);
                end
            end
            
        end
        delete(gcp('nocreate'))
        
        close(1); mrvCleanWorkspace;
    end
end
end