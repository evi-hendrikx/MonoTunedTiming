function [actualValues, evalModels] = monoTuned_load_test_model_data_anonymous(paths, subj_id, save_path, minDuration,maxDuration,minPeriod,maxPeriod, exclude_outside_range)
%% Loads model outputs on ground truth dataset into structs
% paths: where the information about each participant is stored 
% subj_id: index of subject you want (where in paths)
% save_path: general folder where you want results to be stored
% minDuration,maxDuration,minPeriod,maxPeriod: ranges outside which the tuned model is considered invalid


%% Get the necessary data
ROIsizeTuned=362700;
ROIsizeMono=147620;
monoDTOdd=8;
monoDTEven=9;
tunedDTOdd=10:2:16;
tunedDTEven=11:2:17;

datasets = {'Dmono', 'Dtuned'};
respModels = {'Mmono', 'Mtuned'};
cv_runs = {'Odd','Even'};

evalModels = {};
actualValues = {};

crossVal_or_not = {'crossValidated','used_for_crossValidation'};

cd(paths{subj_id})
mrVista 3;
for dset = 1:length(datasets) % the actual designed dataset
    for respModel = 1:length(respModels) % the model that was used for evaluation
        for crossVal = 1:length(crossVal_or_not)

            if dset == 1
                ROIsize = ROIsizeMono;
                if respModel == 1  && crossVal == 1% only needs to be done once
                    load(strcat(paths{subj_id},'/ValidationMonotonic3.mat'), 'params');
                    actualValues.(datasets{dset}).X = params.analysis.x0;
                    actualValues.(datasets{dset}).Y = params.analysis.y0;
                    actualValues.(datasets{dset}).Beta = params.analysis.ratio;
                    actualValues.(datasets{dset}).Noise = params.analysis.noiseSDs;
                end

            elseif dset == 2 && crossVal == 1
                ROIsize = ROIsizeTuned;
                if respModel == 1 % only needs to be done once
                    load(strcat(paths{subj_id},'/ValidationTuned.mat'), 'params');
                    actualValues.(datasets{dset}).X = params.analysis.x0;
                    actualValues.(datasets{dset}).Y = params.analysis.y0;
                    actualValues.(datasets{dset}).sMaj = params.analysis.sigmaMajor;
                    actualValues.(datasets{dset}).sMin = params.analysis.sigmaMinor;
                    actualValues.(datasets{dset}).Theta = params.analysis.theta;
                    actualValues.(datasets{dset}).Exp = params.analysis.exponent;
                    actualValues.(datasets{dset}).Noise = params.analysis.noiseSDs;


                    %Adjust for cases where minor axis is longer than major
                    %axis
                    tmp=actualValues.(datasets{dset}).sMaj;
                    actualValues.(datasets{dset}).sMaj(params.analysis.sigmaMinor>params.analysis.sigmaMajor)=actualValues.(datasets{dset}).sMin(params.analysis.sigmaMinor>params.analysis.sigmaMajor);
                    actualValues.(datasets{dset}).sMin(params.analysis.sigmaMinor>params.analysis.sigmaMajor)=tmp(params.analysis.sigmaMinor>params.analysis.sigmaMajor);
                    actualValues.(datasets{dset}).Theta(params.analysis.sigmaMinor>params.analysis.sigmaMajor) = actualValues.(datasets{dset}).Theta(params.analysis.sigmaMinor>params.analysis.sigmaMajor)+pi/2;
                    actualValues.(datasets{dset}).Theta(actualValues.(datasets{dset}).Theta>pi)=actualValues.(datasets{dset}).Theta(actualValues.(datasets{dset}).Theta>pi)-pi;

                    if exclude_outside_range == 1
                        outsideRange = sum([actualValues.(datasets{dset}).X<=minDuration, actualValues.(datasets{dset}).X>=maxDuration,actualValues.(datasets{dset}).Y<=minPeriod, actualValues.(datasets{dset}).Y>=maxPeriod],2)>0;
                    else
                        outsideRange = [];
                    end

                    actualValues.(datasets{dset}).X(outsideRange) = [];
                    actualValues.(datasets{dset}).Y(outsideRange) = [];
                    actualValues.(datasets{dset}).sMaj(outsideRange) = [];
                    actualValues.(datasets{dset}).sMin(outsideRange) = [];
                    actualValues.(datasets{dset}).Theta(outsideRange) = [];
                    actualValues.(datasets{dset}).Exp(outsideRange) = [];
                    actualValues.(datasets{dset}).Noise(outsideRange) = [];

                end
            end

            for cv = 1:length(cv_runs)
                if cv == 1 && dset == 1
                    DT_run = monoDTOdd;
                elseif cv == 1 && dset == 2
                    DT_run = tunedDTOdd;
                elseif cv == 2 && dset == 1
                    DT_run = monoDTEven;
                elseif cv == 2 && dset == 2
                    DT_run = tunedDTEven;
                end

                %Monotonic model
                for whichDT=1:length(DT_run)
                    if respModel == 1
                        if crossVal == 1
                            folderName=[pwd '/Gray/' dataTYPES(DT_run(whichDT)).name, '/Monotonic/xvalRefit'];
                        else
                            folderName=[pwd '/Gray/' dataTYPES(DT_run(whichDT)).name, '/Monotonic/xval'];
                        end
                        modelFiles=dir([folderName '/*compressiveXcompressiveYNoNormOccupancy-DurFreq-DT0-*']);
                    elseif respModel == 2
                        if crossVal == 1
                            folderName=[pwd '/Gray/' dataTYPES(DT_run(whichDT)).name, '/SearchFitFreeExponent/xvalRefit'];
                        else
                            folderName=[pwd '/Gray/' dataTYPES(DT_run(whichDT)).name, '/SearchFitFreeExponent/xval'];
                        end
                        modelFiles=dir([folderName '/*Lin-2dOvalGaussian-DurationPeriod-DT0-maxValue-2-expIntensity-free-fFit*']);

                    end

                    load([folderName '/' modelFiles.name], 'model');

                    if crossVal == 1
                        if whichDT == 1
                            evalModels.(datasets{dset}).(respModels{respModel}).X.(cv_runs{cv}) = [];
                            evalModels.(datasets{dset}).(respModels{respModel}).Y.(cv_runs{cv}) = [];
                            evalModels.(datasets{dset}).(respModels{respModel}).VE.(cv_runs{cv}) = [];
                            if respModel == 1
                                evalModels.(datasets{dset}).(respModels{respModel}).Beta.(cv_runs{cv})=[];
                            elseif respModel ==2
                                evalModels.(datasets{dset}).(respModels{respModel}).sMaj.(cv_runs{cv})=[];
                                evalModels.(datasets{dset}).(respModels{respModel}).sMin.(cv_runs{cv})=[];
                                evalModels.(datasets{dset}).(respModels{respModel}).Theta.(cv_runs{cv})=[];
                                evalModels.(datasets{dset}).(respModels{respModel}).Exp.(cv_runs{cv})=[];
                            end
                        end

                        get_voxels = 1:ROIsize;

                        if dset == 2 && exclude_outside_range == 1
                            get_voxels(outsideRange(1+ROIsize*(whichDT-1):ROIsize*whichDT))=[];
                        end

                        tmp=rmGet(model{1}, 'x0');
                        evalModels.(datasets{dset}).(respModels{respModel}).X.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).X.(cv_runs{cv}), tmp(get_voxels)];
                        tmp=rmGet(model{1}, 'y0');
                        evalModels.(datasets{dset}).(respModels{respModel}).Y.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).Y.(cv_runs{cv}), tmp(get_voxels)];
                        tmp=rmGet(model{1}, 've');
                        evalModels.(datasets{dset}).(respModels{respModel}).VE.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).VE.(cv_runs{cv}), tmp(get_voxels)];
                        if respModel == 1
                            tmp=rmGet(model{1}, 'b');
                            evalModels.(datasets{dset}).(respModels{respModel}).Beta.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).Beta.(cv_runs{cv}), tmp(1,get_voxels, 1:2)];
                        elseif respModel == 2
                            tmp=rmGet(model{1}, 'sigmamajor');
                            evalModels.(datasets{dset}).(respModels{respModel}).sMaj.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).sMaj.(cv_runs{cv}), tmp(get_voxels)];
                            tmp=rmGet(model{1}, 'sigmaminor');
                            evalModels.(datasets{dset}).(respModels{respModel}).sMin.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).sMin.(cv_runs{cv}), tmp(get_voxels)];
                            tmp=rmGet(model{1}, 'sigmatheta');
                            evalModels.(datasets{dset}).(respModels{respModel}).Theta.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).Theta.(cv_runs{cv}), tmp(get_voxels)];
                            tmp=rmGet(model{1}, 'exponent');
                            evalModels.(datasets{dset}).(respModels{respModel}).Exp.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).Exp.(cv_runs{cv}), tmp(get_voxels)];
                        end
                    else
                        if whichDT == 1
                            evalModels.(datasets{dset}).(respModels{respModel}).not_crossvalidated_VE.(cv_runs{cv}) = [];
                        end

                        get_voxels = 1:ROIsize;

                        if dset == 2 && exclude_outside_range == 1
                            get_voxels(outsideRange(1+ROIsize*(whichDT-1):ROIsize*whichDT))=[];
                        end

                        tmp=rmGet(model{1}, 've');
                        evalModels.(datasets{dset}).(respModels{respModel}).not_crossvalidated_VE.(cv_runs{cv})=[evalModels.(datasets{dset}).(respModels{respModel}).not_crossvalidated_VE.(cv_runs{cv}), tmp(get_voxels)];
                    end
                end
            end
        end
    end
end
save('parameters_validation.mat','evalModels','actualValues');

close all
end