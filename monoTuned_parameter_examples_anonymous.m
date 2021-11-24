function monoTuned_parameter_examples_anonymous(save_path, stat, timing_maps, medianPerSubj, stacked_plots_median, minVE)

save_path_model_properties = strcat(save_path, 'model_property_comparisons/');
cd(save_path_model_properties)

modelFieldNames=fieldnames(stat.xs)
ROIs = fieldnames(stat.xs.(modelFieldNames{1}).data)

if timing_maps ==1
    scatter_order=[4 10 9 8 2 6 5 7 1 3];
elseif timing_maps == 0
    scatter_order=[14:16,10:13,17,1:9];
end
selected_rois = [14 13 4 5 8];

if find(strcmp(modelFieldNames,'TunedLin2d'))
    parameters = {'xs','ys', 'beta_ratio','major','minor','theta','exp','ratio'};
else
    parameters = {'xs','ys','beta_ratio'};
end

for parameter = 1:length(parameters)
    
    if parameter == 1 || parameter == 2
        modelNames = {modelFieldNames{1}, modelFieldNames{2}};
    elseif parameter == 3
        modelNames = {modelFieldNames{1}};
    elseif parameter > 3 && find(strcmp(modelFieldNames,'TunedLin2d'))
        % some analyses can only be done for tuned
        modelNames = {modelFieldNames{2}};
    end
    
    for model =1:length(modelNames)
        % use NaN removed data
        for roi = 1:length(ROIs)
            medians.(modelNames{model}).(ROIs{roi}).(parameters{parameter}) =median(stat.(parameters{parameter}).(modelNames{model}).data.(ROIs{roi}));
        end
    end
end

%% timing info
modelPeriods = [0.0000001 0.005:0.005:1.05 2.055:0.005:2.1];
modelDurations = [0:0.005:1.05 1.955:0.005:2];
modelFrequencies = 1./modelPeriods';

modelDurationsGrid = modelDurations.*ones(size(modelPeriods))';
modelPeriodGrid = ones(size(modelDurations)).*(modelPeriods');
invalidGrid = modelDurationsGrid > modelPeriodGrid;

for roi =  1:length(ROIs)
    %% 2D plane mono
    median_xs =  medians.(modelFieldNames{1}).(ROIs{roi}).xs;
    median_ys =  medians.(modelFieldNames{1}).(ROIs{roi}).ys;
    median_ratio =  medians.(modelFieldNames{1}).(ROIs{roi}).beta_ratio;
    
   
    amplitudeDuration = (modelDurations.^median_xs).*ones(size(modelPeriods))';
    amplitudeFrequency = ones(size(modelDurations)).*((modelFrequencies.^median_ys)./modelFrequencies);
    amplitudeMono = amplitudeDuration* median_ratio + amplitudeFrequency;
    amplitudeMono(invalidGrid)=min(min(amplitudeMono(~invalidGrid)));
    
    figure; imagesc(flipud(amplitudeMono)); axis image
    
    
    savename = 'mono_2D_';
    savename = strcat(savename, modelFieldNames{1}, '_', modelFieldNames{2}, '_');
    savename = strcat(savename, ROIs{roi}, '_');
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
    if timing_maps == 1
        savename = strcat(savename, '_timing_maps');
    elseif timing_maps == 0
        savename = strcat(savename, '_visual_field_maps');
    end
    savename = strcat(savename, '_minVE=', string(minVE));
    savename_eps = strcat(savename, '.eps');
    export_fig(savename_eps,'-eps','-r600','-painters');

    close all

 
    
    %% 2D plane tuned
    median_xs =  medians.(modelFieldNames{2}).(ROIs{roi}).xs;
    median_ys =  medians.(modelFieldNames{2}).(ROIs{roi}).ys;
    median_major =  medians.(modelFieldNames{2}).(ROIs{roi}).major;
    median_minor =  medians.(modelFieldNames{2}).(ROIs{roi}).minor;
    median_theta =  deg2rad(medians.(modelFieldNames{2}).(ROIs{roi}).theta);
    median_exp = medians.(modelFieldNames{2}).(ROIs{roi}).exp;
% %     
    % first image without compressive exponent on frequency
    rfTunedImg1= rfGaussian2d(modelDurations, modelPeriods',...
        median_major, ...
        median_minor, ...
        median_theta, ...
        median_xs , ...
        median_ys);
        
    freq=5./modelPeriods';
    rfTunedImg2 = rfTunedImg1.*((freq.^median_exp)./freq);
    
    rfTunedImg1(invalidGrid) = min(min(rfTunedImg1(~invalidGrid)));
    rfTunedImg2(invalidGrid) = min(min(rfTunedImg2(~invalidGrid)));

    figure; imagesc(flipud(rfTunedImg1)); axis image
    
    savename = 'tuned_2D_';
    savename = strcat(savename, modelFieldNames{1}, '_', modelFieldNames{2}, '_');
    savename = strcat(savename, ROIs{roi}, '_');
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
    if timing_maps == 1
        savename = strcat(savename, '_timing_maps');
    elseif timing_maps == 0
        savename = strcat(savename, '_visual_field_maps');
    end
    savename = strcat(savename, '_minVE=', string(minVE));
    savename_eps = strcat(savename, '.eps');
    export_fig(savename_eps,'-eps','-r600','-painters');

    close all
    
    %% exponent tuned
    frequencies=0.01:0.01:20;
    response=frequencies.^median_exp;
    figure; plot(frequencies, response)
    axis square;
    axis([0 20 0 20]);
    xlabel('Frequency (Hz)')
    ylabel('Response amplitude')
    
    
    savename = 'tuned_exponent_';
    savename = strcat(savename, modelFieldNames{1}, '_', modelFieldNames{2}, '_');
    savename = strcat(savename, ROIs{roi}, '_');
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
    if timing_maps == 1
        savename = strcat(savename, '_timing_maps');
    elseif timing_maps == 0
        savename = strcat(savename, '_visual_field_maps');
    end
    savename = strcat(savename, '_minVE=', string(minVE));
    savename_eps = strcat(savename, '.eps');
    export_fig(savename_eps,'-eps','-r600','-painters');

    close all
    
    %% scale 2D with the frequency component
  
    figure; imagesc(flipud(rfTunedImg2)); axis image
    
    savename = 'tuned_combined_2D_';
    savename = strcat(savename, modelFieldNames{1}, '_', modelFieldNames{2}, '_');
    savename = strcat(savename, ROIs{roi}, '_');
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
    if timing_maps == 1
        savename = strcat(savename, '_timing_maps');
    elseif timing_maps == 0
        savename = strcat(savename, '_visual_field_maps');
    end
    savename = strcat(savename, '_minVE=', string(minVE));
    savename_eps = strcat(savename, '.eps');
    export_fig(savename_eps,'-eps','-r600','-painters');
    
    close all
end
end

