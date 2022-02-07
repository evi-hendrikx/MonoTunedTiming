function monoTuned_make_ecc_plots_anonymous(save_dir_ecc, cv_data, minVE, timing_maps,minDuration,maxDuration,minPeriod,maxPeriod,medianPerSubj)
%% make ecc plots: plots the variance explained of all voxels (with VE>minVE for at least 1 model) against eccentricity (bin width 0.2, range 0-5.6)
% collects all suitable voxels across participants, hemispheres, and runs.
% Bins these voxels in bins with 0.2 width
% fits quadratic and sigmoid function to the datapoints
% outputs plots & fits of these functions to the data
%   save_dir_ecc: folder where eccentricity resuls are stored
%   cv_data: cross validated data created in cv_timing_data from ve_data
%   minVE: threshold for inclusion of the voxels. 0 in eccentricity analyses
%   timing maps: 1 if analyze timing maps; 0 if analyze visual field maps. Always 0 in eccentricity analyses
%   medianPerSubj: whether overall mean (0), median (1), 75th percentile (2) or 90th percentile (2) per bin is plotted
%   minDuration/ maxDuration/ minPeriod/ maxPeriod: most extreme preferred values that we allowed for the tuned model (otherwise its fit was set to 0)


modelFieldNames = fieldnames(cv_data);
subjNames = fieldnames(cv_data.(modelFieldNames{1}));
condNames = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{1})); % even & odd
ROILabels = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{6}).(condNames{1}).crossValidated); %possible, because S7 has all maps
ROILabels(contains(ROILabels,{'VO1','VO2','PHC'})) = []; %skip ventral stream ROIs 
ROIs = unique(erase(ROILabels, ["Right","Left","right","left"]),'stable');


%% combine all voxels and corresponding values from all subjects in an area
for roi = 1: length(ROILabels)
    
    %necessary to combine left and right hemispheres into the same group
    if roi <= length(ROIs)
        new_roi_name = strcat('Both', erase(erase(ROILabels{roi}, 'Left'),'left')); % the left hemisphere is first in ve_data. Change if monoTuned_load_timing_data changes
        group.ecc.(new_roi_name).ecc = [];
        group.cv_data.(modelFieldNames{1}).(new_roi_name) = [];
        group.cv_data.(modelFieldNames{2}).(new_roi_name) = [];
        group.cv_data.difference.(new_roi_name) = [];
    elseif roi > length(ROIs)
        new_roi_name = strcat('Both', erase(erase(ROILabels{roi}, 'Right'),'right'));
    end
    
    % group all voxels in one area together in "group"
    for subj = 1:length(subjNames)
        for condition = 1:length(condNames)
            roi_names_subj = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{subj}).(condNames{condition}).crossValidated); % model, subj, all (automatically also missing in odd and even)
            if any(strcmp(roi_names_subj, ROILabels{roi})) % some timing maps are missing
                
                eccCurrent = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{1})).crossValidated.(ROILabels{roi}).ecc;
                veCurrent_mod1 = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).varianceExplained;
                veCurrent_mod2 = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).varianceExplained;
                
                %necessary for voxel selection
                veCurrent_mod1_pre_cv = cv_data.(modelFieldNames{1}).(subjNames{subj}).(char(condNames{condition})).used_for_crossValidation.(ROILabels{roi}).varianceExplained;
                veCurrent_mod2_pre_cv = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).used_for_crossValidation.(ROILabels{roi}).varianceExplained;
                
                if find(strcmp(modelFieldNames,'TunedLin2d'))
                    % if parameters out of range for the tuned model: set VE to
                    % 0
                    xs_tuned_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).xs';
                    ys_tuned_Current = cv_data.(modelFieldNames{2}).(subjNames{subj}).(char(condNames{condition})).crossValidated.(ROILabels{roi}).ys';
                    outsideRange = find(sum([xs_tuned_Current<=minDuration, xs_tuned_Current>=maxDuration,ys_tuned_Current<=minPeriod,ys_tuned_Current>=maxPeriod],2)>0);
                    veCurrent_mod2(outsideRange)=0;
                    veCurrent_mod2_pre_cv(outsideRange)=0;
                end
                
                veCurrent_diff = veCurrent_mod2 - veCurrent_mod1;
                
                % I want a voxel removed if it performs below threshold for both of the models. It is allowed to be < treshold for one of them
                removeVE = find(sum([veCurrent_mod1_pre_cv' <=minVE,veCurrent_mod2_pre_cv'<= minVE],2)>1);
                
                veCurrent_mod1(removeVE) = [];
                veCurrent_mod2(removeVE) = [];
                veCurrent_diff(removeVE) = [];
                eccCurrent(removeVE) = [];
                
                group.ecc.(new_roi_name).ecc = [group.ecc.(new_roi_name).ecc, eccCurrent];
                group.cv_data.(modelFieldNames{1}).(new_roi_name) = [group.cv_data.(modelFieldNames{1}).(new_roi_name), veCurrent_mod1];
                group.cv_data.(modelFieldNames{2}).(new_roi_name) = [group.cv_data.(modelFieldNames{2}).(new_roi_name), veCurrent_mod2];
                group.cv_data.difference.(new_roi_name) = [group.cv_data.difference.(new_roi_name), veCurrent_diff];
                
            end
        end
    end
end

%% make bins
% very necessary! rest of the script runs with ROINames
ROINames = fieldnames(group.ecc);

% divide all grouped voxels into bins based on eccentricity
group.bins = [];
bins = [];
binsize = 0.20;
thresh.ecc = [0.1 5.5]; % center of the bins
voxelSize=1.77^2;
bins.(modelFieldNames{1}).x = (thresh.ecc(1):binsize:thresh.ecc(2))';
bins.(modelFieldNames{2}).x = (thresh.ecc(1):binsize:thresh.ecc(2))';
value_types = {modelFieldNames{1} modelFieldNames{2} 'difference'};

for roi = 1:length(ROINames)
    bin_id = 0;
    for b = thresh.ecc(1):binsize:thresh.ecc(2)
        
        % find voxels at specific eccentricities
        group_bii = group.ecc.(ROINames{roi}).ecc >=  b-binsize./2 & ...
            group.ecc.(ROINames{roi}).ecc < b+binsize./2;
        bin_id = bin_id + 1;
        
        for valtyp = 1:length(value_types)
            
            if sum(group_bii) >= 50
                % correcting for the fact that we measure voxels at a larger
                % size and do the analyses at 1mm^3
                
                s = wstat(group.cv_data.(value_types{valtyp}).(ROINames{roi})(group_bii),[], voxelSize^2);
                if medianPerSubj == 0
                    group.bins.(value_types{valtyp}).(ROINames{roi}).variance_measure(bin_id) = s.sterr;
                    group.bins.(value_types{valtyp}).(ROINames{roi}).summary_stat(bin_id) = s.mean;
                elseif medianPerSubj == 1
                    group.bins.(value_types{valtyp}).(ROINames{roi}).summary_stat(bin_id) = median(group.cv_data.(value_types{valtyp}).(ROINames{roi})(group_bii));
                    group.bins.(value_types{valtyp}).(ROINames{roi}).variance_measure(:,bin_id) = bootci(1000, {@median, group.cv_data.(value_types{valtyp}).(ROINames{roi})(group_bii)},'alpha',0.05);
                elseif medianPerSubj == 2
                    group.bins.(value_types{valtyp}).(ROINames{roi}).summary_stat(bin_id) = prctile(group.cv_data.(value_types{valtyp}).(ROINames{roi})(group_bii),75);
                    group.bins.(value_types{valtyp}).(ROINames{roi}).variance_measure(:,bin_id) = bootci(1000, {@prctile, group.cv_data.(value_types{valtyp}).(ROINames{roi})(group_bii),75},'alpha',0.05);
                elseif medianPerSubj == 3
                    group.bins.(value_types{valtyp}).(ROINames{roi}).summary_stat(bin_id) = prctile(group.cv_data.(value_types{valtyp}).(ROINames{roi})(group_bii),90);
                    group.bins.(value_types{valtyp}).(ROINames{roi}).variance_measure(:,bin_id) = bootci(1000, {@prctile, group.cv_data.(value_types{valtyp}).(ROINames{roi})(group_bii),90},'alpha',0.05);
                end
            else
                fprintf('Warning:No data in eccentricities %.1f to %.1f for roi %s.\n',...
                    b-binsize./2,b+binsize./2, ROINames{roi});
                group.bins.(value_types{valtyp}).(ROINames{roi}).summary_stat(bin_id) = 0;
                if medianPerSubj == 0
                    group.bins.(value_types{valtyp}).(ROINames{roi}).variance_measure(bin_id) = 0;
                else
                    group.bins.(value_types{valtyp}).(ROINames{roi}).variance_measure(:,bin_id)=[0,0];
                end
            end
        end
    end
    
end

close all

%% MAKE FIGURES
change_path =  [save_dir_ecc, 'ecc_plots'];
y_min_diff = -0.2;
y_max_diff = 0.2;
y_min_mod = 0;
y_max_mod = 0.4;
cd(change_path)

x_min = 0;
x_max = 5.6;
coloursDifference = [0.9290, 0.6940, 0.1250];
coloursmod2 = [205,0,0]./255;coloursmod1 = [0,0,205]./255;

%% quadratic fit
for roi = 1:length(ROINames)
    poly_order = 2;
    
    % prepare data
    x_diff = bins.(modelFieldNames{1}).x;
    x_mod1 = bins.(modelFieldNames{1}).x;
    x_mod2 = bins.(modelFieldNames{2}).x;
    y_difference =group.bins.difference.(ROINames{roi}).summary_stat';
    y_mod1 = group.bins.(modelFieldNames{1}).(ROINames{roi}).summary_stat';
    y_mod2 = group.bins.(modelFieldNames{2}).(ROINames{roi}).summary_stat';
    
    if medianPerSubj==0
        sterr_difference_low =  group.bins.difference.(ROINames{roi}).variance_measure';
        sterr_difference_high =  group.bins.difference.(ROINames{roi}).variance_measure';
        sterr_mod1_low = group.bins.(modelFieldNames{1}).(ROINames{roi}).variance_measure';
        sterr_mod1_high = group.bins.(modelFieldNames{1}).(ROINames{roi}).variance_measure';
        sterr_mod2_low = group.bins.(modelFieldNames{2}).(ROINames{roi}).variance_measure';
        sterr_mod2_high = group.bins.(modelFieldNames{2}).(ROINames{roi}).variance_measure';
    else
        sterr_difference_low =  y_difference - group.bins.difference.(ROINames{roi}).variance_measure(1,:)';
        sterr_difference_high =  group.bins.difference.(ROINames{roi}).variance_measure(2,:)' - y_difference;
        sterr_mod1_low = y_mod1 - group.bins.(modelFieldNames{1}).(ROINames{roi}).variance_measure(1,:)';
        sterr_mod1_high = group.bins.(modelFieldNames{1}).(ROINames{roi}).variance_measure(2,:)' - y_mod1;
        sterr_mod2_low = y_mod2 - group.bins.(modelFieldNames{2}).(ROINames{roi}).variance_measure(1,:)';
        sterr_mod2_high = group.bins.(modelFieldNames{2}).(ROINames{roi}).variance_measure(2,:)' - y_mod2;
    end
    
    
    % remove if entire bins were set to 0 (because of too little data)
    remove_bin = find(sum([y_mod1==0,y_mod2==0],2)>1);
    x_diff(remove_bin) = [];
    x_mod1(remove_bin)=[];
    x_mod2(remove_bin) = [];
    y_mod1(remove_bin) = [];
    y_mod2(remove_bin) = [];
    y_difference(remove_bin) = [];
    sterr_mod1_low(remove_bin) = [];
    sterr_mod1_high(remove_bin) = [];
    sterr_mod2_low(remove_bin) = [];
    sterr_mod2_high(remove_bin) = [];
    sterr_difference_low(remove_bin) = [];
    sterr_difference_high(remove_bin) = [];
    
    %% compute quaddratic values
    for valtyp =1:length(value_types)
        if valtyp == 1
            x_steps = x_mod1;
            y_vals = y_mod1;
        elseif valtyp == 2
            x_steps = x_mod2;
            y_vals = y_mod2;
        elseif valtyp == 3
            x_steps = x_diff;
            y_vals = y_difference;
        end
        
        % returns the coefficients for a polynomial p(x) of degree n that
        % is a best fit for the data y p = polyfit(x,y,n)
        [pQuad,SQuad] = polyfit(x_steps,y_vals,poly_order);
        
        %polyconf executes polynomial p at values in X
        % takes p and S from polyfit
        % calculates confidence intervals
        xvQuad = linspace(min(x_steps), max(x_steps), 1000);
        [data.quad.(value_types{valtyp}).(ROINames{roi}).CI,data.quad.(value_types{valtyp}).(ROINames{roi}).delta] = polyconf(pQuad,xvQuad,SQuad,'predopt','curve');
        
        % calculates values of best fitting polynomial at X (for plot)
        data.quad.(value_types{valtyp}).(ROINames{roi}).xfit = linspace(min(x_steps),max(x_steps),1000)';
        data.quad.(value_types{valtyp}).(ROINames{roi}).yfit = polyval(pQuad,data.quad.(value_types{valtyp}).(ROINames{roi}).xfit,SQuad);
        
        % calculates values of best fitting polynomial for each bin (for R2)
        data.quad.(value_types{valtyp}).(ROINames{roi}).xfitR2 = x_steps;
        data.quad.(value_types{valtyp}).(ROINames{roi}).yfitR2 = polyval(pQuad,data.quad.(value_types{valtyp}).(ROINames{roi}).xfitR2,SQuad);
        
        [R,P] = corrcoef(y_vals,data.quad.(value_types{valtyp}).(ROINames{roi}).yfitR2);
        data.quad.(value_types{valtyp}).(ROINames{roi}).Rsquared = R(2)^2;
        data.quad.(value_types{valtyp}).(ROINames{roi}).Rsquared_p = P(2);
        
        save_name = strcat(modelFieldNames{1}, '_',modelFieldNames{2},'_');
        if timing_maps == 1
            save_name = strcat(save_name, 'timing_maps_');
        end
        save_name = strcat('group_', save_name);
        save_name = strcat(save_name, 'afni2020(cv)_');
        if medianPerSubj == 1
            save_name = strcat(save_name, 'medianPerSubj');
        elseif medianPerSubj == 0
            save_name = strcat(save_name, 'meanPerSubj');
        elseif medianPerSubj == 2
            save_name = strcat(save_name, '75thPerSubj');
        elseif medianPerSubj == 3
            save_name = strcat(save_name, '90thPerSubj');
        end
        save_name = strcat(save_name, '_minVE=', string(minVE));
        save(save_name,'data');
    end
    
    %% plot quadratic figures
    for valtyp =1:length(value_types)
        if valtyp ~=2
            figure('visible','off');
        end
        hold on
        
        if valtyp == 1
            y_min = y_min_mod;
            y_max = y_max_mod;
            x_steps = x_mod1;
            y_vals = y_mod1;
            sterr_low = sterr_mod1_low;
            sterr_high = sterr_mod1_high;
            colours = coloursmod1;
        elseif valtyp == 2
            x_steps = x_mod2;
            y_vals = y_mod2;
            sterr_low = sterr_mod2_low;
            sterr_high = sterr_mod2_high;
            colours = coloursmod2;
        elseif valtyp == 3
            y_min = y_min_diff;
            y_max = y_max_diff;
            x_steps = x_diff;
            y_vals = y_difference;
            sterr_low = sterr_difference_low;
            sterr_high = sterr_difference_high;
            colours = coloursDifference;
        end
        
        line([x_min x_max],[0 0],'Color','black','LineStyle','-')
        plot(data.quad.(value_types{valtyp}).(ROINames{roi}).xfit, data.quad.(value_types{valtyp}).(ROINames{roi}).yfit, '-', 'Color',colours,'LineWidth',1);
        plot(data.quad.(value_types{valtyp}).(ROINames{roi}).xfit,data.quad.(value_types{valtyp}).(ROINames{roi}).CI-data.quad.(value_types{valtyp}).(ROINames{roi}).delta, '--','Color',colours,'LineWidth',0.75); % 95% confidence interval of fit line
        plot(data.quad.(value_types{valtyp}).(ROINames{roi}).xfit,data.quad.(value_types{valtyp}).(ROINames{roi}).CI+data.quad.(value_types{valtyp}).(ROINames{roi}).delta, '--','Color',colours,'LineWidth',0.75); % 95% confidence interval of fit line
   
       if valtyp == 3
            marker = 'o';
            size = 1.5;
        elseif valtyp == 2
            marker = 'd';
            size = 1.7;
        elseif valtyp == 1
            marker = 's';
            size = 1.9;
        end
        
        plot_bins = errorbar(x_steps, y_vals, sterr_low,sterr_high, 'vertical','LineStyle', 'none','Color',colours,'Marker',marker,'MarkerSize',size,'MarkerEdgeColor','k','MarkerFaceColor',colours,'LineWidth',0.5,'Displayname', value_types{valtyp}); % Original data
        
        %for legend
        if valtyp == 1
            mod1_bins = plot_bins;
        end
        
        drawnow;
        plot_bins.MarkerHandle.LineWidth = 0.5;

%         scatter(x_steps,y_vals,2,marker,'MarkerFaceColor',colours,'MarkerEdgeColor','k','HandleVisibility', 'off');
        
        if valtyp ~= 1
            if string(ROIs{roi}) == 'V1' && valtyp == 3
                legend(plot_bins(1),'Location', 'northwest');
            elseif string(ROIs{roi}) == 'V1' && valtyp == 2
                legend([mod1_bins(1), plot_bins(1)],'Location', 'northwest');
            end
            
            text(x_max - x_max/2-0.2, y_max - 0.05, ROIs{roi},'FontWeight', 'bold','FontSize', 8);
            xlabel('Eccentricity','FontWeight','bold');
            if valtyp == 3
                ylabel({'Difference in variance explained'},'FontWeight','bold');
            else
                ylabel({'Variance explained'},'FontWeight','bold');
            end
            yticks(y_min:0.05:y_max);
            yticklabels(y_min:0.05:y_max);
            axis([x_min x_max y_min y_max])
            axis square
            set(gcf,'color','w');
            set(gcf,'Units','centimeters','Position',[0.1 0.1 3.7 3.7]);
            
            box off
            
            save_name = strcat('Quad_',ROINames{roi},'_bootciDefault_');
            save_name = strcat(save_name, modelFieldNames{1}, '_',modelFieldNames{2},'_');
            if valtyp==3
                save_name = strcat(save_name, 'difference_');
            end
            if timing_maps == 1
                save_name = strcat(save_name, 'timing_maps_');
            end
            save_name = strcat('group_', save_name);
            save_name = strcat(save_name, 'afni2020(cv)_');
            if medianPerSubj == 1
                save_name = strcat(save_name, 'medianPerSubj');
            elseif medianPerSubj == 0
                save_name = strcat(save_name, 'meanPerSubj');
            elseif medianPerSubj == 2
                save_name = strcat(save_name, '75thPerSubj');
            elseif medianPerSubj == 3
                save_name = strcat(save_name, '90thPerSubj');
            end
            save_name = strcat(save_name, '_minVE=', string(minVE),'.eps');
            disp(save_name)
            export_fig(save_name,'-eps','-r600','-painters');
            close all
        end
    end
    
    
    %% Sigmoid fits
    %% compute sigmoid values
    for valtyp =1:length(value_types)
        if valtyp == 1
            x_steps = x_mod1;
            y_vals = y_mod1;
        elseif valtyp == 2
            x_steps = x_mod2;
            y_vals = y_mod2;
        elseif valtyp == 3
            x_steps = x_diff;
            y_vals = y_difference;
        end
        if valtyp == 3 && y_vals(1)<0
            %Pass in positive values of variance explained, by subtracting
            %from zero
        [data.sigm.(value_types{valtyp}).(ROINames{roi}).b_lower, data.sigm.(value_types{valtyp}).(ROINames{roi}).b_upper, ...
            data.sigm.(value_types{valtyp}).(ROINames{roi}).xfit, data.sigm.(value_types{valtyp}).(ROINames{roi}).yfit, ...
            data.sigm.(value_types{valtyp}).(ROINames{roi}).xfitR2, data.sigm.(value_types{valtyp}).(ROINames{roi}).yfitR2] = monoTuned_make_sigmoid_CIs_anonymous(x_steps, -y_vals);
            
            %Take the resulting fits and subtract from zero
            data.sigm.(value_types{valtyp}).(ROINames{roi}).b_lower= -data.sigm.(value_types{valtyp}).(ROINames{roi}).b_lower;
            data.sigm.(value_types{valtyp}).(ROINames{roi}).b_upper= -data.sigm.(value_types{valtyp}).(ROINames{roi}).b_upper;
            data.sigm.(value_types{valtyp}).(ROINames{roi}).yfit= -data.sigm.(value_types{valtyp}).(ROINames{roi}).yfit;
            data.sigm.(value_types{valtyp}).(ROINames{roi}).yfitR2= -data.sigm.(value_types{valtyp}).(ROINames{roi}).yfitR2;
        else
        [data.sigm.(value_types{valtyp}).(ROINames{roi}).b_lower, data.sigm.(value_types{valtyp}).(ROINames{roi}).b_upper, ...
            data.sigm.(value_types{valtyp}).(ROINames{roi}).xfit, data.sigm.(value_types{valtyp}).(ROINames{roi}).yfit, ...
            data.sigm.(value_types{valtyp}).(ROINames{roi}).xfitR2, data.sigm.(value_types{valtyp}).(ROINames{roi}).yfitR2] = monoTuned_make_sigmoid_CIs_anonymous(x_steps, y_vals);
        end
        %calculate fit
        [R,P] = corrcoef(y_vals,data.sigm.(value_types{valtyp}).(ROINames{roi}).yfitR2);
        data.sigm.(value_types{valtyp}).(ROINames{roi}).Rsquared = R(2)^2;
        data.sigm.(value_types{valtyp}).(ROINames{roi}).Rsquared_p = P(2);
    end
    
    save_name = strcat(modelFieldNames{1}, '_',modelFieldNames{2},'_');
    if timing_maps == 1
        save_name = strcat(save_name, 'timing_maps_');
    end
    save_name = strcat('group_', save_name);
    save_name = strcat(save_name, 'afni2020(cv)_');
    if medianPerSubj == 1
        save_name = strcat(save_name, 'medianPerSubj');
    elseif medianPerSubj == 0
        save_name = strcat(save_name, 'meanPerSubj');
    elseif medianPerSubj == 2
        save_name = strcat(save_name, '75thPerSubj');
    elseif medianPerSubj == 3
        save_name = strcat(save_name, '90thPerSubj');
    end
    save_name = strcat(save_name, '_minVE=', string(minVE));
    save(save_name,'data');
    
    
    %% plot sigmoid values
    for valtyp =1:length(value_types)
        if valtyp ~=2
            figure('visible','off');
        end
        hold on
        
        if valtyp == 1
            y_min = y_min_mod;
            y_max = y_max_mod;
            x_steps = x_mod1;
            y_vals = y_mod1;
            sterr_low = sterr_mod1_low;
            sterr_high = sterr_mod1_high;
            colours = coloursmod1;
        elseif valtyp == 2
            x_steps = x_mod2;
            y_vals = y_mod2;
            sterr_low = sterr_mod2_low;
            sterr_high = sterr_mod2_high;
            colours = coloursmod2;
        elseif valtyp == 3
            y_min = y_min_diff;
            y_max = y_max_diff;
            x_steps = x_diff;
            y_vals = y_difference;
            sterr_low = sterr_difference_low;
            sterr_high = sterr_difference_high;
            colours = coloursDifference;
        end
        
        line([x_min x_max],[0 0],'Color','black','LineStyle','-')
        plot(data.sigm.(value_types{valtyp}).(ROINames{roi}).xfit, data.sigm.(value_types{valtyp}).(ROINames{roi}).yfit, '-', 'Color',colours,'LineWidth',1);
        plot(data.sigm.(value_types{valtyp}).(ROINames{roi}).xfit, data.sigm.(value_types{valtyp}).(ROINames{roi}).b_upper, '--','Color',colours,'LineWidth',0.75);
        plot(data.sigm.(value_types{valtyp}).(ROINames{roi}).xfit, data.sigm.(value_types{valtyp}).(ROINames{roi}).b_lower,  '--','Color',colours,'LineWidth',0.75);
        
        if valtyp == 3
            marker = 'o';
            size = 1.5;
        elseif valtyp == 2
            marker = 'd';
            size = 1.7;
        elseif valtyp == 1
            marker = 's';
            size = 1.9;
        end
        
        plot_bins = errorbar(x_steps, y_vals, sterr_low,sterr_high, 'vertical','LineStyle', 'none','Color',colours,'Marker',marker,'MarkerSize',size,'MarkerEdgeColor','k','MarkerFaceColor',colours,'LineWidth',0.5,'Displayname', value_types{valtyp}); % Original data
        %for legend
        if valtyp == 1
            mod1_bins = plot_bins;
        end
        drawnow;
        plot_bins.MarkerHandle.LineWidth = 0.5;
%         scatter(x_steps,y_vals,2,marker,'MarkerFaceColor',colours,'MarkerEdgeColor','k','HandleVisibility', 'off');
        
        if valtyp ~= 1
            if string(ROIs{roi}) == 'V1' && valtyp == 3
                legend(plot_bins(1),'Location', 'northwest');
            elseif string(ROIs{roi}) == 'V1' && valtyp == 2
                legend([mod1_bins(1), plot_bins(1)],'Location', 'northwest');
            end
            
            text(x_max - x_max/2-0.2, y_max - 0.05, ROIs{roi},'FontWeight', 'bold','FontSize', 8);
            xlabel('Eccentricity','FontWeight','bold');
            if valtyp == 3
                ylabel({'Difference in variance explained'},'FontWeight','bold');
            else
                ylabel({'Variance explained'},'FontWeight','bold');
            end
            yticks(y_min:0.05:y_max);
            yticklabels(y_min:0.05:y_max);
            axis([x_min x_max y_min y_max])
            axis square
            set(gcf,'color','w');
            set(gcf,'Units','centimeters','Position',[0.1 0.1 3.7 3.7]);
            
            box off
            
            save_name = strcat('Sigm_',ROINames{roi},'_bootciDefault_');
            save_name = strcat(save_name, modelFieldNames{1}, '_',modelFieldNames{2},'_');
            if valtyp ==3
                save_name = strcat(save_name, 'difference_');
            end
            if timing_maps == 1
                save_name = strcat(save_name, 'timing_maps_');
            end
            save_name = strcat('group_', save_name);
            save_name = strcat(save_name, 'afni2020(cv)_');
            if medianPerSubj == 1
                save_name = strcat(save_name, 'medianPerSubj');
            elseif medianPerSubj == 0
                save_name = strcat(save_name, 'meanPerSubj');
            elseif medianPerSubj == 2
                save_name = strcat(save_name, '75thPerSubj');
            elseif medianPerSubj == 3
                save_name = strcat(save_name, '90thPerSubj');
            end
            save_name = strcat(save_name, '_minVE=', string(minVE), '.eps');
            disp(save_name)
            export_fig(save_name,'-eps','-r600','-painters');
            close all
        end
    end
end

end