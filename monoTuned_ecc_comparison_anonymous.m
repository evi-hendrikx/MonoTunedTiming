function monoTuned_ecc_comparison_anonymous(save_path, minVE, timing_maps,medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod)
%% ecc comparisons:  plots the variance explained against eccentricity for tuned & monotonic model and their difference. Also statistically compares near & far eccentricity ranges for tuned & monotonic & difference
%Organizes variance explained data from different models, subjects, conditions
%and ROIs into one structure array.
% plots the variance explained of all available voxels against eccentricity (in bins of 0.2 degree)
% Also, of voxels with a variance explained above the threshold (minVE) for at least one of the models (monotonic/ tuned)
% it compares the variance explained of voxels with a preferred eccentricity <1 to  voxels with a preferred eccentricity >2
% means are created on the level of hemisphere and cross-validation run.
% i.e., 8 participants x 2 hemispheres x 2 runs = 32 datapoints per ROI
% Runs ANOVA, post-hocs, and gives graphs
%   save_path: general folder where you want results to be stored
%   minVE: threshold for inclusion of the voxels: 0 in eccentricity analyses
%   timing maps: 1 if analyze timing maps; 0 if analyze visual field maps --> chose to never do timing maps, since it's only a subset of all possible
%       voxels that have a timing response and a preferred eccentricity
%   medianPerSubj: whether means (0), medians (1), 75th percentile (2) or 90th percentile (2) for each subject are included in the analyses 
%   minDuration/ maxDuration/ minPeriod/ maxPeriod: most extreme preferred values that we allowed for the tuned model (otherwise its fit was set to 0)

%% load data into right format
%% load and prepare data
if timing_maps == 1
    load('parameters_timing_maps.mat')
elseif timing_maps == 0
    load('parameters_visual_field_maps.mat')
end


%choose which models you want to use
modelFieldNames = fieldnames(ve_data);
use_models = {modelFieldNames{1},modelFieldNames{2}}; % fill in models you want

cv_data = monoTuned_cv_timing_data_anonymous(ve_data, use_models);
modelFieldNames = fieldnames(cv_data);
subjNames = fieldnames(cv_data.(modelFieldNames{1}));
condNames = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{1}));
ROILabels = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{6}).(condNames{1}).crossValidated); %possible, because S7 has all maps
ROILabels(contains(ROILabels,{'VO1','VO2','PHC'})) = []; %skip ventral stream ROIs 
ROIs = unique(erase(ROILabels, ["Right","Left","right","left"]),'stable');


save_dir_ecc = [save_path, 'ecc_comparisons/'];


%% make eccentricity plots -takes time (skip if unnecessary, takes a lot of time)
monoTuned_make_ecc_plots_anonymous(save_dir_ecc, cv_data, minVE, timing_maps,minDuration,maxDuration,minPeriod,maxPeriod,medianPerSubj);

%% ANOVAs
near_ecc = 1;
far_ecc = 2;
max_ecc = 5.6; % needs to be 5.6 if last 0.2 bin is centered around 5.5

%% get info in right format for analysis and plots
[meanSub, ANOVAstruc] = monoTuned_build_meanSub_anonymous(cv_data, minVE,1,near_ecc, far_ecc, max_ecc, medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod);

%% how to save the results
save_dir_anova = [save_dir_ecc, 'ecc_ANOVAs'];
cd(save_dir_anova)

save_name = 'meanSub';
save_name = strcat(save_name, '_', modelFieldNames{1}, '_', modelFieldNames{2});
if medianPerSubj == 1
    save_name = strcat(save_name, '_medianPerSubj');
elseif medianPerSubj == 0
    save_name = strcat(save_name, '_meanPerSubj');
elseif medianPerSubj == 2
    save_name = strcat(save_name, '_75thPerSubj');
elseif medianPerSubj == 3
    save_name = strcat(save_name, '_90thPerSubj');
end
if timing_maps == 1
    save_name = strcat(save_name, '_timing_maps');
elseif timing_maps == 0
    save_name = strcat(save_name, '_visual_field_maps');
end
save_name = strcat(save_name, '_minVE=', string(minVE));
save_name = strcat(save_name, '_nearEcc=', string(near_ecc), '_farEcc=', string(far_ecc));
save_name = strcat(save_name, '.mat');
save(save_name, 'meanSub')


%% do the ANOVAs and post-hocs
terms(1,1) = 1; % main subj
terms(2,2) = 1; % main map
terms(3,3) = 1; % main ecc
terms(4,2) = 1; terms(4,3) = 1; % interaction map x ecc

stat = [];

levels = {modelFieldNames{1},modelFieldNames{2},'diff'};
for level = 1:length(levels)
    veANOVA = []; 
    veANOVA = ANOVAstruc.VE.(levels{level});

    % ANOVAstruc had NaNs kicked out, but only for one eccentricity range
    missing_value_near_all = isnan(meanSub.(levels{level}).near(:));
    missing_value_far_all = isnan(meanSub.(levels{level}).far(:));
    
    % NOTE: too complicated programming, to remove values for which either near or far is missing
    % late addition, initially did anova on everything, but not anymore to keep anova & post hocs on the same data
    if sum(missing_value_near_all)>0 && sum(missing_value_far_all)==0
        missing_ids = find(missing_value_near_all);
        veANOVA(missing_ids + length(meanSub.(levels{level}).near(:)) - sum(missing_value_near_all)) = [];
        ANOVAstruc.VE.(levels{level}) = veANOVA;
        if level == 1 % shoud only be done once
        ANOVAstruc.subj(missing_ids + length(meanSub.(levels{level}).near(:)) - sum(missing_value_near_all)) = [];
        ANOVAstruc.roi(missing_ids + length(meanSub.(levels{level}).near(:)) - sum(missing_value_near_all)) = [];
        ANOVAstruc.ecc(missing_ids + length(meanSub.(levels{level}).near(:)) - sum(missing_value_near_all)) = [];
        end
    elseif sum(missing_value_far_all)>0 && sum(missing_value_near_all)==0
        missing_ids = find(missing_value_far_all);
        veANOVA(missing_ids) = [];
        ANOVAstruc.VE.(levels{level}) = veANOVA;
        if level == 1 % shoud only be done once
            ANOVAstruc.subj(missing_ids) = [];
            ANOVAstruc.roi(missing_ids) = [];
            ANOVAstruc.ecc(missing_ids) = [];
        end
    elseif sum(missing_value_far_all)>0 && sum(missing_value_near_all)>0
        only_missing_far = missing_value_far_all & ~missing_value_near_all;
        only_missing_near = missing_value_near_all & ~missing_value_far_all;
        only_missing_far(missing_value_near_all) = [];
        only_missing_near(missing_value_far_all) = [];
        only_missing_near = find(only_missing_near);
        
        veANOVA(only_missing_far) = [];
        if level == 1 % shoud only be done once
            ANOVAstruc.subj(only_missing_far) = [];
            ANOVAstruc.roi(only_missing_far) = [];
            ANOVAstruc.ecc(only_missing_far) = [];
        end
        
        veANOVA(only_missing_near + length(meanSub.(levels{level}).near(:)) - sum(missing_value_near_all) - sum(only_missing_far)) = [];
        ANOVAstruc.VE.(levels{level}) = veANOVA;
        if level == 1 % shoud only be done once
        ANOVAstruc.subj(only_missing_near + length(meanSub.(levels{level}).near(:)) - sum(missing_value_near_all)- sum(only_missing_far)) = [];
        ANOVAstruc.roi(only_missing_near + length(meanSub.(levels{level}).near(:)) - sum(missing_value_near_all)- sum(only_missing_far)) = [];
        ANOVAstruc.ecc(only_missing_near + length(meanSub.(levels{level}).near(:)) - sum(missing_value_near_all)- sum(only_missing_far)) = [];
        end
    
    end
    
    
    % ANOVA
    [stat.ANOVA.p.(levels{level}),stat.ANOVA.tbl.(levels{level}),stat.ANOVA.stats.(levels{level}),stat.ANOVA.terms.(levels{level})] = anovan(veANOVA,{ANOVAstruc.subj,ANOVAstruc.roi,ANOVAstruc.ecc},...
        'model',terms, 'varnames', {'subject', 'map', 'eccentricity'});
    
    % determine normality (incl FDR)
    pnorm_jb = [];
    pnorm_sw = [];
    for roi = 1:length(ROIs)
        % for paired comparisons: load data without NaNs
        % ALSO, if a near or far eccentricity is missing, remove the
        % corresponding value so paired comparisons are still possible
        missing_value_near = isnan(meanSub.(levels{level}).near(:,roi));
        missing_value_far = isnan(meanSub.(levels{level}).far(:,roi));
        missing_value_ids = missing_value_near | missing_value_far;

        stat.data.(levels{level}).(ROIs{roi}).near = meanSub.(levels{level}).near(~missing_value_ids,roi);
        stat.data.(levels{level}).(ROIs{roi}).far = meanSub.(levels{level}).far(~missing_value_ids,roi);
 
        if length(stat.data.(levels{level}).(ROIs{roi}).near) ~= length(stat.data.(levels{level}).(ROIs{roi}).far)
            fprintf('Near and far do not have same level of observations for %s, %s!!! PROBLEM\n', ROIs{roi}, levels{level});
        end
        
        if level == 1 % same for  all levels
            stat.data.(ROIs{roi}).N = length(stat.data.(levels{level}).(ROIs{roi}).near);
        end
        
        [stat.normality_jb.(levels{level}).(ROIs{roi}).h, stat.normality_jb.(levels{level}).(ROIs{roi}).p] = jbtest(stat.data.(levels{level}).(ROIs{roi}).near - stat.data.(levels{level}).(ROIs{roi}).far); 
        pnorm_jb = [pnorm_jb stat.normality_jb.(levels{level}).(ROIs{roi}).p];
        
        [stat.normality_sw.(levels{level}).(ROIs{roi}).h, stat.normality_sw.(levels{level}).(ROIs{roi}).p, stat.normality_sw.(levels{level}).(ROIs{roi}).w] = swtest(stat.data.(levels{level}).(ROIs{roi}).near - stat.data.(levels{level}).(ROIs{roi}).far);
        pnorm_sw = [pnorm_sw stat.normality_sw.(levels{level}).(ROIs{roi}).p];
    end
    [hnorm, crit_pnorm, adj_ci_cvrgnorm, stat.normality_jb.(levels{level}).adj_pnorm]=fdr_bh(pnorm_jb);
    [hnorm, crit_pnorm, adj_ci_cvrgnorm, stat.normality_sw.(levels{level}).adj_pnorm]=fdr_bh(pnorm_sw);
    
    % do appropriate paired post-hoc comparisons
    pvals = [];
    if sum(stat.normality_jb.(levels{level}).adj_pnorm<.05)>0
        for roi = 1:length(ROIs)
            [stat.wilcoxon.(levels{level}).(ROIs{roi}).p, stat.wilcoxon.(levels{level}).(ROIs{roi}).h, stat.wilcoxon.(levels{level}).(ROIs{roi}).stats] = signrank(stat.data.(levels{level}).(ROIs{roi}).near, stat.data.(levels{level}).(ROIs{roi}).far);
            pvals = [pvals stat.wilcoxon.(levels{level}).(ROIs{roi}).p];
        end
        
        % for plots
        stacked_plots_median(level) = 1;
    else
        for roi = 1:length(ROIs)
            [stat.t.(levels{level}).(ROIs{roi}).h, stat.t.(levels{level}).(ROIs{roi}).p, stat.t.(levels{level}).(ROIs{roi}).ci, stat.t.(levels{level}).(ROIs{roi}).stats] = ttest(stat.data.(levels{level}).(ROIs{roi}).near, stat.data.(levels{level}).(ROIs{roi}).far);
            pvals = [pvals stat.t.(levels{level}).(ROIs{roi}).p];
        end
        % for plots
        stacked_plots_median(level) = 0;
        
    end
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);
    
    
    % for plots
    signif.near_05.(levels{level}) = repelem(NaN, length(ROIs));
    signif.near_01.(levels{level}) = repelem(NaN, length(ROIs));
    signif.near_001.(levels{level}) = repelem(NaN, length(ROIs));
    signif.near_0001.(levels{level}) = repelem(NaN, length(ROIs));
    
    signif.far_05.(levels{level}) = repelem(NaN, length(ROIs));
    signif.far_01.(levels{level}) = repelem(NaN, length(ROIs));
    signif.far_001.(levels{level}) = repelem(NaN, length(ROIs));
    signif.far_0001.(levels{level}) = repelem(NaN, length(ROIs));
    
    if level == 3
        if medianPerSubj == 0 ||medianPerSubj == 1
            sig = 0.15;
        elseif medianPerSubj == 2 
            sig = 0.35;
        elseif medianPerSubj == 3
            sig = 0.55;
        end
    else
        if medianPerSubj == 0 ||medianPerSubj == 1
            sig = 0.3;
        elseif medianPerSubj == 2
            sig = 0.5;
        elseif medianPerSubj == 3
            sig = 0.7;
        end
    end
    
    % assign adjusted pvalues and get plot indications
    for roi = 1:length(ROIs)
        if sum(stat.normality_jb.(levels{level}).adj_pnorm<.05)>0
            stat.wilcoxon.(levels{level}).(ROIs{roi}).adj_p = adj_p(roi);
        else
            stat.t.(levels{level}).(ROIs{roi}).adj_p = adj_p(roi);
        end
        
        adjusted_p = adj_p(roi);
        if abs(mean(stat.data.(levels{level}).(ROIs{roi}).near)) > abs(mean(stat.data.(levels{level}).(ROIs{roi}).far))
            if adjusted_p<0.05
                signif.near_05.(levels{level})(roi) = sig - 0.005 * 0;
            end
            if adjusted_p<0.01
                signif.near_01.(levels{level})(roi) = sig - 0.005 * 1;
            end
            if adjusted_p<0.001
                signif.near_001.(levels{level})(roi) = sig - 0.005 * 2;
            end
            if adjusted_p<0.0001
                signif.near_0001.(levels{level})(roi) = sig - 0.005 * 3;
            end
        elseif abs(mean(stat.data.(levels{level}).(ROIs{roi}).near)) < abs(mean(stat.data.(levels{level}).(ROIs{roi}).far))
            if adjusted_p<0.05
                signif.far_05.(levels{level})(roi) = sig - 0.005 * 0;
            end
            if adjusted_p<0.01
                signif.far_01.(levels{level})(roi) = sig - 0.005 * 1;
            end
            if adjusted_p<0.001
                signif.far_001.(levels{level})(roi) = sig - 0.005 * 2;
            end
            if adjusted_p<0.05
                signif.far_0001.(levels{level})(roi) = sig - 0.005 * 3;
            end
        end
    end
    
end


%% plot stacked per eccentricity
save_dir_ecc_stacks = [save_dir_ecc, 'ecc_stacks'];
cd(save_dir_ecc_stacks)

%% calculate means and SEs or medians and CIs
 midpoints_mod1_near = repelem(NaN,length(ROIs));midpoints_mod1_far = repelem(NaN,length(ROIs));midpoints_mod2_near = repelem(NaN,length(ROIs));midpoints_mod2_far = repelem(NaN,length(ROIs));midpoints_diff_near = repelem(NaN,length(ROIs));midpoints_diff_far = repelem(NaN,length(ROIs));
if sum(stacked_plots_median) > 0   
    %use NaN removed data
    for roi =1:length(ROIs)
        midpoints_mod1_near(roi) =median(stat.data.(modelFieldNames{1}).(ROIs{roi}).near, 'omitnan');
        midpoints_mod2_near(roi) =median(stat.data.(modelFieldNames{2}).(ROIs{roi}).near, 'omitnan');
        midpoints_diff_near(roi) =median(stat.data.diff.(ROIs{roi}).near, 'omitnan');
        midpoints_mod1_far(roi) =median(stat.data.(modelFieldNames{1}).(ROIs{roi}).far, 'omitnan');
        midpoints_mod2_far(roi) =median(stat.data.(modelFieldNames{2}).(ROIs{roi}).far, 'omitnan');
        midpoints_diff_far(roi) =median(stat.data.diff.(ROIs{roi}).far, 'omitnan');
        
        barpoints_mod1_near(:,roi) = bootci(1000, {@median, stat.data.(modelFieldNames{1}).(ROIs{roi}).near},'alpha',0.05);
        stat.data.(modelFieldNames{1}).CI.(ROIs{roi}).near = barpoints_mod1_near(:,roi);
        barpoints_mod2_near(:,roi) = bootci(1000, {@median, stat.data.(modelFieldNames{2}).(ROIs{roi}).near},'alpha',0.05);
        stat.data.(modelFieldNames{2}).CI.(ROIs{roi}).near = barpoints_mod2_near(:,roi);
        barpoints_diff_near(:,roi) = bootci(1000, {@median, stat.data.diff.(ROIs{roi}).near},'alpha',0.05);
        stat.data.diff.CI.(ROIs{roi}).near = barpoints_diff_near(:,roi);
        barpoints_mod1_far(:,roi) = bootci(1000, {@median, stat.data.(modelFieldNames{1}).(ROIs{roi}).far},'alpha',0.05);
        stat.data.(modelFieldNames{1}).CI.(ROIs{roi}).far = barpoints_mod1_far(:,roi);
        barpoints_mod2_far(:,roi) = bootci(1000, {@median, stat.data.(modelFieldNames{2}).(ROIs{roi}).far},'alpha',0.05);
        stat.data.(modelFieldNames{2}).CI.(ROIs{roi}).far = barpoints_mod2_far(:,roi);
        barpoints_diff_far(:,roi) = bootci(1000, {@median, stat.data.diff.(ROIs{roi}).far},'alpha',0.05);
        stat.data.diff.CI.(ROIs{roi}).far = barpoints_diff_far(:,roi);

    end
    barpoints_mod1_near_low =  midpoints_mod1_near - barpoints_mod1_near(1,:);
    barpoints_mod1_near_high = barpoints_mod1_near(2,:) -  midpoints_mod1_near;
    barpoints_mod2_near_low =  midpoints_mod2_near - barpoints_mod2_near(1,:);
    barpoints_mod2_near_high = barpoints_mod2_near(2,:)-  midpoints_mod2_near;
    barpoints_diff_near_low = midpoints_diff_near - barpoints_diff_near(1,:);
    barpoints_diff_near_high = barpoints_diff_near(2,:)- midpoints_diff_near;
    barpoints_mod1_far_low = midpoints_mod1_far - barpoints_mod1_far(1,:);
    barpoints_mod1_far_high = barpoints_mod1_far(2,:) - midpoints_mod1_far;
    barpoints_mod2_far_low = midpoints_mod2_far - barpoints_mod2_far(1,:);
    barpoints_mod2_far_high = barpoints_mod2_far(2,:) - midpoints_mod2_far;
    barpoints_diff_far_low = midpoints_diff_far - barpoints_diff_far(1,:);
    barpoints_diff_far_high = barpoints_diff_far(2,:) - midpoints_diff_far;
    
else
    for roi = 1:length(ROIs)
        midpoints_mod1_near(roi) =nanmean(stat.data.(modelFieldNames{1}).(ROIs{roi}).near);
        midpoints_mod2_near(roi) =nanmean(stat.data.(modelFieldNames{2}).(ROIs{roi}).near);
        midpoints_diff_near(roi) =nanmean(stat.data.diff.(ROIs{roi}).near);
        midpoints_mod1_far(roi) =nanmean(stat.data.(modelFieldNames{1}).(ROIs{roi}).far);
        midpoints_mod2_far(roi) =nanmean(stat.data.(modelFieldNames{2}).(ROIs{roi}).far);
        midpoints_diff_far(roi) =nanmean(stat.data.diff.(ROIs{roi}).far);
        
        mod1_std_ve_near=nanstd(stat.data.(modelFieldNames{1}).(ROIs{roi}).near);
        mod2_std_ve_near=nanstd(stat.data.(modelFieldNames{2}).(ROIs{roi}).near);
        diff_std_ve_near=nanstd(stat.data.diff.(ROIs{roi}).near);
        mod1_std_ve_far=nanstd(stat.data.(modelFieldNames{1}).(ROIs{roi}).far);
        mod2_std_ve_far=nanstd(stat.data.(modelFieldNames{2}).(ROIs{roi}).far);
        diff_std_ve_far=nanstd(stat.data.diff.(ROIs{roi}).far);
        
        mod1_se_ve_near(roi)=mod1_std_ve_near/sqrt(stat.data.(ROIs{roi}).N);
        stat.data.(modelFieldNames{1}).SE.(ROIs{roi}).near = mod1_se_ve_near(roi);
        mod2_se_ve_near(roi)=mod2_std_ve_near/sqrt(stat.data.(ROIs{roi}).N);
        stat.data.(modelFieldNames{2}).SE.(ROIs{roi}).near = mod2_se_ve_near(roi);
        diff_se_ve_near(roi)=diff_std_ve_near/sqrt(stat.data.(ROIs{roi}).N);
        stat.data.diff.SE.(ROIs{roi}).near = diff_se_ve_near(roi);
        mod1_se_ve_far(roi)=mod1_std_ve_far/sqrt(stat.data.(ROIs{roi}).N);
        stat.data.(modelFieldNames{1}).SE.(ROIs{roi}).far = mod1_se_ve_far(roi);
        mod2_se_ve_far(roi)=mod2_std_ve_far/sqrt(stat.data.(ROIs{roi}).N);
        stat.data.(modelFieldNames{2}).SE.(ROIs{roi}).far = mod2_se_ve_far(roi);
        diff_se_ve_far(roi)=diff_std_ve_far/sqrt(stat.data.(ROIs{roi}).N);
        stat.data.diff.SE.(ROIs{roi}).far = diff_se_ve_far(roi);
        
    end
    barpoints_mod1_near_low = mod1_se_ve_near;
    barpoints_mod1_near_high = mod1_se_ve_near;
    barpoints_mod2_near_low = mod2_se_ve_near;
    barpoints_mod2_near_high = mod2_se_ve_near;
    barpoints_diff_near_low = diff_se_ve_near;
    barpoints_diff_near_high = diff_se_ve_near;
    barpoints_mod1_far_low = mod1_se_ve_far;
    barpoints_mod1_far_high = mod1_se_ve_far;
    barpoints_mod2_far_low = mod2_se_ve_far;
    barpoints_mod2_far_high = mod2_se_ve_far;
    barpoints_diff_far_low = diff_se_ve_far;
    barpoints_diff_far_high = diff_se_ve_far;
end


savename = 'stats_anova(cv)_bootciDefault';
savename = strcat(savename, '_', modelFieldNames{1}, '_', modelFieldNames{2}, '_');
if medianPerSubj == 1
    savename = strcat(savename, 'medianPerSubj');
elseif medianPerSubj == 0
    savename = strcat(savename, 'meanPerSubj');
elseif medianPerSubj == 2
    savename = strcat(savename, '75thPerSubj');
elseif medianPerSubj == 3
    savename = strcat(savename, '90thPerSubj');
end
if timing_maps == 0
    savename = strcat(savename, '_visual_field_maps');
elseif timing_maps == 1
    savename = strcat(savename, '_timing_maps');
end
savename = strcat(savename, '_minVE=', string(minVE));
savename = strcat(savename, '_nearEcc=', string(near_ecc), '_farEcc=', string(far_ecc));
savename = strcat(savename, '.mat');
save(savename, 'stat')
close all

%% plot stacked per model
if timing_maps ==1
    scatter_order=[4 10 9 8 2 6 5 7 1 3];
elseif timing_maps == 0
    scatter_order=[14:16,18,10:13,17,1:9];% %
end

x_LR=1:length(ROIs);
mod1_barNear=[.56,.56,1];mod1_barFar=[.8,.8,1];
mod2_barNear=[1,.56,.56];mod2_barFar=[1,.8,.8];
diff_barNear = [0.9290, 0.6940, 0.1250] ;diff_barFar = [0.99, 0.89, 0.4];

for level = 1:length(levels)
    fig = figure(5+level);
    if level == 1
        plot_minVe = 0;
        if medianPerSubj == 0 ||medianPerSubj == 1
            plot_maxVe = 0.3;
        elseif medianPerSubj == 2
            plot_maxVe = 0.5;
           elseif medianPerSubj == 3
            plot_maxVe = 0.7; 
        end
    elseif level == 2
        plot_minVe = 0;
        if medianPerSubj == 0 ||medianPerSubj == 1
            plot_maxVe = 0.3;
        elseif medianPerSubj == 2
            plot_maxVe = 0.5;
        elseif medianPerSubj == 3
            plot_maxVe = 0.7;
        end
    elseif level == 3
        plot_minVe = -0.15;
        if medianPerSubj == 0 ||medianPerSubj == 1
            plot_maxVe = 0.15;
        elseif medianPerSubj == 2
            plot_maxVe = 0.35;
        elseif medianPerSubj == 3
            plot_maxVe = 0.55;
        end
    end
    
    %% fill the plots
    hold on
    
    if level == 1
        mod1_far = errorbar(x_LR,midpoints_mod1_far(1,(scatter_order)),barpoints_mod1_far_low(1,(scatter_order)),barpoints_mod1_far_high(1,(scatter_order)), 'vertical','LineStyle', 'none','Color',mod1_barFar,'Marker','s','MarkerSize',5.3,'MarkerEdgeColor','k','MarkerFaceColor',mod1_barFar,'LineWidth',3,'Displayname', [modelFieldNames{1} '>2' char(176)]);
        drawnow;
        mod1_far.MarkerHandle.LineWidth = 0.05;
        
        mod1_near = errorbar(x_LR,midpoints_mod1_near(1,(scatter_order)),barpoints_mod1_near_low(1,(scatter_order)),barpoints_mod1_near_high(1,(scatter_order)),'vertical', 'LineStyle', 'none','Color',mod1_barNear,'Marker','s','MarkerSize',5.3,'MarkerEdgeColor','k','MarkerFaceColor',mod1_barNear,'LineWidth',3,'Displayname', [modelFieldNames{1} '<1' char(176)]);
        drawnow;
        mod1_near.MarkerHandle.LineWidth = 0.05;
        
        % redraw markers on top (otherwised covered by error bars
        scatter(x_LR,midpoints_mod1_far(1,scatter_order),24,'s','MarkerFaceColor',mod1_barFar,'MarkerEdgeColor','k','HandleVisibility', 'off');
        scatter(x_LR,midpoints_mod1_near(1,scatter_order),24,'s','MarkerFaceColor',mod1_barNear,'MarkerEdgeColor','k','HandleVisibility', 'off');
    elseif level == 2
        mod2_far = errorbar(x_LR,midpoints_mod2_far(1,(scatter_order)),barpoints_mod2_far_low(1,(scatter_order)),barpoints_mod2_far_high(1,(scatter_order)), 'vertical','LineStyle', 'none','Color',mod2_barFar,'Marker','d','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',mod2_barFar,'LineWidth',3,'Displayname', [modelFieldNames{2} '>2' char(176)]);
        drawnow;
        mod2_far.MarkerHandle.LineWidth = 0.05;
        
        mod2_near = errorbar(x_LR,midpoints_mod2_near(1,(scatter_order)),barpoints_mod2_near_low(1,(scatter_order)),barpoints_mod2_near_high(1,(scatter_order)),'vertical', 'LineStyle', 'none','Color',mod2_barNear,'Marker','d','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',mod2_barNear,'LineWidth',3,'Displayname',[modelFieldNames{2} '<1' char(176)]);
        drawnow;
        mod2_near.MarkerHandle.LineWidth = 0.05;
        
        % redraw markers on top (otherwised covered by error bars
        scatter(x_LR,midpoints_mod2_far(1,scatter_order),20,'d','MarkerFaceColor',mod2_barFar,'MarkerEdgeColor','k','HandleVisibility', 'off');
        scatter(x_LR,midpoints_mod2_near(1,scatter_order),20,'d','MarkerFaceColor',mod2_barNear,'MarkerEdgeColor','k','HandleVisibility', 'off');
        
    elseif level == 3
        line([0 length(ROIs)+1],[0 0],'Color','black','LineStyle','-','HandleVisibility', 'off');
        
        diff_far = errorbar(x_LR,midpoints_diff_far(1,(scatter_order)),barpoints_diff_far_low(1,(scatter_order)),barpoints_diff_far_high(1,(scatter_order)),'vertical', 'LineStyle', 'none','Color',diff_barFar,'Marker','o','MarkerSize',5.3,'MarkerEdgeColor','k','MarkerFaceColor',diff_barFar,'LineWidth',3,'Displayname', ['Difference >2' char(176)]);
        drawnow;
        diff_far.MarkerHandle.LineWidth = 0.05;
        
        diff_near = errorbar(x_LR,midpoints_diff_near(1,(scatter_order)),barpoints_diff_near_low(1,(scatter_order)),barpoints_diff_near_high(1,(scatter_order)),'vertical', 'LineStyle', 'none','Color',diff_barNear,'Marker','o','MarkerSize',5.3,'MarkerEdgeColor','k','MarkerFaceColor',diff_barNear,'LineWidth',3,'Displayname', ['Difference <1' char(176)]);
        drawnow;
        diff_near.MarkerHandle.LineWidth = 0.05;
        
        % redraw markers on top (otherwised covered by error bars
        scatter(x_LR,midpoints_diff_far(1,scatter_order),24,'o','MarkerFaceColor',diff_barFar,'MarkerEdgeColor','k','HandleVisibility', 'off');
        scatter(x_LR,midpoints_diff_near(1,scatter_order),24,'o','MarkerFaceColor',diff_barNear,'MarkerEdgeColor','k','HandleVisibility', 'off');
        
    end
    
    scatter(x_LR,signif.near_05.(levels{level})(1,(scatter_order)),20,'*','k','HandleVisibility', 'off');
    scatter(x_LR,signif.near_01.(levels{level})(1,(scatter_order)),20,'*','k','HandleVisibility', 'off');
    scatter(x_LR,signif.near_001.(levels{level})(1,(scatter_order)),20,'*','k','HandleVisibility', 'off');
    scatter(x_LR,signif.near_0001.(levels{level})(1,(scatter_order)),20,'*','k','HandleVisibility', 'off');
    
    scatter(x_LR,signif.far_05.(levels{level})(1,(scatter_order)),20,'*', 'b','HandleVisibility', 'off');
    scatter(x_LR,signif.far_01.(levels{level})(1,(scatter_order)),20,'*', 'b','HandleVisibility', 'off');
    scatter(x_LR,signif.far_001.(levels{level})(1,(scatter_order)),20,'*', 'b','HandleVisibility', 'off');
    scatter(x_LR,signif.far_0001.(levels{level})(1,(scatter_order)),20,'*', 'b','HandleVisibility', 'off');
    
    %% finish and save plots
    xticks(1:1:length(ROIs)+1);
    yticks(plot_minVe:0.05:plot_maxVe);
    xticklabels(ROIs(scatter_order(1:end)));
    xtickangle(90);box off
    axis([0,length(ROIs)+1,plot_minVe,plot_maxVe])
    

    set(gcf,'Color','w');        
    xlabel('Visual Field Maps','FontWeight','bold');
    if level == 1
        ylabel({'Variance Explained'},'FontWeight','bold');
         legend([mod1_near(1), mod1_far(1)],'Location','best')
        text(-1.7,.3,'A','FontWeight','bold')
    elseif level == 2
        ylabel({'Variance Explained'},'FontWeight','bold');
         legend([mod2_near(1), mod2_far(1)],'Location','best')
        text(-1.7,.3,'B','FontWeight','bold')
    elseif level == 3
        ylabel({'Difference in Variance Explained'},'FontWeight','bold');
         legend([diff_near(1), diff_far(1)],'Location','best')
        text(-2,.15,'C','FontWeight','bold')
    end
 
    set(gcf,'units','centimeters','position',[0.1 0.1 20 12]);    
    
    hold off
    
    savename = 'stacked_ecc(cv)_bootciDefault_';
    savename = strcat(savename, levels{level},'_');
    if medianPerSubj == 1
        savename = strcat(savename, 'medianPerSubj')
    elseif medianPerSubj == 0
        savename = strcat(savename, 'meanPerSubj')
    elseif medianPerSubj == 2
        savename = strcat(savename, '75thPerSubj')
    elseif medianPerSubj == 3
        savename = strcat(savename, '90thPerSubj')
    end
    if sum(stacked_plots_median) > 0
        savename = strcat(savename, '_median');
    else
        savename = strcat(savename, '_mean');
    end
    savename = strcat(savename, '_minVE=', string(minVE));
    savename = strcat(savename, '_nearEcc=', string(near_ecc), '_farEcc=', string(far_ecc));
    savefig(fig, char(strcat(savename,'.fig')))
    savename_eps = strcat(savename, '.eps');
    export_fig(savename_eps,'-eps','-r600','-painters');
    savename_png = strcat(savename, '.png');
    export_fig(savename_png,'-png','-r600','-painters');
    
    close all
end
end