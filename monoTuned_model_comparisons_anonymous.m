function monoTuned_model_comparisons_anonymous(save_path, minVE, timing_maps,medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod)
%% model_comparisons: compares the fits of the response models in various ROIs
% Compares the summary statistic of two models
% summary statistics are computed in monoTuned_build_meanSub_anonymous
% For a monotonic and tuned model: Runs ANOVA (subj, map, model,
% mapxmodel), post-hocs, and gives graphs for figure 3 and 7
%
% For two monotonic models: Runs wilcoxon (including all voxels in the rois), and gives graphs for these comparisons
% save_path: general folder where you want results to be stored
% minVE: threshold for inclusion of the voxels
% timing maps: 1 if analyze timing maps; 0 if analyze visual field maps
% medianPerSubj: determines the summary statistic: 1 = means; 2 = medians; 3 = 75th percentile; 4 = 90th percentile
% minDuration,maxDuration,minPeriod,maxPeriod: ranges outside which the tuned model is considered invalid

%% load and prepare data
if timing_maps == 1
    load('parameters_timing_maps.mat')
elseif timing_maps == 0
    load('parameters_visual_field_maps.mat')
end

%choose which models you want to use. I now use 2 models at a time, so the
%difference calculation is also easily achievable. You could of course
%compare multiple models at a time (but a couple of scripts, e.g., meanSub
%should then be altered)
modelFieldNames = fieldnames(ve_data);
use_models = {modelFieldNames{1},modelFieldNames{2}}; % fill in models you want, our main analyses are on MonoOcc and TunedLin2d

cv_data = monoTuned_cv_timing_data_anonymous(ve_data, use_models);
modelFieldNames = fieldnames(cv_data);
subjNames = fieldnames(cv_data.(modelFieldNames{1}));
condNames = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{1}));
ROILabels = fieldnames(cv_data.(modelFieldNames{1}).(subjNames{6}).(condNames{1}).crossValidated); %possible, because S7 has all maps
ROIs = unique(erase(ROILabels, ["Right","Left","right","left"]),'stable');

%% get info in right format for analyses and plots
[meanSub, ANOVAstruc] = monoTuned_build_meanSub_anonymous(cv_data, minVE, 0, [],[],[],medianPerSubj,minDuration,maxDuration,minPeriod,maxPeriod);

% save results
save_dir_model = [save_path, 'model_comparisons/'];
cd(save_dir_model)
save_name = 'meanSub';
save_name = strcat(save_name, '_', modelFieldNames{1}, '_', modelFieldNames{2}, '_');
if medianPerSubj == 1
    save_name = strcat(save_name, 'medianPerSubj');
elseif medianPerSubj == 0
    save_name = strcat(save_name, 'meanPerSubj');
elseif medianPerSubj == 2
    save_name = strcat(save_name, '75thPerSubj');
elseif medianPerSubj == 3
    save_name = strcat(save_name, '90thPerSubj');
end
if timing_maps == 1
    save_name = strcat(save_name, '_timing_maps');
elseif timing_maps == 0
    save_name = strcat(save_name, '_visual_field_maps');
end
save_name = strcat(save_name, '_minVE=', string(minVE));
save_name = strcat(save_name, '.mat');
save(save_name, 'meanSub')

%% do the ANOVAs and post-hocs for model comparisons
%% compare with tuned model
if find(strcmp(modelFieldNames,'TunedLin2d'))
    %% ANOVA
    terms(1,1) = 1; % main subj
    terms(2,2) = 1; % main map
    terms(3,3) = 1; % main model
    terms(4,2) = 1; terms(4,3) = 1; % interaction map x model
    
    stat = [];
    [stat.ANOVA.p,stat.ANOVA.tbl,stat.ANOVA.stats,stat.ANOVA.terms] = anovan(ANOVAstruc.VE,{ANOVAstruc.subj,ANOVAstruc.roi,ANOVAstruc.model},...
        'model',terms, 'varnames', {'subject', 'map', 'model'});
    
    %% normality
    % determine normality (incl FDR)
    % sanity check: eventually used jb results, but shapiro wilk would have resulted in the same choices
    pnorm_jb = [];
    pnorm_sw = [];
    for roi = 1:length(ROIs)
        % for paired comparisons: load data without NaNs per ROI
        stat.data.(ROIs{roi}).(modelFieldNames{1}) = meanSub.(modelFieldNames{1})(~isnan(meanSub.(modelFieldNames{1})(:,roi)),roi);
        stat.data.(ROIs{roi}).(modelFieldNames{2}) = meanSub.(modelFieldNames{2})(~isnan(meanSub.(modelFieldNames{2})(:,roi)),roi);
        stat.data.(ROIs{roi}).diff = meanSub.diff(~isnan(meanSub.diff(:,roi)),roi);
        stat.data.(ROIs{roi}).N = length(stat.data.(ROIs{roi}).(modelFieldNames{1}));
        
        try
            [stat.normality_jb.(ROIs{roi}).h, stat.normality_jb.(ROIs{roi}).p] = jbtest(stat.data.(ROIs{roi}).(modelFieldNames{2})-stat.data.(ROIs{roi}).(modelFieldNames{1}));
            pnorm_jb = [pnorm_jb stat.normality_jb.(ROIs{roi}).p];
        catch
            pnorm_jb = [pnorm_jb NaN];
        end
        
        try
            [stat.normality_sw.(ROIs{roi}).h, stat.normality_sw.(ROIs{roi}).p,stat.normality_sw.(ROIs{roi}).w] = swtest(stat.data.(ROIs{roi}).(modelFieldNames{2})-stat.data.(ROIs{roi}).(modelFieldNames{1}));
            pnorm_sw = [pnorm_sw stat.normality_sw.(ROIs{roi}).p];
        catch
            pnorm_sw = [pnorm_sw NaN];
        end
    end
    [hnorm, crit_pnorm, adj_ci_cvrgnorm, stat.normality_jb.adj_pnorm]=fdr_bh(pnorm_jb);
    [hnorm, crit_pnorm, adj_ci_cvrgnorm, stat.normality_sw.adj_pnorm]=fdr_bh(pnorm_sw);
    
    
    %% paired post-hoc comparisons
    % HOWEVER: if either timing maps or visual field maps are not normally
    % distributed, stick with non-parametric test for symmetry
    
    % %     ??? Run non parametric tests, because taking the mean of the medians per subject is just going to be super confusing ???
    pvals = [];
    if sum(stat.normality_jb.adj_pnorm<.05)>0
        for roi = 1:length(ROIs)
            if length(stat.data.(ROIs{roi}).(modelFieldNames{2})) > 1
                [stat.wilcoxon.(ROIs{roi}).p, stat.wilcoxon.(ROIs{roi}).h, stat.wilcoxon.(ROIs{roi}).stats] = signrank(stat.data.(ROIs{roi}).(modelFieldNames{1}), stat.data.(ROIs{roi}).(modelFieldNames{2}), 'method','approximate');
                pvals = [pvals stat.wilcoxon.(ROIs{roi}).p];
            else
                stat.wilcoxon.(ROIs{roi}).p = NaN;
                stat.wilcoxon.(ROIs{roi}).stats = NaN;
                pvals = [pvals NaN];
            end
        end
        
        % for plots
        stacked_plots_median = 1;
    else
        for roi = 1:length(ROIs)  
            % ttest in matlab is paired if you give in two samples
            % already gives NaN if there is 2 vectors to be compared with each only 1 element
            [stat.t.(ROIs{roi}).h, stat.t.(ROIs{roi}).p, stat.t.(ROIs{roi}).ci, stat.t.(ROIs{roi}).stats] = ttest(stat.data.(ROIs{roi}).(modelFieldNames{1}), stat.data.(ROIs{roi}).(modelFieldNames{2}));
            pvals = [pvals stat.t.(ROIs{roi}).p];
        end
        
        % for plots
        stacked_plots_median = 0;
    end
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals);
    
    %% plots: significance markers
    signif_001 = repelem(NaN, length(ROIs));
    signif_01 = repelem(NaN, length(ROIs));
    signif_05 = repelem(NaN, length(ROIs));
    sig = 0.55;
    
    % assign adjusted pvalues and get plot indications
    for roi = 1:length(ROIs)
        if stacked_plots_median == 1
            stat.wilcoxon.(ROIs{roi}).adj_p = adj_p(roi);
        else
            stat.t.(ROIs{roi}).adj_p = adj_p(roi);
        end
        adjusted_p = adj_p(roi);
        
        if adjusted_p<0.001
            signif_001(roi) = sig - 0.01 * 2;
        end
        if adjusted_p<0.01
            signif_01(roi) = sig - 0.01 * 1;
        end
        if adjusted_p<0.05
            signif_05(roi) = sig - 0.01 * 0;
        end
        
    end
    
    
    %% plots: values markers and bars
    % mean & SE values per roi
    if stacked_plots_median > 0
        midpoints_mod1 =median(meanSub.(modelFieldNames{1}), 'omitnan');
        midpoints_mod2 =median(meanSub.(modelFieldNames{2}), 'omitnan');
        midpoints_diff =median(meanSub.diff, 'omitnan');
        
        % use NaN removed data
        for roi = 1:length(ROIs)
            try
            % save CIs, because they're bootstrapped
            barpoints_mod1(:,roi) = bootci(1000, {@median, stat.data.(ROIs{roi}).(modelFieldNames{1})},'alpha',0.05);
            stat.data.(ROIs{roi}).CI.(modelFieldNames{1}) = barpoints_mod1(:,roi);
            barpoints_mod2(:,roi) = bootci(1000, {@median, stat.data.(ROIs{roi}).(modelFieldNames{2})},'alpha',0.05);
            stat.data.(ROIs{roi}).CI.(modelFieldNames{2}) = barpoints_mod2(:,roi);
            barpoints_diff(:,roi) = bootci(1000, {@median, stat.data.(ROIs{roi}).diff},'alpha',0.05);
            stat.data.(ROIs{roi}).CI.diff = barpoints_diff(:,roi);
            catch
                barpoints_mod1(1,roi) = NaN;barpoints_mod1(2,roi) = NaN; barpoints_mod2(1,roi) = NaN;barpoints_mod2(2,roi) = NaN; barpoints_diff(1,roi) = NaN;barpoints_diff(2,roi) = NaN;
                stat.data.(ROIs{roi}).CI.(modelFieldNames{1}) = barpoints_mod1(:,roi);
                stat.data.(ROIs{roi}).CI.(modelFieldNames{2}) = barpoints_mod2(:,roi);
                stat.data.(ROIs{roi}).CI.diff = barpoints_diff(:,roi);
            end
        end
        barpoints_mod1_low =  midpoints_mod1 - barpoints_mod1(1,:);
        barpoints_mod1_high = barpoints_mod1(2,:) -  midpoints_mod1;
        barpoints_mod2_low =  midpoints_mod2 - barpoints_mod2(1,:);
        barpoints_mod2_high = barpoints_mod2(2,:)-  midpoints_mod2;
        barpoints_diff_low = midpoints_diff - barpoints_diff(1,:);
        barpoints_diff_high = barpoints_diff(2,:)- midpoints_diff;
    else
        midpoints_mod1 =nanmean(meanSub.(modelFieldNames{1}));
        midpoints_mod2 =nanmean(meanSub.(modelFieldNames{2}));
        midpoints_diff =nanmean(meanSub.diff);
        
        mod1_std_ve=nanstd(meanSub.(modelFieldNames{1}));
        mod2_std_ve=nanstd(meanSub.(modelFieldNames{2}));
        diff_std_ve=nanstd(meanSub.diff);
        
        for roi = 1:length(ROIs)
            mod1_se_ve(roi) = (mod1_std_ve(roi)/sqrt(stat.data.(ROIs{roi}).N));
            mod2_se_ve(roi)= (mod2_std_ve(roi)/sqrt(stat.data.(ROIs{roi}).N));
            diff_se_ve(roi) = (diff_std_ve(roi)/sqrt(stat.data.(ROIs{roi}).N));
        end
        barpoints_mod1_low = mod1_se_ve;
        barpoints_mod1_high = mod1_se_ve;
        barpoints_mod2_low = mod2_se_ve;
        barpoints_mod2_high = mod2_se_ve;
        barpoints_diff_low = diff_se_ve;
        barpoints_diff_high = diff_se_ve;
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
    savename = strcat(savename, '.mat');
    save(savename, 'stat');
    
    
    %% plots: plot it!
    % set axis and orders
    close all
    if timing_maps ==1
        scatter_order=[4 10 9 8 2 6 5 7 1 3];
    elseif timing_maps == 0
        scatter_order=[14:16,18:21,10:13,17,1:9];
    end
    
    x_LR=1:length(ROIs);
    mod1_bar=[.56,.56,1];
    mod2_bar=[1,.56,.56];
    diff_bar = [0.9290, 0.6940, 0.1250];
    
    plots = {'models', 'difference'};
    legend_names = modelFieldNames;
    
    for plot_kind = 1:length(plots)
        figure(plot_kind)
        hold on
        if plot_kind == 1
            plot_maxVe = 0.6;
            plot_minVe = 0;
        elseif plot_kind == 2
            plot_maxVe = 0.35;
            plot_minVe = -0.25;
            signif_001 = signif_001 - 0.25;
            signif_01 = signif_01 - 0.25;
            signif_05 = signif_05 - 0.25;
        end
        
        % create bars
        if plot_kind == 1
            mod1_err_bar = errorbar(x_LR,midpoints_mod1(1,(scatter_order)),barpoints_mod1_low(1,(scatter_order)),barpoints_mod1_high(1,(scatter_order)), 'vertical','LineStyle', 'none','Color',mod1_bar,'Marker','s','MarkerSize',5.3,'MarkerEdgeColor','k','MarkerFaceColor',mod1_bar,'LineWidth',3,'Displayname', legend_names{1});
            hold on
            drawnow;
            mod1_err_bar.MarkerHandle.LineWidth = 0.05;
            
            %mod2
            mod2_err_bar = errorbar(x_LR,midpoints_mod2(1,scatter_order),barpoints_mod2_low(1,(scatter_order)),barpoints_mod2_high(1,(scatter_order)),'vertical', 'LineStyle', 'none', 'Color',mod2_bar,'LineWidth',3,'Marker','d','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',mod2_bar,'Displayname', legend_names{2});
            drawnow;
            mod2_err_bar.MarkerHandle.LineWidth = 0.05;
            
        elseif plot_kind == 2
            line([0 length(ROIs)+1],[0 0],'Color','black','LineStyle','-','LineWidth',.5)
            hold on
            diff_err_bar = errorbar(x_LR,midpoints_diff(1,scatter_order),barpoints_diff_low(1,(scatter_order)),barpoints_diff_high(1,(scatter_order)),'vertical', 'LineStyle', 'none', 'Color',diff_bar,'LineWidth',3, 'Marker','o','MarkerSize',5.3,'MarkerEdgeColor','k','MarkerFaceColor',diff_bar,'Displayname', 'Difference');
            drawnow;
            diff_err_bar.MarkerHandle.LineWidth = 0.05;
            
        end
        
        % create significance markers
        scatter(x_LR,signif_001(1,(scatter_order)),24,'*','k', 'HandleVisibility', 'off');
        scatter(x_LR,signif_01(1,(scatter_order)),24,'*','k', 'HandleVisibility', 'off');
        scatter(x_LR,signif_05(1,(scatter_order)),24,'*','k', 'HandleVisibility', 'off');
        
        % finish and save plots
        if plot_kind == 1
            legend([mod1_err_bar(1), mod2_err_bar(1)],'Location', 'southeast');
        elseif plot_kind == 2
            legend(diff_err_bar(1),'Location', 'southeast');
        end
        
        set(gcf,'Color','w');
        if timing_maps ==0
            xlabel('Visual Field Maps','FontWeight','bold');
        elseif timing_maps == 1
            xlabel('Timing Maps','FontWeight','bold');
        end
        if plot_kind == 1
            ylabel({'Variance Explained'},'FontWeight','bold');
        elseif plot_kind == 2
            ylabel({'Difference in Variance Explained'},'FontWeight','bold');
        end
        
        set(gcf,'units','centimeters','position',[0.1 0.1 20 12]);
        if timing_maps == 0
            x_pos = -1.4;
        elseif timing_maps == 1
            x_pos = -1;
        end
        if plot_kind == 1
            text(x_pos,.55,'A','FontWeight','bold')
        elseif plot_kind == 2
            text(x_pos,.3,'B','FontWeight','bold')
        end
        
        xticks(1:1:length(ROIs)+1);
        xticklabels([ROIs(scatter_order(1:end))]);
        xtickangle(90);box off
        
        yticks(plot_minVe:0.05:plot_maxVe);
        yticklabels(plot_minVe:0.05:plot_maxVe);
        
        
        axis([0,length(ROIs)+1,plot_minVe,plot_maxVe])
        
        
        savename = 'stacked_models(cv)_bootciDefault_';
        savename = strcat(savename, modelFieldNames{1}, '_', modelFieldNames{2}, '_');
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
        if plot_kind == 1
            savename = strcat(savename, '_onlymodels');
        elseif plot_kind == 2
            savename = strcat(savename, '_onlydifference');
        end
        if timing_maps == 1
            savename = strcat(savename, '_timing_maps');
        elseif timing_maps == 0
            savename = strcat(savename, '_visual_field_maps');
        end
        savename = strcat(savename, '_minVE=', string(minVE));
        savename_eps = strcat(savename, '.eps');
        %         savename_png = strcat(savename, '.png');
        export_fig(savename_eps,'-eps','-r600','-painters');
        %         export_fig(savename_png,'-png','-r600','-painters');
        disp(savename)
        close all
    end
    
    %% compare two monotonic models
else
    
    % for paired comparisons: load data without NaNs
    meanSub_all_1 = meanSub.(modelFieldNames{1})(:);
    stat.data.(modelFieldNames{1}) = meanSub_all_1(~isnan(meanSub_all_1));
    meanSub_all_2 = meanSub.(modelFieldNames{2})(:);
    stat.data.(modelFieldNames{2}) = meanSub_all_2(~isnan(meanSub_all_2));
    meanSub_all_diff = meanSub.diff(:);
    stat.data.diff = meanSub_all_diff(~isnan(meanSub_all_diff));
    stat.data.N = length(stat.data.(modelFieldNames{1}));
    
    %% normality
    [stat.normality_jb.h, stat.normality_jb.p] = jbtest(stat.data.(modelFieldNames{2})-stat.data.(modelFieldNames{1}));
    [stat.normality_sw.h, stat.normality_sw.p, stat.normality_sw.w] = swtest(stat.data.(modelFieldNames{2})-stat.data.(modelFieldNames{1}));
    
    %% paired comparisons
    if stat.normality_jb.p<.05
        [stat.wilcoxon.p, stat.wilcoxon.h, stat.wilcoxon.stats] = signrank(stat.data.(modelFieldNames{1}), stat.data.(modelFieldNames{2}));
        % for plots
        stacked_plots_median = 1;
    else
        % ttest in matlab is paired if you give in two samples
        [stat.t.(ROIs{roi}).h, stat.t.(ROIs{roi}).p, stat.t.(ROIs{roi}).ci, stat.t.(ROIs{roi}).stats] = ttest(stat.data.(ROIs{roi}).(modelFieldNames{1}), stat.data.(ROIs{roi}).(modelFieldNames{2}));
        % for plots
        stacked_plots_median = 0;
    end
    
% %     %% plots: significance markers
% %     signif_001 = NaN; signif_01 = NaN; signif_05 = NaN;
% %     sig = 0.32;
% %     
% %     if stat.normality_jb.p<.05
% %         p_value = stat.wilcoxon.p;
% %     else
% %         p_value = stat.t.p;
% %     end
% %     
% %     if p_value<0.001
% %         signif_001 = sig - 0.01 * 2;
% %     end
% %     if p_value<0.01
% %         signif_01 = sig - 0.01 * 1;
% %     end
% %     if p_value<0.05
% %         signif_05 = sig - 0.01 * 0;
% %     end
% %     
% %     
% %     %% plots: values markers and bars
% %     if stacked_plots_median > 0
% %         midpoints_mod1 =median(stat.data.(modelFieldNames{1}));
% %         midpoints_mod2 =median(stat.data.(modelFieldNames{2}));
% %         midpoints_diff =median(stat.data.diff);
% %         
% %         barpoints_mod1 = bootci(1000, {@median, stat.data.(modelFieldNames{1})},'alpha',0.05);
% %         stat.data.CI.(modelFieldNames{1}) = barpoints_mod1;
% %         barpoints_mod2 = bootci(1000, {@median, stat.data.(modelFieldNames{2})},'alpha',0.05);
% %         stat.data.CI.(modelFieldNames{2}) = barpoints_mod2;
% %         barpoints_diff = bootci(1000, {@median, stat.data.diff},'alpha',0.05);
% %         stat.data.CI.diff = barpoints_diff;
% %         
% %         barpoints_mod1_low =  midpoints_mod1 - barpoints_mod1(1,:);
% %         barpoints_mod1_high = barpoints_mod1(2,:) -  midpoints_mod1;
% %         barpoints_mod2_low =  midpoints_mod2 - barpoints_mod2(1,:);
% %         barpoints_mod2_high = barpoints_mod2(2,:)-  midpoints_mod2;
% %         barpoints_diff_low = midpoints_diff - barpoints_diff(1,:);
% %         barpoints_diff_high = barpoints_diff(2,:)- midpoints_diff;
% %     else
% %         midpoints_mod1 =nanmean(stat.data.(modelFieldNames{1}));
% %         midpoints_mod2 =nanmean(stat.data.(modelFieldNames{2}));
% %         midpoints_diff =nanmean(stat.data.diff);
% %         
% %         mod1_std_ve=nanstd(stat.data.(modelFieldNames{1}));
% %         mod2_std_ve=nanstd(stat.data.(modelFieldNames{2}));
% %         diff_std_ve=nanstd(stat.data.diff);
% %         
% %         
% %         mod1_se_ve = (mod1_std_ve/sqrt(stat.data.N));
% %         mod2_se_ve= (mod2_std_ve/sqrt(stat.data.N));
% %         diff_se_ve = (diff_std_ve/sqrt(stat.data.N));
% %         
% %         barpoints_mod1_low = mod1_se_ve;
% %         barpoints_mod1_high = mod1_se_ve;
% %         barpoints_mod2_low = mod2_se_ve;
% %         barpoints_mod2_high = mod2_se_ve;
% %         barpoints_diff_low = diff_se_ve;
% %         barpoints_diff_high = diff_se_ve;
% %     end
% %     
% %     
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
    savename = strcat(savename, '.mat');
    save(savename, 'stat')
% %     
% %     
% %     %% plots: plot it!
% %     % set axis and orders
% %     x_LR=1:length(modelFieldNames);
% %     figure(99);
% %     mod1_bar=[.56,.56,1];
% %     mod2_bar=[.56,.56,.56];
% %     
% %     legend_names = modelFieldNames;
% %     plot_maxVe = 0.3;
% %     plot_minVe = 0;
% %     
% %     mod1_err_bar = errorbar(1,midpoints_mod1,barpoints_mod1_low,barpoints_mod1_high, 'vertical','LineStyle', 'none','Color',mod1_bar,'Marker','s','MarkerSize',5.3,'MarkerEdgeColor','k','MarkerFaceColor',mod1_bar,'LineWidth',3,'Displayname', legend_names{1});
% %     hold on
% %     drawnow;
% %     mod1_err_bar.MarkerHandle.LineWidth = 0.05;
% %     
% %     %mod2
% %     mod2_err_bar = errorbar(2,midpoints_mod2,barpoints_mod2_low,barpoints_mod2_high,'vertical', 'LineStyle', 'none', 'Color',mod2_bar,'LineWidth',3,'Marker','d','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',mod2_bar,'Displayname', legend_names{2});
% %     drawnow;
% %     mod2_err_bar.MarkerHandle.LineWidth = 0.05;
% %     
% %     % create significance markers
% %     scatter(1.5,signif_001,24,'*','k', 'HandleVisibility', 'off');
% %     scatter(1.5,signif_01,24,'*','k', 'HandleVisibility', 'off');
% %     scatter(1.5,signif_05,24,'*','k', 'HandleVisibility', 'off');
% %     
% %     % finish and save plots
% %     legend([mod1_err_bar(1), mod2_err_bar(1)],'Location', 'southeast');
% %     
% %     
% %     set(gcf,'Color','w');
% %     xlabel('Neural response models','FontWeight','bold');
% %     ylabel({'Variance Explained'},'FontWeight','bold');
% %     set(gcf,'units','centimeters','position',[0.1 0.1 20 12]);
% %     
% %     if timing_maps == 0
% %         x_pos = -1.4;
% %     elseif timing_maps == 1
% %         x_pos = -1;
% %     end
% %     
% %     text(x_pos,.3,'A','FontWeight','bold')
% %     
% %     xticks(1:1:length(modelFieldNames)+1);
% %     xticklabels([modelFieldNames]);
% %     xtickangle(0);box off
% %     
% %     yticks(plot_minVe:0.05:plot_maxVe);
% %     yticklabels(plot_minVe:0.05:plot_maxVe);
% %     
% %     axis([0,length(modelFieldNames)+1,plot_minVe,plot_maxVe])
% %     
% %     savename = 'stacked_models(cv)_bootciDefault_';
% %     savename = strcat(savename, modelFieldNames{1}, '_', modelFieldNames{2}, '_');
% %     if medianPerSubj == 1
% %         savename = strcat(savename, 'medianPerSubj');
% %     elseif medianPerSubj == 0
% %         savename = strcat(savename, 'meanPerSubj');
% %     elseif medianPerSubj == 2
% %         savename = strcat(savename, '75thPerSubj');
% %     elseif medianPerSubj == 3
% %         savename = strcat(savename, '90thPerSubj');
% %     end
% %     if sum(stacked_plots_median) > 0
% %         savename = strcat(savename, '_median');
% %     else
% %         savename = strcat(savename, '_mean');
% %     end
% %     
% %     savename = strcat(savename, '_onlymodels');
% %     
% %     if timing_maps == 1
% %         savename = strcat(savename, '_timing_maps');
% %     elseif timing_maps == 0
% %         savename = strcat(savename, '_visual_field_maps');
% %     end
% %     savename = strcat(savename, '_minVE=', string(minVE));
% %     savename_eps = strcat(savename, '.eps');
% %     export_fig(savename_eps,'-eps','-r600','-painters');
% %     disp(savename)
% %     close all
end


end
