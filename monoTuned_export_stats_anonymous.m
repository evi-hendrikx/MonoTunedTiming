function monoTuned_export_stats_anonymous(timing_maps, model_comparison, save_path)
%% export stats: writes structure with data to csv file so they can be used in JASP and reported in tables
% timing_maps: whether you want timing maps (1) or visual field maps (0)
% model_comparison: whether you want model comparison results (1) or eccentricity results (0)
% save_path: 

if timing_maps == 1
    scatter_order=[4 10 9 8 2 6 5 7 1 3];
elseif timing_maps == 0 && model_comparison == 1
    scatter_order=[14:16,18:21,10:13,17,1:9];
elseif timing_maps == 0 && model_comparison == 0
    scatter_order=[14:16,18,10:13,17,1:9];
end

%% model comparisons
if model_comparison == 1
    cd([save_path, 'model_comparisons'])
    if timing_maps == 1
        % timing maps model comparisons
        load('meanSub_MonoOcc_TunedLin2d_medianPerSubj_timing_maps_minVE=0.2.mat')
        load('stats_anova(cv)_bootciDefault_MonoOcc_TunedLin2d_medianPerSubj_timing_maps_minVE=0.2.mat')
    elseif timing_maps == 0
        load('meanSub_MonoOcc_TunedLin2d_medianPerSubj_visual_field_maps_minVE=0.2.mat')
        load('stats_anova(cv)_bootciDefault_MonoOcc_TunedLin2d_medianPerSubj_visual_field_maps_minVE=0.2.mat')
    end
    
    % load in like this so you know it is the same order as ve_data, meanSub, and stat
    fieldnames_mS = fieldnames(meanSub);
    modelFieldNames{1} = fieldnames_mS{1};
    modelFieldNames{2} = fieldnames_mS{2};
    
    ROIs = meanSub.roi(1,:); % could've taken any row, they're all the same
    
    %% JASP file
    model_tbl = table();
    model_tbl.pp = meanSub.subj(:,1); % could've taken any column, they're all the same
    for roi = 1:length(ROIs)
        for model =1:length(modelFieldNames)
            col_mod = [ROIs{roi},'_', modelFieldNames{model}];
            model_tbl.(col_mod) = meanSub.(modelFieldNames{model})(:,roi);
        end
    end
    
    if timing_maps == 1
        save_name = 'model_within_TM.csv'
    elseif timing_maps == 0
        save_name = 'model_within_VFM.csv';
    end
    writetable(model_tbl, save_name)
    
    %% summary stats
    stats_model_tbl = table();
    stats_model_tbl.ROI = ROIs(scatter_order)';
    mono_medians = median(meanSub.(modelFieldNames{1}),'omitnan');
    tuned_medians = median(meanSub.(modelFieldNames{2}),'omitnan');
    diff_medians = median(meanSub.diff,'omitnan');
    
    CI_mono = repelem("", length(ROIs)); CI_tuned = CI_mono; CI_diff = CI_mono; Z_val = CI_mono; eff_size = CI_mono; n = CI_mono;p_val=CI_mono;
    for roi = 1:length(ROIs)
        CI_mono(roi) = strcat(num2str(mono_medians(scatter_order(roi)),'%.2f'), " [", num2str(stat.data.(ROIs{scatter_order(roi)}).CI.(modelFieldNames{1})(1),'%.2f'), ", ", num2str(stat.data.(ROIs{scatter_order(roi)}).CI.(modelFieldNames{1})(2),'%.2f'), "]");
        CI_tuned(roi) = strcat(num2str(tuned_medians(scatter_order(roi)),'%.2f'), " [", num2str(stat.data.(ROIs{scatter_order(roi)}).CI.(modelFieldNames{2})(1),'%.2f'), ", ", num2str(stat.data.(ROIs{scatter_order(roi)}).CI.(modelFieldNames{2})(2),'%.2f'), "]");
        CI_diff(roi) = strcat(num2str(diff_medians(scatter_order(roi)),'%.2f'), " [", num2str(stat.data.(ROIs{scatter_order(roi)}).CI.diff(1),'%.2f'), ", ", num2str(stat.data.(ROIs{scatter_order(roi)}).CI.diff(2),'%.2f'), "]");

        n(roi) = num2str(stat.data.(ROIs{scatter_order(roi)}).N);
        try
        Z_val(roi) = num2str(stat.wilcoxon.(ROIs{scatter_order(roi)}).stats.zval,'%.2f');
        eff_size(roi) = num2str(stat.wilcoxon.(ROIs{scatter_order(roi)}).stats.zval/sqrt(stat.data.(ROIs{scatter_order(roi)}).N),'%.2f');
        p_val(roi) = num2str(stat.wilcoxon.(ROIs{scatter_order(roi)}).adj_p);
        catch
            Z_val(roi)="";
            eff_size(roi)="";
            p_val(roi) = "";
        end
    end
    
    stats_model_tbl.n = n';
    stats_model_tbl.monotonic = CI_mono';
    stats_model_tbl.tuned = CI_tuned';
    stats_model_tbl.difference = CI_diff';
    stats_model_tbl.Z_stat = Z_val';
    stats_model_tbl.Effect_size = eff_size';
    stats_model_tbl.p = p_val';
    
    save_name = 'stats_table_with_diff_';
    if timing_maps == 1
        save_name = strcat(save_name, 'timing_maps');
    elseif timing_maps == 0
        save_name = strcat(save_name, 'visual_field_maps');
    end
    save_name = strcat(save_name, '.csv');
    writetable(stats_model_tbl, save_name)
    
elseif model_comparison == 0
    %% ecc comparisons
    cd([save_path, 'ecc_comparisons/ecc_stacks'])
    load('stats_anova(cv)_bootciDefault_MonoOcc_TunedLin2d_meanPerSubj_visual_field_maps_minVE=0_nearEcc=1_farEcc=2.mat')
    cd([save_path, 'ecc_comparisons/ecc_ANOVAs'])
    load('meanSub_MonoOcc_TunedLin2d_meanPerSubj_visual_field_maps_minVE=0_nearEcc=1_farEcc=2.mat')

    % load in like this so you know it is the same order as ve_data, meanSub, and stat
    fieldnames_mS = fieldnames(meanSub);
    modelFieldNames{1} = fieldnames_mS{1};
    modelFieldNames{2} = fieldnames_mS{2};
    modelFieldNames{3} = fieldnames_mS{3};
    
    ROIs = meanSub.roi(1,:); % could've taken any row, they're all the same
    
    %% JASP file
    ecc_tbl = table();
    ecc_tbl.pp = meanSub.subj(:,1); % could've taken any column, they're all the same
    for roi = 1:length(ROIs)
        for model =1:length(modelFieldNames)
            % if one eccentricity range was missing, both were removed in the analyses. 
            % This only happens in one measurement of the iPCS
            near_values = meanSub.(modelFieldNames{model}).near(:,roi);
            far_values = meanSub.(modelFieldNames{model}).far(:,roi);
            far_values(isnan(near_values)) = NaN;
            near_values(isnan(far_values)) = NaN;
            col_near = [ROIs{roi},'_near_', modelFieldNames{model}];
            ecc_tbl.(col_near) = near_values;
            col_far = [ROIs{roi},'_far_',modelFieldNames{model}];
            ecc_tbl.(col_far) = far_values;
        end
    end
    
    writetable(ecc_tbl, 'ecc_within_VFM.csv')
    
    %% summary stats
    
    for model = 1:length(modelFieldNames)
        stats_ecc_tbl = table();
        stats_ecc_tbl.ROI = ROIs(scatter_order)';
 
        CI_near = repelem("", length(ROIs)); CI_far = CI_near; Z_val = CI_near; eff_size = CI_near; p_val = CI_near; n=CI_near;
        for roi = 1:length(ROIs)
            near_median = num2str(median(stat.data.(modelFieldNames{model}).(ROIs{scatter_order(roi)}).near),'%.2f');
            far_median = num2str(median(stat.data.(modelFieldNames{model}).(ROIs{scatter_order(roi)}).far),'%.2f');
            CI_near(roi) = strcat(num2str(near_median,'%.2f'), " [", num2str(stat.data.(modelFieldNames{model}).CI.(ROIs{scatter_order(roi)}).near(1),'%.2f'), ", ", num2str(stat.data.(modelFieldNames{model}).CI.(ROIs{scatter_order(roi)}).near(2),'%.2f'), "]");
            CI_far(roi) = strcat(num2str(far_median,'%.2f')," [", num2str(stat.data.(modelFieldNames{model}).CI.(ROIs{scatter_order(roi)}).far(1),'%.2f'), ", ", num2str(stat.data.(modelFieldNames{model}).CI.(ROIs{scatter_order(roi)}).far(2),'%.2f'), "]");
            
            Z_val(roi) = num2str(stat.wilcoxon.(modelFieldNames{model}).(ROIs{scatter_order(roi)}).stats.zval,'%.2f');
            n(roi) = num2str(stat.data.(ROIs{scatter_order(roi)}).N);
            eff_size(roi) = num2str(stat.wilcoxon.(modelFieldNames{model}).(ROIs{scatter_order(roi)}).stats.zval/sqrt(stat.data.(ROIs{scatter_order(roi)}).N),'%.2f');
            p_val(roi) = num2str(stat.wilcoxon.(modelFieldNames{model}).(ROIs{scatter_order(roi)}).adj_p);
        end
        
        stats_ecc_tbl.n = n';
        stats_ecc_tbl.near = CI_near';
        stats_ecc_tbl.far = CI_far';
        stats_ecc_tbl.Z_stat = Z_val';
        stats_ecc_tbl.Effect_size = eff_size';
        stats_ecc_tbl.p = p_val';
        
        
        save_name = 'stats_table_ecc_';
        save_name = strcat(save_name, modelFieldNames{model}, '_');
        if timing_maps == 1
            save_name = strcat(save_name, 'timing_maps');
        elseif timing_maps == 0
            save_name = strcat(save_name, 'visual_field_maps');
        end
        save_name = strcat(save_name, '.csv');
        writetable(stats_ecc_tbl, save_name)
        
    end
    
end
end