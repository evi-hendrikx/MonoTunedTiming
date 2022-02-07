function monoTuned_load_time_series_test_model_anonymous(paths, subj_id, save_path)
%% Loads model outputs on ground truth dataset into structs
% paths: where the information about each participant is stored 
% subj_id: index of subject you want (where in paths)
% save_path: general folder where you want results to be stored


%% Get the necessary data
datasets = {'Dmono', 'Dtuned'};
ROIsizeTuned=362700;
ROIsizeMono=147620;

cv_runs = {'Odd','Even'};
monoDTOdd=8;
monoDTEven=9;
tunedDTOdd=10:2:16;
tunedDTEven=11:2:17;

time_series = {};

scan_nr_timing = 4;

cd(paths{subj_id})
mrVista 3;

for dset = 1:length(datasets) % the actual designed dataset

    if dset == 1
        ROIsize = ROIsizeMono;
    elseif dset == 2
        ROIsize = ROIsizeTuned;
    end
    get_voxels = 1:ROIsize;
    
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
            for scan = 1:scan_nr_timing
                scan_name = ['Scan', char(num2str(scan))];
                cd(strcat(paths{subj_id}, '/Gray/', dataTYPES(DT_run(whichDT)).name,  '/TSeries/', scan_name))
                load('tSeries1');
                
                if whichDT == 1
                    time_series.(datasets{dset}).(cv_runs{cv}).(scan_name) = [];
                end
                time_series.(datasets{dset}).(cv_runs{cv}).(scan_name) = [time_series.(datasets{dset}).(cv_runs{cv}).(scan_name) tSeries(:,get_voxels)];
            end
        end

    end
end

cd(save_path);
save('time_series_validation.mat','time_series');

