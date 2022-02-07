function fitting_figures_anonymous(paths, subj_id, whichVoxel, datapoint, roi_name, cv_run)
%% fitting figures: makes elements of Figure 1 (model fitting) and Supplementary Figure 1. 
% paths: where the information about each participant is stored 
% subj_id: index of subject you want (where in paths)
% whichVoxel: voxel coordinate (like in iCoords of ve_data)
% datapoint: index of voxel coordinates within the ROI
% roi_name: filename for roi in which the voxel is located
% cv_run: for which run (TimingSweeps, EvenScans, OddScans) you want the info


%% info on monotonic model for selected voxel
cd(fullfile(strcat(paths{subj_id},'/Gray/', cv_run, ' (L1)/Monotonic/')));
if string(cv_run) ~= "TimingSweeps"
    cd('xvalRefit');
end
model_name = '*202109*Lin-compressiveXcompressiveYNoNormOccupancy-DurFreq-20-gFit*';
load(strcat(ls(model_name)))

%% compressive exponents
%Plot amplitude as a function of frequency
frequencies=0.01:0.01:20;
response=frequencies.^model{1}.y0(whichVoxel);
%response=frequencies.^1; %to check everything works, use linear function
figure; plot(frequencies, response)
axis square;
axis([0 20 0 20]);
xlabel('Frequency (Hz)')
ylabel('Response amplitude')

%Plot amplitude per event as a function of frequency
figure; plot(frequencies, response./frequencies)
axis square;
axis([0 20 0 1]);
xlabel('Frequency (Hz)')
ylabel('Response amplitude per event')

%Plot amplitude per event as a function of duration
durations=0.05:0.05:1;
responseDuration=durations.^model{1}.x0(whichVoxel);
figure; plot(durations, responseDuration)
axis square;
axis([0 1 0 1]);
xlabel('Duration (s)')
ylabel('Response amplitude per event')


%% calculate ratio duration and frequency component
% get betas
grayROI=load(fullfile(strcat(paths{subj_id},'/Gray/ROIs/gray-Layer1.mat')));
cd(paths{subj_id})

VOLUME{1} = initHiddenGray;
grayCoords=grayROI.ROI.coords;
[tmp, betaCrds] = intersectCols(grayCoords, VOLUME{1}.coords(:,whichVoxel));
betas = squeeze(model{1}.beta(1,betaCrds,:));
beta1=betas(1);
beta2=betas(2);
close all

%% 2D plane images combined
modelPeriods = [0:0.005:1.05 2.055:0.005:2.1];
modelDurations = [0:0.005:1.05 1.955:0.005:2];
modelFrequencies = 1./modelPeriods'; % do I have to make this 5 like tuned??

amplitudeDuration = (modelDurations.^model{1}.x0(whichVoxel)).*ones(size(modelPeriods))';
amplitudeFrequency = ones(size(modelDurations)).*((modelFrequencies.^model{1}.y0(whichVoxel))./modelFrequencies);
amplitudeMono = amplitudeDuration * betas(1) + amplitudeFrequency * betas(2);
figure; imagesc(flipud(amplitudeMono)); axis image
max_amp_mono = caxis;

%% conditions info
% constant luminanace
durations1=[0.05:0.05:1 repmat(2, [1,6])];
periods1=[0.05:0.05:1 repmat(2.1, [1,6])];
% constant duration
durations2=repmat(0.05, [1,26]);
periods2=[0.05:0.05:1, repmat(2.1, [1,6])];
% constant period
durations3=[0.05:0.05:1 repmat(2, [1,6])];
periods3=[repmat(1, [1,20]), repmat(2.1, [1,6])];
% gaps
durations4=[0.05:0.05:0.5, repmat(0.05, [1,3]), 0.05:0.05:0.5, repmat(0.05, [1,3])];
periods4=[0.95:-0.05:0.5, repmat(2.1, [1,3]),0.55:0.05:1,repmat(2.1, [1,3])];

frequencies1=1./periods1;
frequencies2=1./periods2;
frequencies3=1./periods3;
frequencies4=1./periods4;


%% predictions of neural responses per condition combined
responsesMono1=(durations1.^model{1}.x0(whichVoxel)*betas(1)+((frequencies1.^model{1}.y0(whichVoxel))./frequencies1)*betas(2))/max_amp_mono(2);
responsesMono2=(durations2.^model{1}.x0(whichVoxel)*betas(1)+((frequencies2.^model{1}.y0(whichVoxel))./frequencies2)*betas(2))/max_amp_mono(2);
responsesMono3=(durations3.^model{1}.x0(whichVoxel)*betas(1)+((frequencies3.^model{1}.y0(whichVoxel))./frequencies3)*betas(2))/max_amp_mono(2);
responsesMono4=(durations4.^model{1}.x0(whichVoxel)*betas(1)+((frequencies4.^model{1}.y0(whichVoxel))./frequencies4)*betas(2))/max_amp_mono(2);

%Plotting this is more complex
colStr{1}='k';
colStr{2}='b';
colStr{3}='r';
colStr{4}='g';

RF=[responsesMono1; responsesMono1];
RF=RF(:);
Xaxis=[1:26;1:26];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{1})
xticks(0.5:2:25.5)
xticklabels(durations1(1:2:end))
xlabel('Duration (s)')
ylabel('Response amplitude per event')
axis square
xlim([0,26])
ylim([0,1])

RF=[responsesMono2; responsesMono2];
RF=RF(:);
Xaxis=[1:26;1:26];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{2})
xticks(0.5:2:25.5)
xticklabels(durations2(1:2:end))
xlabel('Duration (s)')
ylabel('Response amplitude per event')
axis square;
ylim([0 1]);
xlim([0,26]);

RF=[responsesMono3; responsesMono3];
RF=RF(:);
Xaxis=[1:26;1:26];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{3})
xticks(0.5:2:25.5)
xticklabels(durations3(1:2:end))
xlabel('Duration (s)')
ylabel('Response amplitude per event')
axis square;
ylim([0 1]);
xlim([0,26]);

RF=[responsesMono4; responsesMono4];
RF=RF(:);
Xaxis=[1:26;1:26];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{4})
xticks(0.5:2:25.5)
xticklabels(durations4(1:2:end))
xlabel('Duration (s)')
ylabel('Response amplitude per event')
axis square;
ylim([0 1]);
xlim([0,26]);
close all

%% get actual voxel data
% predicted fMRI over time
combinedDT=5:7;
dt = ["TimingSweeps","OddScans","EvenScans"]==string(cv_run);
cd(paths{subj_id});
mrVista 3
VOLUME{1} = initHiddenGray;

%IMPORTANT!!! Load the model you want, then put a breakpoint in rmGridFit line 170 %% IMPORTANT!!!
rmRunDurFreq2d(VOLUME{1},combinedDT(dt),roi_name,13,{'1g'},[],[],[],0);
assignin('base','params',params)
assignin('base','allstimimages',allstimimages)
% CONTINUE to rmGridFit line 363
assignin('base','data',data) % voxel data
dbquit
voxelData=data(:,datapoint);

%Generate prediction for duration response (monotonic)
rf=params.analysis.X.^model{1}.x0(whichVoxel);
predictionDuration = allstimimages*rf;

%Generate prediction for frequency response (monotonic)
rf2   = (params.analysis.Y.^model{1}.y0(whichVoxel))./params.analysis.Y;
predictionFrequency = allstimimages*rf2;

% fMRI predictions
%In resulting models, to make TSeries plots
colStr{1}='k';
colStr{2}='b';
colStr{3}='r';
colStr{4}='g';

predictionMono = predictionDuration*betas(1)+predictionFrequency*betas(2);


%% load tuned for selected voxel
close all

cd(fullfile(strcat(paths{subj_id},'/Gray/', cv_run, ' (L1)/SearchFitFreeExponent/')));
if string(cv_run) ~= "TimingSweeps"
    cd('xvalRefit');
end
model_name = '*2019*Lin-2dOvalGaussian-DurationPeriod-DT0.5-maxValue-2-expIntensity-free-fFit-fFit.mat';
load(strcat(ls(model_name)))
cd(paths{subj_id})

mrVista 3

%% make 2d for model
% first image without compressive exponent on frequency
rfTunedImg1= rfGaussian2d(modelDurations, modelPeriods',...
    model{1}.sigma.major(whichVoxel), ...
    model{1}.sigma.minor(whichVoxel), ...
    model{1}.sigma.theta(whichVoxel), ...
    model{1}.x0(whichVoxel), ...
    model{1}.y0(whichVoxel));
figure; imagesc(flipud(rfTunedImg1)); axis image
max_tuned1 = caxis;

% scale this image with the frequency component
freq=5./modelPeriods'; 
rfTunedImg2 = rfTunedImg1.*((freq.^model{1}.exp(whichVoxel))./freq);
figure; imagesc(flipud(rfTunedImg2)); axis image
max_tuned2 = caxis;

%% look at predictions Tuned per condition
rfTuned1= rfGaussian2d(durations1, periods1,...
    model{1}.sigma.major(whichVoxel), ...
    model{1}.sigma.minor(whichVoxel), ...
    model{1}.sigma.theta(whichVoxel), ...
    model{1}.x0(whichVoxel), ...
    model{1}.y0(whichVoxel));
freq1=5./periods1;
scale1=1./(freq1.^model{1}.exp(whichVoxel)); % compressive exponent thingie
scale1=scale1.*freq1;
rfTuned1=rfTuned1./scale1;
rfTuned1=rfTuned1/max_tuned2(2);

RF=[rfTuned1; rfTuned1];
RF=RF(:);
Xaxis=[1:26;1:26];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{1})
xticks(0.5:2:25.5)
xticklabels(durations1(1:2:end))
xlabel('Duration (s)')
ylabel('Response amplitude per event')
axis square;
xlim([0 26]);
ylim([0 1]);


rfTuned2= rfGaussian2d(durations2, periods2,...
    model{1}.sigma.major(whichVoxel), ...
    model{1}.sigma.minor(whichVoxel), ...
    model{1}.sigma.theta(whichVoxel), ...
    model{1}.x0(whichVoxel), ...
    model{1}.y0(whichVoxel));
freq2=5./periods2;
scale2=1./(freq2.^model{1}.exp(whichVoxel)); % compressive exponent thingie
scale2=scale2.*freq2;
rfTuned2=rfTuned2./scale2;
rfTuned2=rfTuned2/max_tuned2(2);


RF=[rfTuned2; rfTuned2];
RF=RF(:);
Xaxis=[1:26;1:26];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{2})
xticks(0.5:2:25.5)
xticklabels(durations2(1:2:end))
xlabel('Duration (s)')
ylabel('Response amplitude per event')
axis square;
xlim([0 26]);
ylim([0 1]);


rfTuned3= rfGaussian2d(durations3, periods3,...
    model{1}.sigma.major(whichVoxel), ...
    model{1}.sigma.minor(whichVoxel), ...
    model{1}.sigma.theta(whichVoxel), ...
    model{1}.x0(whichVoxel), ...
    model{1}.y0(whichVoxel));

freq3=5./periods3;
scale3=1./(freq3.^model{1}.exp(whichVoxel)); % compressive exponent thingie
scale3=scale3.*freq3;
rfTuned3=rfTuned3./scale3;
rfTuned3=rfTuned3/max_tuned2(2);

RF=[rfTuned3; rfTuned3];
RF=RF(:);
Xaxis=[1:26;1:26];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{3})
xticks(0.5:2:25.5)
xticklabels(durations3(1:2:end))
xlabel('Duration (s)')
ylabel('Response amplitude per event')
axis square;
xlim([0 26]);
ylim([0 1]);


rfTuned4= rfGaussian2d(durations4, periods4,...
    model{1}.sigma.major(whichVoxel), ...
    model{1}.sigma.minor(whichVoxel), ...
    model{1}.sigma.theta(whichVoxel), ...
    model{1}.x0(whichVoxel), ...
    model{1}.y0(whichVoxel));

freq4=5./periods4;
scale4=1./(freq4.^model{1}.exp(whichVoxel)); % compressive exponent thingie
scale4=scale4.*freq4;
rfTuned4=rfTuned4./scale4;
rfTuned4=rfTuned4/max_tuned2(2);

RF=[rfTuned4; rfTuned4];
RF=RF(:);
Xaxis=[1:26;1:26];
Xaxis=Xaxis(:);
Xaxis=[0; Xaxis(1:(end-1))];
figure; plot(Xaxis,RF, colStr{4})
xticks(0.5:2:25.5)
xticklabels(durations4(1:2:end))
xlabel('Duration (s)')
ylabel('Response amplitude per event')
axis square;
xlim([0 26]);
ylim([0 1]);



%% Tuned predictions total
rfTuned= rfGaussian2d(params.analysis.X, params.analysis.Y,...
    model{1}.sigma.major(whichVoxel), ...
    model{1}.sigma.minor(whichVoxel), ...
    model{1}.sigma.theta(whichVoxel), ...
    model{1}.x0(whichVoxel), ...
    model{1}.y0(whichVoxel));

freq=5./params.analysis.Y;
scale=1./(freq.^model{1}.exp(whichVoxel)); % compressive exponent thingie
scale=scale.*freq;
rfTuned=rfTuned./scale;

predictionTuned = allstimimages*rfTuned; % like *D in figure

% best scaling for the data
X=[predictionTuned ones(size(predictionTuned))];
betas=pinv(X)*voxelData; 

%224 (4x56)
predictionTuned=X*betas;


%% Rescale combined monotonic prediction
X=[predictionMono ones(size(predictionMono))];
betas=pinv(X)*voxelData; 

%224 (4x56)
predictionMono=X*betas;


%% Plot compressive exponent tuned
frequencies=0.01:0.01:20;
response=frequencies.^model{1}.exp(whichVoxel);
%response=frequencies.^1; %to check everything works, use linear function
figure; plot(frequencies, response)
axis square;
axis([0 20 0 20]);
xlabel('Frequency (Hz)')
ylabel('Response amplitude')

%% plot predictions and data
close all
scanTime=0:2.1:115.5;
smallest=0;

if max(predictionMono)>=max(voxelData) && max(predictionMono)>=max(predictionTuned)
    maxScaler=predictionMono;
elseif max(predictionTuned)>=max(voxelData) && max(predictionTuned)>=max(predictionMono)
    maxScaler=predictionTuned;
elseif max(voxelData)>=max(predictionMono) && max(voxelData)>=max(predictionTuned) 
    maxScaler=voxelData;
end

figure; hold on;
for n=1:4
    
        monoData=predictionMono(((n-1)*56+1):(n*56));
        dataData=voxelData(((n-1)*56+1):(n*56));
        tunedData=predictionTuned(((n-1)*56+1):(n*56));
        
        monoData=monoData-mean(maxScaler);
        monoData=monoData./(max(maxScaler)-mean(maxScaler));
        
        dataData=dataData-mean(maxScaler);
        dataData=dataData./(max(maxScaler)-mean(maxScaler));
        
        tunedData=tunedData-mean(maxScaler);
        tunedData=tunedData./(max(maxScaler)-mean(maxScaler));
        
        smallest=min([smallest, min(monoData) min(dataData) min(tunedData)]);
        
        plot(scanTime+1.05, monoData, [colStr{n} '-']);
        plot(scanTime+1.05, dataData, [colStr{n} 'o']);
        plot(scanTime+1.05, tunedData, [colStr{n} '--']);

end
axis square
set(gca, 'XTick', [0 29.4 58.8 88.2 117.6])
axis([0 117.6 smallest 1])
