function MakeValidationPredictions_anonymous(paths, subj_id)
%% Makes ground truth dataset
% paths_background: paths to where params, allstimimages, and TimingModelParamsMix are stored
% paths: where the information about each participant is stored 
% subj_id: index of subject you want (where in paths)

load('params.mat');
load('allstimimages.mat');
cd(paths{subj_id});

%% Monotonic
xs=[0.0000001 0.1:0.1:1];
ys=xs;
[xs ys]=meshgrid(xs, ys);
figure; plot(xs, ys,'.')
params.analysis.x0=xs(:);
params.analysis.y0=ys(:);
figure; plot(params.analysis.x0, params.analysis.y0,'.')
n = numel(params.analysis.x0);
s = [[1:ceil(n./1000):n-2] n+1]; %#ok<NBRAK>
allstimimages2=allstimimages;
%Duration component
prediction = zeros(size(allstimimages,1),n,'single');
fprintf(1,'[%s]:Making %d model samples:',mfilename,n);
drawnow;tic;
for n=1:numel(params.analysis.x0);
    rf=params.analysis.X.^params.analysis.x0(n);
    %         rf=rf-min(rf(:));
    %         rf=rf./max(rf(:));
    % store
    prediction(:,n) = allstimimages2*rf;
    %rf=repmat(rf', [224 1]);
end
%Frequency component
prediction2=zeros(size(prediction));
for n=1:numel(params.analysis.x0);
    rf2   = (params.analysis.Y.^params.analysis.y0(n))./params.analysis.Y;
    % convolve with stimulus
    %         rf2=rf2-min(rf2(:));
    %         rf2=rf2./max(rf2(:));
    prediction2(:,n) = allstimimages2*rf2;
end

params.analysis.ratio=[0:0.1:1 1./(0.9:-0.1:0.1)]';

predictions=[];
for n=1:length(params.analysis.ratio)
    predictions=[predictions prediction*params.analysis.ratio(n)+prediction2];
end

params.analysis.x0=repmat(params.analysis.x0, [size(params.analysis.ratio, 1) 1]);
params.analysis.y0=repmat(params.analysis.y0, [size(params.analysis.ratio, 1) 1]);
tmp=repmat(params.analysis.ratio, [1 length(xs(:))]);
tmp=tmp';
tmp=tmp(:);
params.analysis.ratio=tmp;

% %Normalise response by maximum response amplitude.
% for n=1:size(predictions, 2)
%     predictions(:,n)=predictions(:,n)./max(predictions(:,n));
% end

%Normalise by standard deviation of response amplitudes, so we can
%calculate noise level easily
for n=1:size(predictions, 2)
    predictions(:,n)=predictions(:,n)-mean(predictions(:,n));
    predictions(:,n)=predictions(:,n)./std(predictions(:,n));
end

noiseSDs=[0:0.1:6]';
predictionsNoiseOdd=[];
predictionsNoiseEven=[];
for n=1:length(noiseSDs)
    predictionsNoiseOdd=[predictionsNoiseOdd predictions+noiseSDs(n).*randn(size(predictions))];
    predictionsNoiseEven=[predictionsNoiseEven predictions+noiseSDs(n).*randn(size(predictions))];
end

%To test expected variance explained for a particular SD
whichNoiseSD=1;
%Expected variance explained in non-cross-validated data
where=size(predictions, 2)*whichNoiseSD; corrcoef(predictionsNoiseOdd(:,(where-2420)+1:where), predictions)
%Correlation between odd and even
where=size(predictions, 2)*whichNoiseSD; corrcoef(predictionsNoiseOdd(:,(where-2420)+1:where), predictionsNoiseEven(:,(where-2420)+1:where))

%Rescale predicitons to roughly follow data amplitudes
PercentageBOLD=2/3; %Standard deviation
MeanBOLD=60000;
predictionsNoiseOdd=predictionsNoiseOdd.*(MeanBOLD*PercentageBOLD/100)+MeanBOLD;
predictionsNoiseEven=predictionsNoiseEven.*(MeanBOLD*PercentageBOLD/100)+MeanBOLD;

tmp=repmat(noiseSDs, [1 length(params.analysis.x0(1:size(predictions, 2)))]);
tmp=tmp';
tmp=tmp(:);
params.analysis.noiseSDs=tmp;
params.analysis.x0=repmat(params.analysis.x0(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);
params.analysis.y0=repmat(params.analysis.y0(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);
params.analysis.ratio=repmat(params.analysis.ratio(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);

%% Put in right format
%Duplicate a data type and load tSeries1.mat of first run
load(fullfile(paths{subj_id},'/Gray/ValidationMonotonicOdd/TSeries/Scan1/tSeries1.mat'));
tSeries(:, 1:size(predictionsNoiseOdd, 2))=predictionsNoiseOdd(1:56, :);
save(strcat(paths{subj_id},'/Gray/ValidationMonotonicOdd/TSeries/Scan1/tSeries1.mat'),'tSeries');
%load tSeries1.mat of other runs run, modify and save
load(fullfile(paths{subj_id},'/Gray/ValidationMonotonicOdd/TSeries/Scan2/tSeries1.mat'));
tSeries(:, 1:size(predictionsNoiseOdd, 2))=predictionsNoiseOdd(57:112, :);
save(strcat(paths{subj_id},'/Gray/ValidationMonotonicOdd/TSeries/Scan2/tSeries1.mat'),'tSeries');

load(fullfile(paths{subj_id},'/Gray/ValidationMonotonicOdd/TSeries/Scan3/tSeries1.mat'));
tSeries(:, 1:size(predictionsNoiseOdd, 2))=predictionsNoiseOdd(113:168, :);
save(strcat(paths{subj_id},'/Gray/ValidationMonotonicOdd/TSeries/Scan3/tSeries1.mat'),'tSeries');

load(fullfile(paths{subj_id},'/Gray/ValidationMonotonicOdd/TSeries/Scan4/tSeries1.mat'));
tSeries(:, 1:size(predictionsNoiseOdd, 2))=predictionsNoiseOdd(169:224, :);
save(strcat(paths{subj_id},'/Gray/ValidationMonotonicOdd/TSeries/Scan4/tSeries1.mat'),'tSeries')

%Now for the other validation split
load(fullfile(paths{subj_id},'/Gray/ValidationMonotonicEven/TSeries/Scan1/tSeries1.mat'));
tSeries(:, 1:size(predictionsNoiseEven, 2))=predictionsNoiseEven(1:56, :);
save(strcat(paths{subj_id},'/Gray/ValidationMonotonicEven/TSeries/Scan1/tSeries1.mat'),'tSeries');
%load tSeries1.mat of other runs run, modify and save
load(fullfile(paths{subj_id},'/Gray/ValidationMonotonicEven/TSeries/Scan2/tSeries1.mat'));
tSeries(:, 1:size(predictionsNoiseEven, 2))=predictionsNoiseEven(57:112, :);
save(strcat(paths{subj_id},'/Gray/ValidationMonotonicEven/TSeries/Scan2/tSeries1.mat'),'tSeries');

load(fullfile(paths{subj_id},'/Gray/ValidationMonotonicEven/TSeries/Scan3/tSeries1.mat'));
tSeries(:, 1:size(predictionsNoiseEven, 2))=predictionsNoiseEven(113:168, :);
save(strcat(paths{subj_id},'/Gray/ValidationMonotonicEven/TSeries/Scan3/tSeries1.mat'),'tSeries');

load(fullfile(paths{subj_id},'/Gray/ValidationMonotonicEven/TSeries/Scan4/tSeries1.mat'));
tSeries(:, 1:size(predictionsNoiseEven, 2))=predictionsNoiseEven(169:224, :);
save(strcat(paths{subj_id},'/Gray/ValidationMonotonicEven/TSeries/Scan4/tSeries1.mat'),'tSeries');

%Now make and ROI to run this in
load(fullfile(paths{subj_id},'/Gray/ROIs/gray-Layer1.mat'));
load(fullfile(paths{subj_id},'/Gray/coords.mat'));
ROI.coords=coords(:,1:size(predictionsNoiseOdd, 2));
ROI.name='ValidationMonotonic';
save(strcat(paths{subj_id},'/Gray/ROIs/ValidationMonotonic.mat'), 'ROI');

%Now make an ROI for only the no-noise data to run this in
load(strcat(paths{subj_id},'/Gray/ROIs/gray-Layer1.mat'));
load(strcat(paths{subj_id},'/Gray/coords.mat'));
ROI.coords=coords(:,1:2420);
ROI.name='ValidationMonotonicNoNoise';
save(strcat(paths{subj_id},'/Gray/ROIs/ValidationMonotonicNoNoise.mat'), 'ROI');


%% Tuned model. First get predictions
load('TimingModelParamsMix.mat');
combinedDT=5:7;
setAllRetParams(paramsDurationPeriod, combinedDT);
VOLUME{1} = initHiddenGray;
%Run to breakpoint in rmGridFit line 265
rmRunDuration2dOval(VOLUME{1},combinedDT, 'gray-Layer1',5,{'1g'},[],[],[],0,'free'); %This is the best fitting model. Setting the 4th input arguement to '5' fits the HRF parameters too.

%Now choose some values to make predictions for
xy=unique([params.analysis.x0(:) params.analysis.y0(:)], 'rows');
% xs=xy(:,1)';
% ys=xy(:,2)';
sMaj=[0.05 0.5 1 1.5 2];
sMin=sMaj;
theta=unique(params.analysis.theta)';
exp=unique(params.analysis.exponent)';

comb=combvec(xy', sMaj, sMin, theta, exp);
params.analysis.x0=comb(1,:)';
params.analysis.y0=comb(2,:)';
params.analysis.sigmaMajor=comb(3,:)';
params.analysis.sigmaMinor=comb(4,:)';
params.analysis.theta=comb(5,:)';
params.analysis.exponent=comb(6,:)';

%Now run rmGridFit to line 313 to generate predictions

%Then:
predictions=prediction;
%Normalise by standard deviation of response amplitudes, so we can
%calculate  easily
for n=1:size(predictions, 2)
    predictions(:,n)=predictions(:,n)-mean(predictions(:,n));
    predictions(:,n)=predictions(:,n)./std(predictions(:,n));
end

noiseSDs=[0:0.2:6]';
predictionsNoiseOdd=[];
predictionsNoiseEven=[];
for n=1:length(noiseSDs)
    predictionsNoiseOdd=[predictionsNoiseOdd predictions+noiseSDs(n).*randn(size(predictions))];
    predictionsNoiseEven=[predictionsNoiseEven predictions+noiseSDs(n).*randn(size(predictions))];
end

%To test expected variance explained for a particular SD
whichNoiseSD=31;
%Expected variance explained in non-cross-validated data
where=size(prediction, 2)*whichNoiseSD; corrcoef(predictionsNoiseOdd(:,(where-size(prediction, 2)+1):where), predictions)
%Correlation between odd and even
where=size(prediction, 2)*whichNoiseSD; corrcoef(predictionsNoiseOdd(:,(where-size(prediction, 2)+1):where), predictionsNoiseEven(:,(where-size(prediction, 2)+1):where))

%Rescale predicitons to roughly follow data amplitudes
PercentageBOLD=2/3; %Standard deviation
MeanBOLD=60000;
predictionsNoiseOdd=predictionsNoiseOdd.*(MeanBOLD*PercentageBOLD/100)+MeanBOLD;
predictionsNoiseEven=predictionsNoiseEven.*(MeanBOLD*PercentageBOLD/100)+MeanBOLD;

tmp=repmat(noiseSDs, [1 length(params.analysis.x0(1:size(predictions, 2)))]);
tmp=tmp';
tmp=tmp(:);
params.analysis.noiseSDs=tmp;
params.analysis.x0=repmat(params.analysis.x0(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);
params.analysis.y0=repmat(params.analysis.y0(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);
params.analysis.sigmaMajor=repmat(params.analysis.sigmaMajor(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);
params.analysis.sigmaMinor=repmat(params.analysis.sigmaMinor(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);
params.analysis.theta=repmat(params.analysis.theta(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);
params.analysis.exponent=repmat(params.analysis.exponent(1:size(predictions, 2)), [size(noiseSDs, 1) 1]);

%% Save everything
save('ValidationTuned.mat', 'predictions', 'params', 'predictionsNoiseOdd', 'predictionsNoiseEven', 'prediction')

%Exit rmGridFit
dbquit;

load('ValidationTuned.mat')
howManyDTs=ceil(size(predictionsNoiseEven, 2)/size(VOLUME{1}.coords, 2));
voxelsPerDT=ceil(size(predictionsNoiseEven, 2)/howManyDTs);
lastDT=length(dataTYPES);
for n=1:howManyDTs
    duplicateDataType(VOLUME{1}, ['ValidationTunedOdd' num2str(n)], lastDT)
    duplicateDataType(VOLUME{1}, ['ValidationTunedEven' num2str(n)], lastDT)
end

%% put in right format 
%Duplicate a data type and load tSeries1.mat of first run
TRsPerRun=56;
for n=1:howManyDTs
    for scan=1:4
        load(['Gray/ValidationTunedOdd', num2str(n),'/TSeries/Scan', num2str(scan), '/tSeries1.mat']);
        tSeries(:, 1:voxelsPerDT)=predictionsNoiseOdd(((scan-1)*TRsPerRun+1):(scan*TRsPerRun), ((n-1)*voxelsPerDT+1):(n*voxelsPerDT));
        save(['Gray/ValidationTunedOdd', num2str(n),'/TSeries/Scan', num2str(scan), '/tSeries1.mat'],'tSeries');
        load(['Gray/ValidationTunedEven', num2str(n),'/TSeries/Scan', num2str(scan), '/tSeries1.mat']);
        tSeries(:, 1:voxelsPerDT)=predictionsNoiseEven(((scan-1)*TRsPerRun+1):(scan*TRsPerRun), ((n-1)*voxelsPerDT+1):(n*voxelsPerDT));
        save(['Gray/ValidationTunedEven', num2str(n),'/TSeries/Scan', num2str(scan), '/tSeries1.mat'],'tSeries');
    end
end

load(strcat(paths{subj_id},'/Gray/ROIs/gray-Layer1.mat'));
load(strcat(paths{subj_id},'/Gray/coords.mat'));
ROI.coords=coords(:,1:voxelsPerDT);
ROI.name='ValidationTuned';
save(strcat(paths{subj_id},'/Gray/ROIs/ValidationTuned.mat'), 'ROI');
end
%% RAN MODELS USING CROSSVALIDATE2COMPONENTS
