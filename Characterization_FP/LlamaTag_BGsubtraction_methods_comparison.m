% projectName = 'CentralDogma'
% Last updated : 6/2/2020
% Description : For a given Prefix (single embryo), subtract the background
% fluo from the nuclear fluo, for the LlamaTag datasets. 
% Note that the cytoplasmic/nuclear ratio of TF can be TF specific.

% Inputs : 1) single embryo (Prefix), 2) DataType (LoadMS2Sets) or
% AveragedDatasets, 3) individual embryos (optional...)

%% Description
% After this script is somewhat finalized, I want to either rename/split
% the code with analysis part, and plot generation part. Perhaps saving the
% results of this analysis into a .mat structure for all different
% approaches would be great for plotting in different ways.

%% Load a dataset with LlamaTag
FilePath = 'S:\YangJoon\Dropbox\CentralDogmaResults';
Prefix = '2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-09-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1-3';
cn = load([FilePath,filesep,...
            Prefix,filesep,'CompiledNuclei.mat']);

%% Extract useful fields
% mean nuclear fluo, sd, number of nuclei, time
% time info
% Time = cn.ElapsedTime;
tLength = length(cn.MeanVectorAP(:,1));
nc12 = cn.nc12;
nc13 = cn.nc13;
nc14 = cn.nc14;
nucTime = cn.ElapsedTime;

% mean, sd, number of nuclei
nucfluo_mean = cn.MeanVectorAP;
nucfluo_sd = cn.SDVectorAP;
num_nuclei = cn.NParticlesAP;

nucfluo_sem = nucfluo_sd./num_nuclei;

%% Step1. Subtract the nuc fluo from the most posterior AP bins
% Here, I'll just pick the simplest way, which is picking the most posterior bins (50%-60%)
% that doesn't have Nans for all time points. 
% Perhaps we can average 2-3 APbins from the most posterior bin.

% the first time point, find the APbins with NaNs. Then, pick the posterior
% one that is larger than 9th bin.
tPoint = 1;
nanAPbins = find(isnan(nucfluo_mean(tPoint,:)));
% among the posterior bins, pick the most anterior.
nanAPIndex = min(find(nanAPbins>9));
nanAPbinIndex = nanAPbins(nanAPIndex) - 1;

nanAPbins = [nanAPbinIndex-1, nanAPbinIndex];

% At each time point, subtract the nuc fluo of those most posterior APbins
for i=1:tLength
    nucfluo_mean_BGsubt(i,:) = nucfluo_mean(i,:) - nanmean(nucfluo_mean(i,nanAPbins));
end

%% generate plots for check (before BG subtraction)
APaxis = 0:0.025:1;

% take every tInterval minutes
tInterval = 10; % min
tRes = median(diff(nucTime));
tStep = floor(tInterval/tRes);

hold on
for t=[nc14+floor(5/tRes):tStep:tLength]
    %errorbar(APaxis, nucfluo_mean_BGsubt(t,:), nucfluo_sem(t,:))
    errorbar(APaxis, nucfluo_mean(t,:), nucfluo_sem(t,:))
end
xlim([0.2 0.6])
ylim([0 max(nucfluo_mean(tLength,:))+100])
xticks([0.2 0.3 0.4 0.5 0.6])
xticklabels([20 30 40 50 60])
title('protein-measured(LlamaTag)')
xlabel('embryo length (%)')
ylabel('protein concentration(AU)')
legend('5 min','15 min','25 min','35 min')
StandardFigure(gcf,gca)

%% generate plots for check (after BG subtraction)
hold on
for t=[nc14+floor(5/tRes):tStep:tLength]
    errorbar(APaxis, nucfluo_mean_BGsubt(t,:), nucfluo_sem(t,:))
end
xlim([0.2 1])
ylim([0 max(nucfluo_mean_BGsubt(tLength,:))+100])
xticks([0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
xticklabels([20 30 40 50 60 70 80 90 100])
% xticks([0.2 0.3 0.4 0.5 0.6])
% xticklabels([20 30 40 50 60])
title('protein-measured(LlamaTag)')
xlabel('embryo length (%)')
ylabel('protein concentration(AU)')
legend('5 min','15 min','25 min','35 min')
StandardFigure(gcf,gca)

%% Step2. Subtract the nuc fluo from the NoNB datasets
%% Load the dataset from NoNB
% Note. This dataset is processed using the AverageDatasets.m
% that averages multiple embryos, with synchronization from the nc13.
% also, the time resolution is 1 min, thus we need to interpolate.
NoNBdata = load('S:\YangJoon\Dropbox\CentralDogmaResults\NoNB-eGFP\NoNB-Averaged.mat');

% extract the useful fields
NoNB_time = NoNBdata.ElapsedTime;
nucfluo_NoNB = NoNBdata.MeanVectorAP;
nucfluo_error_NoNB = NoNBdata.SEVectorAP; % write dowm how this SEM is calculated somewhere.

nc13_NoNB = NoNBdata.nc13;
nc14_NoNB = NoNBdata.nc14

% generate the time series with tRes from the LlamaTag data
NewTime = NoNB_time(1):tRes:NoNB_time(end);

% interpolate the NoNB datsets with the new Time vector
%% Step3. Subtract the free eGFP using the cytoplasmic fluo

