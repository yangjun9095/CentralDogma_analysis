% projectName = 'CentralDogma'
% Last updated : 
% Description : For a given Prefix (single embryo),
% generate plots of raw data, raw MS2 time traces, spatial profile, as well
% as LlamaTag data.

% Inputs : 1) single embryo (Prefix), 2) DataType (LoadMS2Sets) or
% AveragedDatasets, 3) individual embryos (optional...)

% Outputs : 

function plt_rawData_MS2_LlamaTag(Prefix, varargin)
%% Load datasets
FilePath = 'S:\YangJoon\Dropbox\CentralDogmaResults';
Prefix = '2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-09-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1-3';
cp = load([FilePath,filesep,...
            Prefix,filesep,'CompiledParticles.mat']);
        
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
%% Extract useful fields from the cp
% For now, we will think about nc13, and nc14. This can be (and
% should be) relaxed to include specific NCs

% time info
Time = cp.ElapsedTime;
nc12 = cp.nc12;
nc13 = cp.nc13;
nc14 = cp.nc14;

% mean spot fluorescence info
% Note that the compiledparticles has fields inside a cell in some cases.
if iscell(cp.MeanVectorAP)
    spotfluo_mean = cell2mat(cp.MeanVectorAP);
    spotfluo_sd = cell2mat(cp.SDVectorAP);
    num_particles = cell2mat(cp.NParticlesAP);
    spotfluo_sem =  spotfluo_sd./num_particles;
else
    spotfluo_mean = (cp.MeanVectorAP);
    spotfluo_sd = (cp.SDVectorAP);
    num_particles = (cp.NParticlesAP);
    spotfluo_sem =  spotfluo_sd./num_particles;
end

%% generate plots of MS2 time traces
hold on
plot(Time(nc13:end) - Time(nc13), spotfluo_mean(nc13:end,:))

%% generate plots of MS2 time traces for one AP bin
APbin = 11; % 
errorbar(Time(nc13:end) - Time(nc13), ...
        spotfluo_mean(nc13:end,APbin),...
        spotfluo_sem(nc13:end,APbin))
title('MS2 fluorescence - 25% EL')
xlabel('time (min)')
ylabel('MS2 spot fluorescence (AU)')
StandardFigure(gcf,gca)

% save the plot
saveas(gcf, [figPath, filesep, 'MS2_spotfluo_time_25%','.tif'])
saveas(gcf, [figPath, filesep, 'MS2_spotfluo_time_25%','.pdf'])

%% generate plots of MS2 spatial profile
APaxis = 0:0.025:1;
% 0,10,20,30, minutes into nc14
tStep = floor(10/0.6);
tvec = [nc14:tStep:length(Time),length(Time)]
hold on
for t=tvec
errorbar(APaxis, spotfluo_mean(t,:), spotfluo_sem(t,:) )
end
xlim([0.2 0.6])
ylim([0 350])
xticks([0.2 0.3 0.4 0.5 0.6])
xticklabels([20 30 40 50 60])
title('MS2 fluorescence')
xlabel('embryo length (%)')
ylabel('MS2 spot fluorescence (AU)')
legend('0 min','10 min','20 min','30 min','35 min')
StandardFigure(gcf,gca)

%% import the estimated cytomRNA and Protein
APbin = 11; % 20%
tWindow = 1:60:length(Time_interp);
errorbar(Time_interp(tWindow),...
        AccumulatedmRNA(tWindow,APbin),...
        ErrorAccumulatedmRNA(tWindow,APbin),'r')
title('Accumulated mRNA - 25% EL')
xlabel('time into nc 13(min)')
ylabel('Accumulated mRNA (AU)')
StandardFigure(gcf,gca)
%% protein - prediction
APbin = 11; % 20%
tWindow = 1:60:length(Time_interp);
errorbar(Time_interp(tWindow),...
        Protein(tWindow,APbin),...
        ErrorProtein(tWindow,APbin))
title(' protein-prediction- 25% EL')
xlabel('time into nc 13(min)')
ylabel('protein (AU)')
StandardFigure(gcf,gca)

%% import/load the protein datasets
% Load the dataset
cn = load([FilePath,filesep,...
            Prefix,filesep,'CompiledNuclei.mat']);

% Extract useful fields
% mean nuclear fluo, sd, number of nuclei, time
% time info
% Time = cn.ElapsedTime;
tLength = length(cn.MeanVectorAP(:,1));
nc12 = cn.nc12;
nc13 = cn.nc13;
nc14 = cn.nc14;
Time_nuc = cn.ElapsedTime;
tRes = median(diff(Time_nuc));

% mean, sd, number of nuclei
nucfluo_mean = cn.MeanVectorAP;
nucfluo_sd = cn.SDVectorAP;
num_nuclei = cn.NParticlesAP;

nucfluo_sem = nucfluo_sd./num_nuclei;

%% subtract the Background fluo (free eGFP)
% Note that there are multiple ways to do so and it's summarized in my
% OpposingGradients notes. Here, I'll just pick the simplest way, which is
% picking the most posterior bins (50%-60%) that doesn't have Nans for all
% time points. Perhaps we can average 2-3 APbins from the most posterior
% bin.

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

%%
APbin = 11;%

errorbar(Time(nc13:end) - Time(nc13),...
            nucfluo_mean_BGsubt(nc13:length(Time),APbin),...
            nucfluo_sem(nc13:length(Time),APbin))
ylim([0 700])
title('protein-LlamaTag - 25% EL')
xlabel('time (min)')
ylabel('protein (AU)')
StandardFigure(gcf,gca)

% save the plot
saveas(gcf, [figPath, filesep, 'LlamaTag_time_25%','.tif'])
saveas(gcf, [figPath, filesep, 'LlamaTag_time_25%','.pdf'])
end