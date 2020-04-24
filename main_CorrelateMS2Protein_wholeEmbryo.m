% DESCRIPTION
% Script to analyze the cross-correlation of mRNA and Protein for the whole
% embryo, for a given Prefix (Maybe we'd ultimately want to expand this for
% multiple embryos.

%
% ARGUMENTS
% DataStatusTab: master ID variable (should match a tab name in the Data Status
% sheet)
% DropboxFolder: full file path to folder containing compiled imaging
% results

% OPTIONS
% dropboxFolder: Pass this option, followed by the path to data folder 
%                where you wish to save
%                pipeline-generated data sets and figures. If this
%                var is not specified, output will be saved one level
%                above git repo in the folder structure
% first_nc: script defaults to taking only nc14 traces. If you want
%           earlier traces, pass 'first_nc', followed by desired nuclear cycle
%           number
%
% OUTPUT: nucleus_struct: compiled data set contain key nucleus and
% particle attributes


function main_CorrelateMS2Protein_wholeEmbryo(Prefix,varargin)
addpath('./utilities')
% set defaults
% firstNC = 14;
% minDP = 15;
% pctSparsity = 50;
% two_spot_flag = contains(DataStatusTab, '2spot');
% min_time = 0*60; % take no fluorescence data prior to this point
% TresInterp = 20; 
% calculatePSF = false;
% project = DataStatusTab;
% for i = 1:numel(varargin)
%     if ischar(varargin{i}) && i < numel(varargin) && mod(i,2)==1
%         eval([varargin{i} '=varargin{i+1};']);
%     end
% end
DropboxFolder = 'S:\YangJoon\Dropbox\CentralDogmaResults';

figPath = 'S:\YangJoon\Dropbox\CentralDogmaResults\CentralDogma_plots';

% Load data from a single embryo (Prefix)
processed_spot_data = load([DropboxFolder, filesep, Prefix, filesep, 'CompiledParticles.mat' ]); % processed particles  
processed_nucleus_data = load([DropboxFolder, filesep, Prefix, filesep, 'CompiledNuclei.mat' ]); % processed nuclei  
processed_schnitz_data = load([DropboxFolder, filesep, Prefix, filesep, Prefix, '_lin.mat' ]);        
frame_info = load([DropboxFolder, filesep, Prefix, filesep,'FrameInfo.mat' ]);   

nc12 = processed_spot_data.nc12;
nc13 = processed_spot_data.nc13;
nc14 = processed_spot_data.nc14;
%% Here, I want to connect the MeanVectorAll from particles(spots) and nuclei (schnitz)

% Note that I need to process the nuclear fluo by background-subtraction.
% there are at least a couple of ways for background-subtraction
% 1) subtracting the most posterior nuclear fluo at each time point
% 2) subtracting the cytoplasmic fluo * 1/K_G (1/0.8) for nuclear free
% eGFP.
% 3) For now, let's just see whether there's any steady-state

% Extract useful fields 
% Time
Time =  extractfield(frame_info.FrameInfo, 'Time');

% Particle fluo (over the whole imaging window : 20-60% of the EL)
MeanParticleFluoAll = processed_spot_data.MeanVectorAll{1,1};
SDParticleFluoAll = processed_spot_data.SDVectorAll{1,1};
NParticlesAll = processed_spot_data.NParticlesAll{1,1};
SEParticleFluoAll = SDParticleFluoAll./sqrt(NParticlesAll);

% Nuclear fluo (over the whole imaging window)
MeanNuclearFluoAll = processed_nucleus_data.MeanVectorAll;
SDNuclearFluoAll = processed_nucleus_data.SDVectorAll;
NNucleiAll = processed_nucleus_data.NParticlesAll;
SDNuclearFluoAll = SDNuclearFluoAll./sqrt(NNucleiAll);
%% Plot to check the spot fluo trend
% averagaed spot fluo trend
%errorbar(Time/60, MeanParticleFluoAll, SEParticleFluoAll)

% number of spots over time
%plot(Time/60, NParticlesAll)

% spot fluo * number of spots
errorbar(Time/60-8, MeanParticleFluoAll.*NParticlesAll, SEParticleFluoAll.*NParticlesAll)

title('')
xlabel(' time (min)')
ylabel('fluorescence (AU)')
xlim([0 max(Time/60)-8])

StandardFigure(gcf,gca)

saveas(gcf, [figPath, filesep, 'ParticleFluo_times_NSpots_wholeEmbryo_', Prefix, '.tif'])
saveas(gcf, [figPath, filesep, 'ParticleFluo_times_NSpots_wholeEmbryo_', Prefix, '.pdf'])
%% Integrate the spot fluo over time (area under the curve)
% First, we need to convert all NaNs to zeros
MeanParticleFluoAll(isnan(MeanParticleFluoAll))=0;
SDParticleFluoAll(isnan(SDParticleFluoAll)) = 0;
NParticlesAll(isnan(NParticlesAll))=0;

% Use the trapz for the time-integral
% Integrate from the beginning of NC13
tBegin = nc13;
clear integratedmRNA

for i=tBegin:length(Time)
    integratedmRNA(i) = trapz(Time(tBegin-1:i), MeanParticleFluoAll(tBegin-1:i).*NParticlesAll(tBegin-1:i));
    integratedmRNA_SD(i) = sqrt(trapz(Time(tBegin-1:i), SDParticleFluoAll(tBegin-1:i).^2.*NParticlesAll(tBegin-1:i)));
end

errorbar((Time(tBegin:end)-Time(tBegin))/60, integratedmRNA(tBegin:end), integratedmRNA_SD(tBegin:end))
xlabel(' time (min)')
ylabel(' integrated mRNA (AU)')

xlim([0 (max(Time)-Time(tBegin))/60])

StandardFigure(gcf,gca)

% saveas(gcf, [figPath, filesep, 'IntegratedmRNA_Allspots_wholeEmbryo_nc14_', Prefix, '.tif'])
% saveas(gcf, [figPath, filesep, 'IntegratedmRNA_Allspots_wholeEmbryo_nc14_', Prefix, '.pdf'])

%% Note that, in the whole embryo, we don't have to think about diffusion as the whole thing is somewhat confined.
% Then, the only parameters are half-life of mRNA and protein, and
% production rate of protein (averaged).
% In this case, we can also use 60 minutes as half-life of hb mRNA (Little, 2013)
% Then, we have only two free parameters, protein production and
% degradation rate.
%% Predict the integrated mRNA with mRNA half-life

% Predict the mRNA
tBegin = nc13;
clear integratedmRNA
integratedmRNA = zeros(1,length(Time));

Tm = 60; % mRNA half-life (Little, 2013)
gammaM = log(2)/Tm; % degradation constant

for i=tBegin+2:length(Time)
    dt = (Time(i) - Time(i-1))/60; % converting to minutes
    integratedmRNA(i) = trapz(Time(tBegin:i-1)/60,...
                MeanParticleFluoAll(tBegin:i-1).*NParticlesAll(tBegin:i-1)) - ...
                gammaM*integratedmRNA(i-1)*dt;
    integratedmRNA_SD(i) = sqrt(trapz(Time(tBegin:i-1)/60,...
                SDParticleFluoAll(tBegin:i-1).^2.*NParticlesAll(tBegin:i-1)));
end

errorbar((Time(tBegin:end)-Time(tBegin))/60, integratedmRNA(tBegin:end), integratedmRNA_SD(tBegin:end))
xlabel(' time (min)')
ylabel(' integrated mRNA (AU)')

xlim([0 (max(Time)-Time(tBegin))/60])

StandardFigure(gcf,gca)

% saveas(gcf, [figPath, filesep, 'IntegratedmRNA_Allspots_wholeEmbryo_nc14_', Prefix, '.tif'])
% saveas(gcf, [figPath, filesep, 'IntegratedmRNA_Allspots_wholeEmbryo_nc14_', Prefix, '.pdf'])

%% Predict the protein level from the integrated mRNA
% We define two parameters here which we can tune further.
% Ultimately, for a set of parameters, we need to find a range of
% parameters that has good correlations with the actual measured protein at
% each time point. 

% production rate: rp, half-life: Tp
rp = 2; % 2 protein molecules / (mrna * time)
Tp = 50; % 50 min (from Bicoid, Drocco paper)

gammaP = log(2)/Tp;

% Initialize the vector
PredictedProtein = zeros(1,length(Time));
PredictedProtein_SD = zeros(1,length(Time));

for i=tBegin+2:length(Time)
    dt = (Time(i) - Time(i-1))/60; % converting to minutes
    PredictedProtein(i) = PredictedProtein(i-1) + ...
                rp*integratedmRNA(i-1)*dt - ...
                gammaP*PredictedProtein(i-1)*dt;
    PredictedProtein_SD(i) = sqrt(trapz(Time(tBegin:i-1)/60, ...
                integratedmRNA_SD(tBegin:i-1).^2)*rp);
end

%% plot to check
errorbar((Time(tBegin:end)-Time(tBegin))/60, PredictedProtein(tBegin:end), PredictedProtein_SD(tBegin:end))
xlabel(' time (min)')
ylabel(' predicted protein (AU)')

xlim([0 (max(Time)-Time(tBegin))/60])

StandardFigure(gcf,gca)

saveas(gcf, [figPath, filesep, 'PredictedProtein_Allspots_wholeEmbryo_nc14_rp=2_Tp=50_', Prefix, '.tif'])
saveas(gcf, [figPath, filesep, 'PredictedProtein_Allspots_wholeEmbryo_nc14_rp=2_Tp=50_', Prefix, '.pdf'])

%% Plot the prediction and data on top
% First, I'll just subtract the protein fluo from the beginning of nc13 as
% an approximation for the background (assuming that the background doesn't
% fluctuate that much over time).
hold on
MeasuredProtein = MeanNuclearFluoAll.*NNucleiAll - ...
                    MeanNuclearFluoAll(tBegin).*NNucleiAll(tBegin);
MeasuredProtein_SD = sqrt(SDNuclearFluoAll.^2.*NNucleiAll);
scaling = 1/3*10^3;

errorbar((Time(tBegin:end)-Time(tBegin))/60, PredictedProtein(tBegin:end), PredictedProtein_SD(tBegin:end))
errorbar((Time(tBegin:end)-Time(tBegin))/60, scaling*MeasuredProtein(tBegin:end), scaling*MeasuredProtein_SD(tBegin:end))

xlabel(' time (min)')
ylabel(' protein (AU)')
xlim([0 (max(Time)-Time(tBegin))/60])
legend('Prediction','Measurement','Location','NorthWest')

StandardFigure(gcf,gca)
%saveas(gcf, [figPath, filesep, 'Protein_Time_Prediction_Measurement_wholeEmbryo_nc13_rp=2_Tp=50_', Prefix, '.tif'])
%saveas(gcf, [figPath, filesep, 'Protein_Time_Prediction_Measurement_wholeEmbryo_nc13_rp=2_Tp=50_', Prefix, '.pdf'])

%% Plot the prediction and measurement 
index = nc14+10:length(Time);
plothandle = errorbarxy(PredictedProtein(index), MeasuredProtein(index), PredictedProtein_SD(index), MeasuredProtein_SD(index))

% plothandle.hMain.Marker = 'o'
% plothandle.hMain.MarkerSize = 2
% plothandle.hMain.Parent.XScale = 'log'
% plothandle.hMain.Parent.YScale = 'log'
xlabel('predicted protein')
ylabel('measured protein')
StandardFigure(gcf,gca)
saveas(gcf, [figPath, filesep, 'Protein_Prediction_vs_Measurement_wholeEmbryo_nc14_5-40min_rp=2_Tp=50_', Prefix, '.tif'])
saveas(gcf, [figPath, filesep, 'Protein_Prediction_vs_Measurement_wholeEmbryo_nc14_5-40min_rp=2_Tp=50_', Prefix, '.pdf'])
end