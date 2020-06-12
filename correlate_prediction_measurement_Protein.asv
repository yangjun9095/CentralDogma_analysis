% projectName = 'CentralDogma'
% Last updated : 6/5/2020
% Description : This script is for calculating the correlation between the
% prediction and measurement of protein at different time points.

% Inputs : 1) single embryo (Prefix), 2) DataType (LoadMS2Sets) or
% AveragedDatasets, 3) individual embryos (optional...)

% Caveat : time-delay between the prediction and measurement?

%% Load the dataset
FilePath = 'S:\YangJoon\Dropbox\CentralDogmaResults';
Prefix = '2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1';
% cp = load([FilePath,filesep,...
%             Prefix,filesep,'CompiledParticles.mat']);

%% Step1. Calculate the measured protein (background subtracted)
% To be edited : 
% Make this as a sub-function that has the inputs of Prefix, and filepath, then
% returns Time, NucFluo, NucFluoError (after BG subtraction, with a given
% method).

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

% method. Subtract the nuc fluo from the most posterior AP bins
% subtract the Background fluo (free eGFP)
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

% There shouldn't be negative TF concentration, we will deal with those
% better in BG subtraction. For now, let's put NaNs to those AP bins.
nucfluo_mean_BGsubt(nucfluo_mean_BGsubt<=0) = nan;

%% Trim the measured protein with the startNC
nucfluo_mean_BGsubt = nucfluo_mean_BGsubt(nc13:end,:);
nucfluo_sem = nucfluo_sem(nc13:end,:);
% redefine the NCs
nc13 = 1;
nc14 = cn.nc14 - cn.nc13;

%% Step2. Calculate the predicted protein

% Diffusion
Dm = [0,1,10,100]; % diffusion constant (um^2/sec)
Dp = [0,1,10,100];% diffusion constant (um^2/sec)

% half-life
T_half_mRNA = [1,10,100,1000];% min
T_half_protein = [1,10,100,1000]; % min


params = [Dm, T_half_mRNA, Dp, T_half_protein];

% This step of defining the parameter should be done in a more systematic
% way. For example, for picking values one by one from n-dimensional matrix

for i=1:length(Dm)
    for j=1:length(Dp)
        for k=1:length(T_half_mRNA)
            for l=1:length(T_half_protein)
                % define the parameters
                params = [Dm(i), T_half_mRNA(k), Dp(j), T_half_protein(l)];
                % calculate the correlation coefficient (Pearson)
                [Pearson_corr,~,~,~,~] = calculate_PearsonCorr_prediction_measurement_Protein(Prefix,...
                                                params,nucfluo_mean_BGsubt,nucfluo_sem);
                Corr_Coeff(i,j,k,l) = Pearson_corr;
            end
        end
    end
end

%% Gunawardena/Depace, Estrada 2016 approach?
% can we visualize the results in 2D?

%% Plot the Correlation Coefficient for a combination of 3 parameters
% pick the mRNA half life of 100 min,
Corr_hb = squeeze(Corr_Coeff(:,:,4,:))
% [X,Y,Z] = ndgrid(1:size(Corr_hb,1), 1:size(Corr_hb,2), 1:size(Corr_hb,3));
hm_cm = flipud(brewermap(1000,'Spectral'));


hold on
colorbar(hm_cm);
for i=1:length(Dm)
    for j=1:length(Dp)
        for l=1:length(T_half_protein)
            scatter3(Dm(i), Dp(j), T_half_protein(l),'o','MarkerFaceColor',hm_cm(floor(Corr_hb(i,j,l)*1000)+1,:),'MarkerEdgeColor','black')
        end
    end
end
view(-30,30)

grid on
xlabel('D_{m}')
ylabel('D_{p}')
zlabel('T_{1/2, protein}')
set(gca,'XScale','log')
set(gca,'YScale','log')
set(gca,'ZScale','log')

%% plot in 2D
[X,Y] = meshgrid(Dp,T_half_protein)
Z = squeeze(Corr_Coeff(1,3,:,:));
surf(X,Y,squeeze(Corr_Coeff(1,3,:,:)))



%% (Optional) Plotting the Prediciton/Measurement/Correlation for a couple of set of parameters as an example. 
%% generate plots from a set of parameters (hb and/or eve)
% first, define the parameters for hb and eve
params_hb = [0, 60, 7, 50]; % [Dm, T_half_mRNA, Dp, T_half_protein]
params_eve = [0, 7, 7, 7]; % [Dm, T_half_mRNA, Dp, T_half_protein]

% Second, calculate the Corr, Prediction, Measurement for these sets of
% parameters.
[Corr_Coeff_hb, PredictedProtein_hb, PredictedProtein_error_hb,...
                MeasuredProtein, MeasuredProtein_error] = ...
                                calculate_PearsonCorr_prediction_measurement_Protein(Prefix,...
                                params_hb,nucfluo_mean_BGsubt,nucfluo_sem);

[Corr_Coeff_eve, PredictedProtein_eve, PredictedProtein_error_eve,~, ~] =...
                                calculate_PearsonCorr_prediction_measurement_Protein(Prefix,...
                                params_eve,nucfluo_mean_BGsubt,nucfluo_sem);
%% Plotting module
addpath('../utilities');

nDataPoints = 10;
% color map
hm_cm = flipud(brewermap(nDataPoints,'Spectral'));
colormap(hm_cm);

%% Plot the Prediction and Measurement over AP axis

% Define the AP axis
APaxis = 0:0.025:1;

% plot every 10 minutes from 5 minutes, so, 5, 15, 25, and 25 minutes.
% framePoints

% First, the predicted protein - hb
predict_fig = figure(1)
hold on
for frame = 1:length(framePoints)
    errorbar(APaxis, PredictedProtein_hb(frame,:), PredictedProtein_error_hb(frame,:))
end
xlim([0.2 1])
ylim([0 max(max(PredictedProtein))+30000])
xticks([0.2 0.4 0.6 0.8  1])
xticklabels([20 40 60 80 100])
% xticks([0.2 0.3 0.4 0.5 0.6])
% xticklabels([20 30 40 50 60])
title('protein-predicted')
xlabel('embryo length (%)')
ylabel('protein concentration(AU)')
legend('5 min','15 min','25 min','35 min')
StandardFigure(predict_fig,predict_fig.CurrentAxes)

figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
saveas(predict_fig, [figPath, filesep, 'protein_predicted_hbCondition','.tif'])
saveas(predict_fig, [figPath, filesep, 'protein_predicted_hbCondition','.pdf'])

% Second, the predicted protein - eve
predict_fig_eve = figure(2)
hold on
for frame = 1:length(framePoints)
    errorbar(APaxis, PredictedProtein_hb(frame,:), PredictedProtein_error_hb(frame,:))
end
xlim([0.2 1])
ylim([0 max(max(PredictedProtein))+30000])
xticks([0.2 0.4 0.6 0.8  1])
xticklabels([20 40 60 80 100])
% xticks([0.2 0.3 0.4 0.5 0.6])
% xticklabels([20 30 40 50 60])
title('protein-predicted')
xlabel('embryo length (%)')
ylabel('protein concentration(AU)')
legend('5 min','15 min','25 min','35 min')
StandardFigure(predict_fig_eve,predict_fig_eve.CurrentAxes)

figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
saveas(predict_fig_eve, [figPath, filesep, 'protein_predicted_eveCondition','.tif'])
saveas(predict_fig_eve, [figPath, filesep, 'protein_predicted_eveCondition','.pdf'])

% Third, the measured protein (This has to be revisited after a more
% rigorous characterization of BG subtraction).
% First, the predicted protein
measure_fig = figure(3)
hold on
for frame = 1:length(framePoints)
    errorbar(APaxis, MeasuredProtein(frame,:), MeasuredProtein_error(frame,:))
end
xlim([0.2 1])
ylim([0 max(max(MeasuredProtein))+max(max(MeasuredProtein))/10])
xticks([0.2 0.4 0.6 0.8  1])
xticklabels([20 40 60 80 100])
% xticks([0.2 0.3 0.4 0.5 0.6])
% xticklabels([20 30 40 50 60])
title('protein-measured (LlamaTag)')
xlabel('embryo length (%)')
ylabel('protein concentration(AU)')
legend('5 min','15 min','25 min','35 min')
StandardFigure(measure_fig,measure_fig.CurrentAxes)

figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
saveas(measure_fig, [figPath, filesep, 'protein_measured','.tif'])
saveas(measure_fig, [figPath, filesep, 'protein_measured','.pdf'])
%% Plot the Prediction and Measurement to see the correlation

% first, color for time points
% color map
hm_cm = flipud(brewermap(4,'Spectral'));

% hb condition
hold on
for tpoint = 1:length(framePoints)
    plot(PredictedProtein_hb(tpoint,:), MeasuredProtein(tpoint,:),...
            'o','MarkerFaceColor',hm_cm(tpoint,:),'MarkerEdgeColor','black', 'MarkerSize', 10)
end

% labels
title('protein : prediction vs measurement')
xlabel('prediction (AU)')
ylabel('measurement (AU)')
legend('5','15','25','35','Location','NorthWest')

% Axes Ticks
xticks([0 0.4 0.8 1.2 1.6]*10^6)
yticks([100 200 300 400 500])

StandardFigure(gcf,gca)
% Save the plot
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
saveas(gcf, [figPath, filesep, 'correlation_prediction_measurement_timepoints_nc14_hbCondition','.tif'])
saveas(gcf, [figPath, filesep, 'correlation_prediction_measurement_timepoints_nc14_hbCondition','.pdf'])

%% Plot the Prediction and Measurement to see the correlationeve condition - eve
hold on
for tpoint = 1:length(framePoints)
    plot(PredictedProtein_eve(tpoint,:), MeasuredProtein(tpoint,:),...
            'o','MarkerFaceColor',hm_cm(tpoint,:),'MarkerEdgeColor','black', 'MarkerSize', 10)
end

% labels
title('protein : prediction vs measurement')
xlabel('prediction (AU)')
ylabel('measurement (AU)')
legend('5','15','25','35','Location','NorthWest')

% Axes Ticks
%xticks([0 0.4 0.8 1.2 1.6 2.0 2.4]*10^6)
yticks([100 200 300 400 500])

StandardFigure(gcf,gca)
% Save the plot
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
saveas(gcf, [figPath, filesep, 'correlation_prediction_measurement_timepoints_nc14_eveCondition','.tif'])
saveas(gcf, [figPath, filesep, 'correlation_prediction_measurement_timepoints_nc14_eveCondition','.pdf'])
%% 