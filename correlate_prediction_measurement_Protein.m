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

% First, define the set of parameters for prediction
% first, hb condition
% Dm = 0; % diffusion constant (um^2/sec)
% T_half_mRNA = 60; % min
% Dp = 7; % diffusion constant (um^2/sec)
% T_half_protein = 50; % min

% second, eve condition
% Dm = 0; % diffusion constant (um^2/sec)
% T_half_mRNA = 7; % min
% Dp = 7; % diffusion constant (um^2/sec)
% T_half_protein = 7; % min

Dm = [0,1,10,100];
T_half_mRNA = [1,10,100,1000];%[1,5,10,50,100,1000];
Dp = [0,1,10,100];
T_half_protein = [1,10,100,1000];


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
                Corr_Coeff(i,j,k,l) = calculate_correlate_prediction_measurement_Protein(Prefix,...
                                        params,nc14,nucfluo_mean_BGsubt,nucfluo_sem,tRes);
            end
        end
    end
end
          
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


%% Define a function to calculate the correlation coefficient from a set of parameters.                            
function [Corr_Coeff] = calculate_correlate_prediction_measurement_Protein(Prefix,...
                                params,nc14,nucfluo_mean_BGsubt,nucfluo_sem,tRes)   
%% Read out essential fields
% FilePath = 'S:\YangJoon\Dropbox\CentralDogmaResults';
% cn = load([FilePath,filesep,...
%             Prefix,filesep,'CompiledNuclei.mat']);
% nc13 = cn.nc13;
% nc14 = cn.nc14;
% 
% % redefine the nc14 for synchronization of datasets
% nc14 = nc14-nc13;
%% read parameters
Dm = params(1);
T_half_mRNA = params(2);
Dp = params(3);
T_half_protein = params(4);

params = [Dm, T_half_mRNA, Dp, T_half_protein];

% This step of defining the parameter should be done in a more systematic
% way. For example, for picking values one by one from n-dimensional matrix

[Time,AccumulatedmRNA, ErrorAccumulatedmRNA, ...
            Protein,ErrorProtein,startNC,dt,Time_interp] = ...
                                estimate_cytomRNA_protein(Prefix,'NC',13,'params',params);
%% Correlate the prediction and measurement (2D plot for a given set of parameters)
% First, pick values from a subset of time points (in the future, we might be
% able to get all time points, which would be helpful in calculating the
% Correlation Coefficient, but for now, for convenience, I'll just grab 4
% time window, T = 10, 20, 30, and 40 minutes into nc14.
% Alternatively, we can do 5, 15, 25, 35 minutes
% The reason I'm not taking into account of the 0 min is because it's
% confounded with the import of TFs at the beginning of nuclear cycle.
% 
timePoints = [5,15,25,35]; % for this specific dataset (Prefix, as I don't have the 40 minutes long).

framePoints = nc14 + round(timePoints/tRes); % time frames from the beginning of nc14 to the time window

% Extract the Predicted and Measured protein at selected time points (frame
% points)
PredictedProtein = Protein(framePoints,:);
PredictedProtein_error = ErrorProtein(framePoints,:);

MeasuredProtein = nucfluo_mean_BGsubt(framePoints,:);
MeasuredProtein_error = nucfluo_sem(framePoints,:);

%% Plotting module
% addpath('../utilities');
% 
% nDataPoints = 10;
% % color map
% hm_cm = flipud(brewermap(nDataPoints,'Spectral'));
% colormap(hm_cm);

%% (Optional) Plot the Prediction and Measurement over AP axis

% % Define the AP axis
% APaxis = 0:0.025:1;
% 
% % plot every 10 minutes from 5 minutes, so, 5, 15, 25, and 25 minutes.
% % framePoints
% 
% % First, the predicted protein
% predict_fig = figure(1)
% hold on
% for frame = 1:length(framePoints)
%     errorbar(APaxis, PredictedProtein(frame,:), PredictedProtein_error(frame,:))
% end
% xlim([0.2 1])
% ylim([0 max(max(PredictedProtein))+30000])
% xticks([0.2 0.4 0.6 0.8  1])
% xticklabels([20 40 60 80 100])
% % xticks([0.2 0.3 0.4 0.5 0.6])
% % xticklabels([20 30 40 50 60])
% title('protein-predicted')
% xlabel('embryo length (%)')
% ylabel('protein concentration(AU)')
% legend('5 min','15 min','25 min','35 min')
% StandardFigure(predict_fig,predict_fig.CurrentAxes)
% 
% figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
% saveas(predict_fig, [figPath, filesep, 'protein_predicted_eveCondition','.tif'])
% saveas(predict_fig, [figPath, filesep, 'protein_predicted_eveCondition','.pdf'])
% 
% % Second, the measured protein (This has to be revisited after a more
% % rigorous characterization of BG subtraction).
% % First, the predicted protein
% measure_fig = figure(2)
% hold on
% for frame = 1:length(framePoints)
%     errorbar(APaxis, MeasuredProtein(frame,:), MeasuredProtein_error(frame,:))
% end
% xlim([0.2 1])
% ylim([0 max(max(MeasuredProtein))+max(max(MeasuredProtein))/10])
% xticks([0.2 0.4 0.6 0.8  1])
% xticklabels([20 40 60 80 100])
% % xticks([0.2 0.3 0.4 0.5 0.6])
% % xticklabels([20 30 40 50 60])
% title('protein-measured (LlamaTag)')
% xlabel('embryo length (%)')
% ylabel('protein concentration(AU)')
% legend('5 min','15 min','25 min','35 min')
% StandardFigure(measure_fig,measure_fig.CurrentAxes)
% 
% figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
% saveas(measure_fig, [figPath, filesep, 'protein_measured','.tif'])
% saveas(measure_fig, [figPath, filesep, 'protein_measured','.pdf'])
%% Plot the Prediction and Measurement to see the correlation

% % first, color for time points
% % color map
% hm_cm = flipud(brewermap(4,'Spectral'));
% hold on
% for tpoint = 1:length(framePoints)
%     plot(PredictedProtein(tpoint,:), MeasuredProtein(tpoint,:), 'o','MarkerFaceColor',hm_cm(tpoint,:),'MarkerEdgeColor','black')
% end
% 
% % labels
% title('protein : prediction vs measurement')
% xlabel('prediction (AU)')
% ylabel('measurement (AU)')
% legend('5','15','25','35','Location','NorthWest')
% 
% % Axes Ticks
% xticks([0 0.4 0.8 1.2 1.6]*10^6)
% yticks([100 200 300 400 500])
% 
% StandardFigure(gcf,gca)
% 
% % Save the plot
% figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
% saveas(gcf, [figPath, filesep, 'correlation_prediction_measurement_timepoints_nc14_eveCondition','.tif'])
% saveas(gcf, [figPath, filesep, 'correlation_prediction_measurement_timepoints_nc14_eveCondition','.pdf'])

%% Calculate the correlation coefficient
% Pearson's correlation coefficient

% first, we have to get rid of NaNs.
% Here, we will trim the matrix for the APbins that has NaN, either in the
% prediction or in our measurement. Usually, the measurement has more NaNs,
% than the prediction.
sum(PredictedProtein);
nanIndex_Predicted = isnan(sum(PredictedProtein));

sum(MeasuredProtein);
nanIndex_measured = isnan(sum(MeasuredProtein));

nanIndex_sum = nanIndex_Predicted + nanIndex_measured;

% This is a filter for the columns that does not have NaNs
filter_index = find(nanIndex_sum==0);

% Trim the PredictedProtein and MeasuredProtein for the APbins with Nans
PredictedProtein_filtered = PredictedProtein(:,filter_index);
MeasuredProtein_filtered = MeasuredProtein(:,filter_index);

% Now, our Prediction and Measurements have two dimensions, (frame points)
% x (APbins), we will linearize it for an easier calculation as all time
% points should be aligned in a linear curve in an ideal case.

Prediction = reshape(PredictedProtein_filtered, [length(framePoints)*length(filter_index), 1]);
Measurement = reshape(MeasuredProtein_filtered, [length(framePoints)*length(filter_index), 1]);

% calculate the covariance
COV = cov(Prediction, Measurement);

STD_predicted = nanstd(Prediction);
STD_measured = nanstd(Measurement);

Corr_Coeff = COV(1,2)/(STD_predicted*STD_measured);
end
%% plot the correlation coefficient
%% 