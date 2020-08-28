%% Define a function to calculate the correlation coefficient from a set of parameters.                            
function [Corr_Coeff, PredictedProtein, PredictedProtein_error,...
          MeasuredProtein, MeasuredProtein_error] = ...
                    calculate_PearsonCorr_prediction_measurement_Protein(Prefix,...
                                            params,nucfluo_mean_BGsubt,nucfluo_sem)                               
%% Description
% This function takes inputs of Prefix, and parameters for the
% Reaction-Diffusion model, then calculates the Pearson correlation
% coefficient between the prediction (R-D model) and measurement.
% For now, the measurement values (LlamaTag) are plugged into this function
% as inputs, because we have not settled down for how to subtract the
% Bakcground Fluo yet.
% In the future, this function should only take the Prefix and parameters
% as inputs.
% Written by : Yang Joon Kim
% Last updated : 06/11/2020
%% Read out essential fields for the time vector
FilePath = 'S:\YangJoon\Dropbox\CentralDogmaResults';
cn = load([FilePath,filesep,...
            Prefix,filesep,'CompiledNuclei.mat']);
nc13 = cn.nc13;
nc14 = cn.nc14;

% redefine the nc14 for synchronization of datasets
nc14 = nc14-nc13;

time = cn.ElapsedTime;
% time resolution
tRes = median(diff(time)); % min
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

%% clear up the naming of data fields for saving/output
%Corr_Coeff
end