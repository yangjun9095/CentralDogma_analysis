function main_CorrelateMS2Protein_APbins(DataType)
%% Description
% Script to analyze the cross-correlation of mRNA and Protein for each AP bin, 
% for a given Prefix (Maybe we'd ultimately want to expand this for
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
% OUTPUT: 

%% Load datasets
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

%% Define the fields for the instantaneous MS2 fluorescence
end