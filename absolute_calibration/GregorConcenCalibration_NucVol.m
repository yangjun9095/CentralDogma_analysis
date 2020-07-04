
% Gregor concentration calibration
APpos = [0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
GregorBcdConcen_nM = [53 40 23 14 12 6 4 2.5 1 0.5 0];  %Bcd gradient in nM

%1nM conversion to molecules per cubic micron (molecules per femtoliter)
oneNanoMolar = (10^(-6)) * (10^(-6))^3 * 6.022*(10^23); 

% Bcd gradient in molecules per cubic micron
GregorBcdConcen_molecPerCubeMicron = GregorBcdConcen_nM * oneNanoMolar

% 1 fL = 1 um^3 ~ excitation volume ~ 1 pixel
% We are integrating over 9 pixels to get the protein signal at the locus
% numPixels = 9;
% Gregor_totalBcdPerDisk = GregorBcdConcen_molecPerCubeMicron * numPixels

%% getting the number of Bcd molecules per nucleus
% 1 nucleus : a sphere with a diameter of 10um (This should be revisited as
% the diameter gets smaller as the cycles goes). Perhaps I can use the
% radius from the nuclear masks (schnitzcells)
% load(['S:\YangJoon\Dropbox\CentralDogmaResults\2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1',...
%             filesep,'2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1_lin.mat'])

% Load the Ellipses to read the radius of nuclei at each cycle
load(['S:\YangJoon\Dropbox\CentralDogmaResults\2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1',...
    filesep,'Ellipses.mat'])

% Note.
% the diameter of nucleus at each nuclear cycle is defined in
% getDefaultParameters.m which is called in getDiameters in putCircleOnNuclei(TrackNuclei).
%     {space_resolution, 'space resolution', 'spaceResolution', 'space', 's'}...% pixel size in micrometers. 
%     {7.26, 'diameter nc7', 'd7'}... % diameter of nuclei during nc7 in micrometers.
%     {7.26, 'diameter nc8', 'd8'}... % diameter of nuclei during nc8 in micrometers.
%     {7.26, 'diameter nc9', 'd9'}... % diameter of nuclei during nc9 in micrometers.
%     {7.26, 'diameter nc10', 'd10'}... % diameter of nuclei during nc10 in micrometers.
%     {6.16, 'diameter nc11', 'd11'}... % diameter of nuclei during nc11 in micrometers.
%     {5.72, 'diameter nc12', 'd12'}... % diameter of nuclei during nc12 in micrometers.
%     {4.84, 'diameter nc13', 'd13'}... % diameter of nuclei during nc13 in micrometers.
%     {3.96, 'diameter nc14', 'd14'}... % diameter of nuclei during nc14 in micrometers.
% These are converted to the number of pixels in getDiameter using the
% space_resolution.
% In short, the radius of nuclei is calculated in putCircleOnNuclei, then
% saved into 3rd and 4th columns of Ellipses.mat as unit of [pixels]

nucDiameter = 3.96; % [um] in nc14

nucVolume = 4/3*pi*nucDiameter^3; %[um^3]

Gregor_totalBcdPerNucleus = GregorBcdConcen_molecPerCubeMicron * nucVolume;

