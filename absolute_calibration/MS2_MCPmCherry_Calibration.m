%% MS2-MCP-mCherry calibration to absolute number of RNAPs
% Use HG's 2013 paper's smFISH calibration method
% The detailed description coudl be found on OneNote 
% (HGLab-Protocols > Misc.Protocols > Absolute Calibration of MS2 signals)
% or Lammers, et al., PNAS, 2020 supplementary section

%% smFISH data
% FISH absolute count: This is the number of mRNA molecules per nucleus produced in nc13 using single molecule FISH. 
% It is obtained by measuring the number of mRNA molecules at mitosis 13 and subtracting it from the number of mRNA molecules in nc12. 
% This is the mean for all nuclei anterior of 0.35AP. The calculation of this value can be found in LacZFISHDepth.m in the cell 
% "Figure S3G, S3H (deep sampling with manually selected regions)". Note that this value can also be called MeanOnRegion. 
% The factor of two accounts for the fact that we had two copies of the P2P in the flies we used for FISH.

% From the Lammers, 2020
OnRegionFISHFit = 220; %438.5733/2;
SEOnRegionFISHFit = 30; %55.2832/2;
SDOnRegionFISHFit = 135.1962/2;


%% MS2-MCP (mCherry) data
% DataStatus.xlsx tab : "MCPmCh_Calibration"
MS2Data = LoadMS2Sets('MCPmCh_Calibration');

% Note that the MCP-mCherry datasets with P2P.V1 does not have a
% well-defined APDivision.mat. Thus, I went around to manually define the
% NCs using the NCs in CompiledParticles.mat this could be revisited if
% needed. That version is IntegratemRNA_MCPmCh

[TotalProd,TotalProdError,TotalProdN,...
    MeanTotalProd,SDTotalProd,SETotalProd]=IntegratemRNA_MCPmCh(MS2Data,2,2);

%% Calculate the average in the anterior bins (20-30% of EL)
numEmbryos = length(MS2Data);

% Pick the AP bins anterior to 30% of the embryo length
APrange = 0.2:0.025:0.3;
APbinRange = int8(APrange/0.025) + 1;

MeanOnRegion13 = nanmean(MeanTotalProd(APbinRange,13));
SDOnRegion13 = nanmean(SDTotalProd(APbinRange,13));
SEOnRegion13 = SDOnRegion13./sqrt(numEmbryos);

%% Calculate the fluorescence of 1 RNAP using the equations in the OneNote
L1 = 1.275; % [kb], MS2 seq. length
L2 = 4.407; % [kb], from the end of MS2 seq. to the end of 3' sequence
Velong = 1.5; % [kb/min]

F_RNAP = Velong*MeanOnRegion13/OnRegionFISHFit*(1/(L1 + L2)); % [AU/RNAP]

F_RNAP_error = Velong*SEOnRegion13/OnRegionFISHFit*(1/(L1 + L2)); %[AU/RNAP]

laserPower = 25; % [uW]

%% Laser power correction
% Since we used 20uW for our data acquisition (all the other imaging
% condition are exactly the same)

%% save the result

% filePath
filePath = 'S:\YangJoon\Dropbox\CentralDogmaResults\MCPmCh_Absolute_Calibration';

% MS2MCPmCherry_RNAP_Calibration
save([filePath,filesep,'MS2MCPmCherry_RNAP_Calibration.mat'],...
                    'OnRegionFISHFit', 'SEOnRegionFISHFit', 'SDOnRegionFISHFit',...
                    'MeanOnRegion13', 'SDOnRegion13','SEOnRegion13',...
                    'L1','L2','Velong', 'F_RNAP', 'F_RNAP_error', 'laserPower', '-v7.3')

%% check the saved file
abs_calibration_result = load([filePath, filesep, 'MS2MCPmCherry_RNAP_Calibration.mat'])