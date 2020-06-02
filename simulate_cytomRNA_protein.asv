% projectName = 'CentralDogma'
% Last updated : 5/22/2020
% Description : For a given Prefix (single embryo), estimate the
% cytoplasmic mRNA, and protein based on a set of parameters (Dm, Dp,
% rp, GammaP)

% Inputs : 1) single embryo (Prefix), 2) DataType (LoadMS2Sets) or
% AveragedDatasets, 3) individual embryos (optional...)

% Outputs : 1) cytomRNA : Time, MeanVectorAP, SDVectorAP, SEVectorAP, NParticlesAP, etc.
%           2) protein :  MeanVectorAP, SDVectorAP, SEVectorAP, NParticlesAP, etc.
%           3) parameters used to make this prediction/estimation, also
%               with dt
%           4) Time, interpolated time, interpolated MS2 meanvectorap, etc.
% Save those into two different structures, cytomRNA_prediction_params,
% protein_prediciton_params

function [Time,AccumulatedmRNA, ErrorAccumulatedmRNA, ...
            Protein,ErrorProtein,startNC,dt,Time_interp] = ...
                                simulate_cytomRNA_protein(Prefix, varargin)
%% Options

% Default options
startNC = 13; % which NC to start calculating

% parameters for prediction/estimation
% KEEP THE UNITS as "um" and "minutes"!
Dm = 0; % diffusion constant (um^2/sec)
T_half_mRNA = 60 % min
Dp = 7 % diffusion constant (um^2/sec)
T_half_protein = 50 % min
params = [Dm, T_half_mRNA, Dp, T_half_protein];

for i=1:length(varargin)
    if strcmpi(varargin{i},'params')
        params = varargin{i+1}; % [Dm, T_half_mRNA, Dp, T_half_protein];
    elseif strcmpi(varargin{i},'NC')
        startNC=varargin{i+1};
    end
end

% Re-define the parameters after passing the option filter
Dm = params(1);
T_half_mRNA = params(2);
Dp = params(3);
T_half_protein = params(4);
%% Load datasets
FilePath = 'S:\YangJoon\Dropbox\CentralDogmaResults';
Prefix = '2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-09-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1-3';
cp = load([FilePath,filesep,...
            Prefix,filesep,'CompiledParticles.mat']);

%% Calibration of MS2 spot fluo with the # of RNAPs (to be done later)  
% Sketchy calibration by just comparing fluorescence (AU) with Hernan's
% 2013 paper. The AU maximum is ~1200, and the corresponding # of active
% % RNAP loaded is ~100. 
% % So, I will divide the fluorescence intensity (of my construct) with 12. 
% % This calibration should be done more rigorously later.
% spotfluo_mean = cp.MeanVectorAP(cp.nc14:end,:)/12;
% spotfluo_sd = cp.SDVectorAP(cp.nc14:end,:)/12;
% NParticles = cp.NParticlesAP(cp.nc14:end,:);
% 
% NBTime = cp.ElapsedTime(cp.nc14:end);

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

% AP axis
APaxis = cp.APbinID;
APbins = 1:length(APaxis);
nAPbins = length(APbins); % number of APbins
EL = 500; % embryo length



%% Convert the parameters with the same units for calculation (for both mRNA and protein)
% Parameters (um for length, and minutes for time)
% params = [Dm, T_half_mRNA, Dp, T_half_protein] ;
Dm = params(1);
T_half_mRNA = params(2);
Dp = params(3);
T_half_protein = params(4);
% Note that the units are minutes for the half-lives, and um^2/sec for
% diffusion coefficients.

% translation rate : Petkova,2014 on Bicoid
rp = 2; % [protein / (minute * mRNA)]
% Note that this could be measured using the translation reporter, or by
% careful calibration of our MS2 signal and nanobodies

% convert the units properly for the mRNA parameters
%T_half_mRNA = 60; %mRNA half-life = 60min, Little, 2013
GammaM=log(2)./T_half_mRNA; 
%Dm = 0;
Dm=Dm*60; % %um^2/min.
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
km=Dm/((dx)^2); %1/min. 

% convert the units properly for the protein parameters
%T_half_mRNA = 60; %mRNA half-life = 60min, Little, 2013
GammaP=log(2)./T_half_protein; 
%Dp = 7;
Dp=Dp*60; % %um^2/min.
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
kp=Dp/((dx)^2); %1/min.

%% Define the dt (time-step for numerical simulation).
% As we learned from PBoC classes, the jump rate k should meet the
% condition that dt << 1/k (=dx^2/D) to capture the jump event reliably.
% This means that dt is affected by the D as we fixed the "dx".
% Thus, we have to fix the maximum value of D. For now, I'll assume that
% the Dmax = 100 [um^2/sec] as this is order of magnitude larger than other
% proteins, freely diffusing eGFP protein (~25um^2/sec). I can also just
% the Dextrin molecule's diffusion coefficients as a reference for this
% justification.

Dmax = 100; %[um^2/sec]
Dmax = Dmax * 60; %[um^2/min]
k_max=Dmax/((dx)^2); %1/min
dt_ub = 1/k_max;
dt = 0.01; %[min], half of dt_ub (upper bound)

%% Interpolate the vectors with "dt" time interval
% For accumulated mRNA, we're taking the spot fluo 2 frames before, thus we
%'d need to start from the 3rd frame at least.

% initial time point
if startNC==12 && nc12>3
    startFrame = nc12;
elseif startNC==13 && nc13>3
    startFrame = nc13;
elseif startNC==14 && nc14>3
    startFrame = nc14;
else 
    startFrame = 3; % 3rd frame
end

% Trim down the vectors from the startFrame to the end frame
% Note that all the vectors now start from the 1st frame as we defined the
% startFrame. The rest of the useless info are omitted.
Time_trimmed = Time(startFrame:end) - Time(startFrame); % initialize the time points to 0
spotfluo_mean_trimmed = spotfluo_mean(startFrame:end,:);
spotfluo_sd_trimmed = spotfluo_sd(startFrame:end,:);
num_particles_trimmed = num_particles(startFrame:end,:);

Time_interp = Time_trimmed(1):dt:Time_trimmed(end);

% interpolate the 2D vectors for all AP bins
for i=1:nAPbins
    spotfluo_mean_interp(:,i) = interp1(Time_trimmed, spotfluo_mean_trimmed(:,i), Time_interp);
    spotfluo_sd_interp(:,i) = interp1(Time_trimmed, spotfluo_sd_trimmed(:,i), Time_interp);
    num_particles_interp(:,i) = interp1(Time_trimmed, num_particles_trimmed(:,i), Time_interp);
end


%% Calculate the accumulated mRNA using numerical simulation

% Convert NaNs to zeros for calculation (which we might need to convert
% back again at the end).
% MeanVectorAP
spotfluo_mean_interp(isnan(spotfluo_mean_interp))=0;
spotfluo_sd_interp(isnan(spotfluo_sd_interp))=0;
% number of MS2 particles per frame, AP bin
num_particles_interp(isnan(num_particles_interp)) = 0;

frameLength=length(Time_interp); %# of time frames.
AccumulatedmRNA=zeros(frameLength,nAPbins);
ErrorAccumulatedmRNA = zeros(frameLength,nAPbins);

%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
APstart = 20; % APbin start : 20%
APend = 60; % APbin end : 60%

APbinStart = APstart/2.5 + 1;
APbinEnd = APend/2.5 + 1;

tDelay = 1; % min
tDelayIndex = tDelay/dt;
% calculate the accumulated mRNA and Error
for i=2+tDelayIndex:frameLength % Integration starts from the 3rd time point (as we pull the MS2 fluo from 1 min before) % Bothma, 2014 paper supple.
    % For AP bins in the middle
    for j=APbinStart+1:nAPbins-1
        % master equation
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+... % mRNA from the previous time point
            dt*spotfluo_mean_interp(i-1-tDelayIndex,j).*num_particles_interp(i-1-tDelayIndex,j)-... % mRNA synthesis, Et=2min, thus dNp/dt ~Fluo(t-Et/2)=Fluo(t-1min)
            GammaM*AccumulatedmRNA(i-1,j)*dt+...% mRNA degradation
            km*dt*(AccumulatedmRNA(i-1,j-1)+AccumulatedmRNA(i-1,j+1))-2*km*AccumulatedmRNA(i-1,j)*dt; %mRNA diffusion
        ErrorAccumulatedmRNA(i,j) = sqrt((ErrorAccumulatedmRNA(i-1,j).^2 +...
                                            spotfluo_sd_interp(i-1-tDelayIndex,j).^2.*num_particles_interp(i--1-tDelayIndex,j) * dt));
    end
    
    % For the last AP bin(nAPbins), we assume that the diffusion to and from the
    % right is zero.
        AccumulatedmRNA(i,nAPbins) = AccumulatedmRNA(i-1,nAPbins)+...% mRNA from the previous time point
            dt*spotfluo_mean_interp(i-1-tDelayIndex,nAPbins).*num_particles_interp(i-1-tDelayIndex,j)-... % mRNA synthesis, Et=2min, thus dNp/dt ~Fluo(t-Et/2)=Fluo(t-1min)
            GammaM*AccumulatedmRNA(i-1,nAPbins)*dt+...% mRNA degradation
            km*dt*(AccumulatedmRNA(i-1,nAPbins-1))-km*dt*AccumulatedmRNA(i-1,nAPbins);
        ErrorAccumulatedmRNA(i,nAPbins) = sqrt(ErrorAccumulatedmRNA(i-1,nAPbins).^2 +...
                                                (spotfluo_sd_interp(i-1-tDelayIndex,nAPbins).^2.*num_particles_interp(i-1-tDelayIndex,nAPbins) * dt));

    % For the first AP bin(APbinStart), we assume that the diffusion to and from the
    % left is same.
        AccumulatedmRNA(i,APbinStart) = AccumulatedmRNA(i-1,APbinStart)+... % mRNA from the previous time point
            dt*spotfluo_mean_interp(i-1-tDelayIndex,APbinStart).*num_particles_interp(i-1-tDelayIndex,j)-... % mRNA synthesis
            GammaM*AccumulatedmRNA(i-1,APbinStart)*dt+...% mRNA degradation
            km*dt*(AccumulatedmRNA(i-1,APbinStart+1)-AccumulatedmRNA(i-1,APbinStart));
        ErrorAccumulatedmRNA(i,APbinStart) = sqrt(ErrorAccumulatedmRNA(i-1,APbinStart).^2 + ...
                                                    (spotfluo_sd_interp(i-1-tDelayIndex,APbinStart).^2 .*num_particles_interp(i-1-tDelayIndex,APbinStart) * dt));
    % For the last AP bin, we assume that the diffusion from the right is
    % zero (so, I included APbinEnd into for for loop above)
    
end

%% plot the accumulated mRNA profile over AP
% Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(Time_interp); %End time (cp.nc14:end)
    colormap(jet(256));
    cmap=colormap ;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

figure(1)
hold on
for i=2:100:frameLength
    errorbar(0:0.025:1,AccumulatedmRNA(i,:),...
                sqrt(ErrorAccumulatedmRNA(i,:).^2),...
                'color',Color(i-iStart+1,:))
     %pause
end
hold off
title('Accumulated mRNA along the AP axis')
xlabel('AP axis')
ylabel('Accumulated mRNA (Number of mRNA)')
xlim([0.2 0.6])
set(gca,'fontsize',30)
%colorbar

%% Calculate the Protein using numerical simulation
% Make sure to check the dimension of AccumulatedmRNA (depending on
% which interpolation I did)

% Note that I didn't consider the delay between the transcription and
% translation. Transcripts should be shuttled to the cytoplasm, then
% translated by the ribosome. Thus, there could be delay, which I'm not
% sure how to address now.

clear Protein
clear ErrorProtein

Protein=zeros(frameLength,nAPbins);
ErrorProtein = zeros(frameLength,nAPbins);

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%P(i,1)=P(i-1,1)+rp(1)*dt-GammaP*P(i-1,1)*dt+k*P(i-1,2)*dt-k*P(i-1,1)*dt;
%P(i,41)=...
%But, right now, we don't have to worry too much, because we don't care
%about the boundary now.

% Make a Protein reaction-diffusion equation
% I will assume that the translation occurs during the whole time( nc13 and
% nc14). From the Dufourt and Lagha paper, this seems to be a fair assumption.
% I will use the same argument as I calculated the mRNA,
for i=startFrame:frameLength  %for all Timepoints
    for j=APbinStart+1:nAPbins-1 %for all APbins
        
        Protein(i,j)=Protein(i-1,j)+...
            rp*AccumulatedmRNA(i-1,j)*dt-...
            GammaP*Protein(i-1,j)*dt+...
            kp*dt*(Protein(i-1,j-1)+Protein(i-1,j+1))-2*kp*Protein(i-1,j)*dt;
        ErrorProtein(i,j) = sqrt(ErrorProtein(i-1,j).^2 + rp*ErrorAccumulatedmRNA(i-1,j).^2*dt);
    end
    % the last AP bin (nAPbins) : no jump to and from the right
        Protein(i,nAPbins)=Protein(i-1,nAPbins)+...
            rp*AccumulatedmRNA(i-1,nAPbins)*dt-...
            GammaP*Protein(i-1,nAPbins)*dt+...
            kp*dt*Protein(i-1,nAPbins-1)-kp*Protein(i-1,nAPbins)*dt;
        ErrorProtein(i,nAPbins) = sqrt(ErrorProtein(i-1,nAPbins).^2 +...
                            rp*ErrorAccumulatedmRNA(i-1,nAPbins).^2*dt);
    % APbinStart bin, assuming the influx from the left side and outflux
    % from the APbinStart bin to the left is the same
        Protein(i,APbinStart)=Protein(i-1,APbinStart)+...
        rp*AccumulatedmRNA(i-1,APbinStart)*dt-...
        GammaP*Protein(i-1,APbinStart)*dt+...
        kp*dt*Protein(i-1,APbinStart+1)-kp*dt*Protein(i-1,APbinStart);

        ErrorProtein(i,APbinStart) = sqrt(ErrorProtein(i-1,APbinStart).^2 +...
                            rp*ErrorAccumulatedmRNA(i-1,APbinStart).^2*dt);
end

figure(2)
hold on
for i=3:100:frameLength
    errorbar(0:0.025:1,Protein(i,:),ErrorProtein(i,:),'color',Color(i-iStart+1,:))
end
hold off
title('Predicted protein along the AP axis')
xlabel('AP axis')
ylabel('Predicted Protein (Number of protein molecules)')
xlim([0.2 0.6])
set(gca,'fontsize',30)
colorbar

%% Plot the cyto mRNA and protein for 
%% cytomRNA - prediction
APbin = 11; % 20%
% 0.66 min is our temporal resolution, so we will just plot roughly similar
% time points
tRes = median(diff(Time));
tSteps = floor(tRes/dt); 

tWindow = 1:tSteps:length(Time_interp);
errorbar(Time_interp(tWindow),...
        AccumulatedmRNA(tWindow,APbin),...
        ErrorAccumulatedmRNA(tWindow,APbin),'r')
title('Accumulated mRNA - 25% EL')
xlabel('time into nc 13(min)')
ylabel('Accumulated mRNA (AU)')
StandardFigure(gcf,gca)
%% protein - prediction
APbin = 11; % 20%
tWindow = 1:tSteps:length(Time_interp);

errorbar(Time_interp(tWindow),...
        Protein(tWindow,APbin),...
        ErrorProtein(tWindow,APbin))
title(' protein-prediction- 25% EL')
xlabel('time into nc 13(min)')
ylabel('protein (AU)')
StandardFigure(gcf,gca)

%% Clean up the variables to save for further analysis
Time = Time_trimmed;
% nc13 and nc14 should be adjusted with ncStart

% savedVariables = {};
% savedVariables = [savedVariables,'Time','startNC','dt','Time_interp',...
%                     'AccumulatedmRNA', 'ErrorAccumulatedmRNA', ...
%                     'Protein','ErrorProtein','params'];
%                
% DropboxFolder = 'S:\YangJoon\Dropbox\CentralDogmaResults';                
% save([DropboxFolder,filesep,'cytomRNA_protein_estimated.mat'],'savedVariables', '-v7.3', '-nocompression')

                
end