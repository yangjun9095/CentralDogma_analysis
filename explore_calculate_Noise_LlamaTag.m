function explore_calculate_Noise_MS2(DataType)
% Description
% I'm curious to see the noise (variability) at the level of  instantaneous
% Txn which we monitor with MS2 technique.
%
% Input : 
% DataType (potentially use somewhat processed structure of particles that
% are synchronized already)
% 
% Output : Variability at each time point
% 
% main script 
% Notes .
% 1) which NC to look at? Let's start with NC13, and NC14 for now.
% 2) Some quality control steps are needed meaning that if there are
% particlees with less than Nf number of frames, we cut off those as it's
% likely that those are actually false-positives.
% 3) This whole analysis should be revisited for multiple datasets, also
% after careful curation.
% 4) As the CV^2 is a unitless metric, it'd be interesting to see this for
% MCP-GFP vs MCP-mCherry.

%% Load datasets
FilePath = 'S:\YangJoon\Dropbox\CentralDogmaResults';
Prefix = '2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-09-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1-3';
cn = load([FilePath,filesep,...
            Prefix,filesep,'CompiledNuclei.mat']);

%% Define some options/parameters
startNC = 13;
%% Extract useful fields
% mean nuclear fluo, sd, number of nuclei, time
% time info
Time = cn.ElapsedTime;
tLength = length(cn.MeanVectorAP(:,1));
% nc12 = cn.nc12;
% nc13 = cn.nc13;
% nc14 = cn.nc14;

% mean, sd, number of nuclei
nucfluo_mean = cn.MeanVectorAP;
nucfluo_sd = cn.SDVectorAP;
num_nuclei = cn.NParticlesAP;

nucfluo_sem = nucfluo_sd./num_nuclei;

%% generate plots for check

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

%% generate plots for check (before BG subtraction)
APaxis = 0:0.025:1;
tStep = floor(10/0.6);
hold on
for t=[nc14:tStep:tLength,tLength]
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
legend('0 min','10 min','20 min','30 min','35 min')
StandardFigure(gcf,gca)

%% generate plots for check (after BG subtraction)
% hold on
% for t=[nc14:tStep:tLength,tLength]
%     errorbar(APaxis, nucfluo_mean_BGsubt(t,:), nucfluo_sem(t,:))
% end
% xlim([0.2 1])
% ylim([0 max(nucfluo_mean_BGsubt(tLength,:))+100])
% xticks([0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
% xticklabels([20 30 40 50 60 70 80 90 100])
% % xticks([0.2 0.3 0.4 0.5 0.6])
% % xticklabels([20 30 40 50 60])
% title('protein-measured(LlamaTag)')
% xlabel('embryo length (%)')
% ylabel('protein concentration(AU)')
% legend('0 min','10 min','20 min','30 min','35 min')
% StandardFigure(gcf,gca)

%% Calculate the CV^2 for each frame and APbin
% Use the sortedFluo processed above (Frame x APbin)
% caveat : The number of particles could be something we have to be
% careful. If the number of particles are too small, then we might
% undersample that APbin and that time frame...This is something to think
% about.

% For now, let's go through each time frame and APbin to calculate the mean
% and std of the fluo.
Mean_fluo = nan(length(Time), length(cp.APbinID));
STD_fluo = nan(length(Time), length(cp.APbinID));
NParticles_fluo = nan(length(Time), length(cp.APbinID));
CV_fluo = nan(length(Time), length(cp.APbinID));
CV_sqrd_fluo = nan(length(Time), length(cp.APbinID));

for i=1:length(Time)
    for j=1:length(cp.APbinID)
        % Extract the fluorescence of all particles in that specific frame
        % and APbin (from sortedFluo)
        fluo = sortedFluo{i,j};
        
        if ~isempty(fluo)
            % calculate the mean, std, number of particles, etc.
            mean = nanmean(fluo);
            nParticles = length(fluo);
            std = nanstd(fluo,0,2);
            cv = std./mean;
            cv_sqrd = cv.^2;
        
            % Save Mean, STD, nParticles, CV, CV^2 into arrays.
            Mean_fluo(i,j) = mean;
            STD_fluo(i,j) = std;
            NParticles_fluo(i,j) = nParticles;
            CV_fluo(i,j) = cv;
            CV_sqrd_fluo(i,j) = cv_sqrd;
        end
    end
end

%% Visualize the CV^2 over time and AP axis

%% First, let's see how the CV^2 changes over time (this is basically a
% recap of Shawn Little 2013 analysis).
% Let's pick an AP bin, say in the anterior
APbin = 11;
plot(Time, CV_sqrd_fluo(:,APbin))

%% Second, CV^2 over AP axis
% at which time point?
for i=1:length(Time)
    plot(0:0.025:1, CV_fluo(i,:))
    ylim([0 1])
    pause
    % Put a legend with the time point (either in nc12,13,or 14, with min)
    
end

% It does seem that the CV^2 is mostly within [0,0.2] range.
% Note that the posterior positions might need another filter for the
% particle number ,etc.

%% Exploratory
% Plot the CV^2 versus mean
figure(2)
nThresh = 5;
hold on
for j=1:length(cp.APbinID)
    for i=1:length(Time)
        if NParticles_fluo(i,j) > 5
            plot(Mean_fluo(i,j), CV_sqrd_fluo(i,j),'o')
            pause
        end
    end
end

%% Reshape the vectors to make 1D for gramm plotting
% Plot from nc12
tBegin = nc12;

% Strip out the unnecessary time points
Time = Time(tBegin:end);
Mean_fluo = Mean_fluo(tBegin:end,:);
CV_fluo = CV_fluo(tBegin:end,:);
CV_sqrd_fluo = CV_sqrd_fluo(tBegin:end,:);
NParticles_fluo = NParticles_fluo(tBegin:end,:);

% Define APbin vector for color
APbinID = zeros(length(Time),length(cp.APbinID));
for i=1:length(cp.APbinID)
    APbinID(:,i) = i;
end

% Define Time vector for color
tRes = median(diff(cp.ElapsedTime));
Time_assign = zeros(length(Time),length(cp.APbinID));
for i=1:length(Time)
    Time_assign(i,:) = i*tRes;
end

mean_fluo = reshape(Mean_fluo,[length(Time)*length(cp.APbinID),1]);
cv_sqrd_fluo = reshape(CV_sqrd_fluo,[length(Time)*length(cp.APbinID),1]);
apbinID = reshape(APbinID,[length(Time)*length(cp.APbinID),1]);
time = reshape(Time_assign,[length(Time)*length(cp.APbinID),1]);
nParticles = reshape(NParticles_fluo,[length(Time)*length(cp.APbinID),1]);
%% Use the gramm for plotting
clear g
g(1,1)=gramm('x',mean_fluo,'y',cv_sqrd_fluo,'subset',nParticles>nThresh);
g(1,1).geom_point();
g(1,1).set_names('x','mean (AU)','y','CV^2');
g(1,1).set_title('No groups');


g(1,2)=gramm('x',mean_fluo,'y',cv_sqrd_fluo,'subset',nParticles>nThresh,'color',apbinID);
g(1,2).geom_point();
g(1,2).set_names('x','mean (AU)','y','CV^2','color','APbin');
g(1,2).set_title('AP bin');

g(1,3)=gramm('x',mean_fluo,'y',cv_sqrd_fluo,'subset',nParticles>nThresh,'color',time);
g(1,3).geom_point();
g(1,3).set_names('x','mean (AU)','y','CV^2','color','Time');
g(1,3).set_title('Time (min)');

%figure('Position',[100 100 800 800]);
g.draw();

% Save the plot
figPath = 'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures';
saveas(gcf, [figPath, filesep, 'Noise_MS2_instant','.tif'])
saveas(gcf, [figPath, filesep, 'Noise_MS2_instant','.pdf'])
end