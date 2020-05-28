% Goal : generate plots of averaged MS2 fluo over each nc (over AP axis)
% for nc12, nc13, and nc14.

function plt_timeAveraged_MS2Fluo(Prefix)
%% Load datasets
FilePath = 'S:\YangJoon\Dropbox\CentralDogmaResults';
Prefix = '2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1';
%Prefix = '2018-08-09-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1-3';
cp = load([FilePath,filesep,...
            Prefix,filesep,'CompiledParticles.mat']);

        %% Extract useful fields from the cp
% For now, we will think about nc13, and nc14. This can be (and
% should be) relaxed to include specific NCs

% time info
Time = cp.ElapsedTime;
nc12 = cp.nc12;
nc13 = cp.nc13;
nc14 = cp.nc14

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
EL = 500; % embryo length

%% time-averaging with certain windows (this could be fed into the options)

%% plot the time-averaged MS2 fluo over AP axis
%% Plot the averaged MS2 fluo over AP axis (for different time points)
% % plot
% frameLength=length(Time); %# of time points.
% % Color(gradation from blue to red) for each time point
% iStart=2; %Start time
% iEnd=length(Time); %End time (cp.nc14:end)
%     colormap(jet(256));
%     cmap=colormap ;
%     Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
% 
% figure(1)
% hold on
% for i=2:frameLength
%     errorbar(APaxis,spotfluo_mean(i,:),spotfluo_sem(i,:),'color',Color(i-iStart+1,:))
%     pause
%     xlim([0.2 0.6])
%     ylim([0 500])
% end
% hold off
% title('RNAP loading rate over AP')
% xlabel('AP axis')
% ylabel('averaged MS2 spot fluorescence')
% %ylabel('RNAP loading rate (number of RNAPs/ min)')
% 
% set(gca,'fontsize',30)
% colorbar
end
