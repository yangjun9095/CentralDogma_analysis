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
% file path
DropboxPath = 'S:\YangJoon\Dropbox\CentralDogmaResults';

% Load nucleus_struct.mat (compiled over multiple embryos)
%load(['S:\YangJoon\Dropbox\CentralDogmaResults',filesep,'nucleus_struct.mat']);
% Note. For these nucleus_struct.mat format, there's QC flag field which I
% can use to filter out false-positives.

% For now, we will pick hb P2P-MS2-NB datasets
Prefix = '2018-08-20-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1'

cp = load([DropboxPath, filesep, Prefix, filesep, 'CompiledParticles'])

particles = cp.CompiledParticles{1,1};

Time = cp.ElapsedTime;
nc12 = cp.nc12;
nc13 = cp.nc13;
nc14 = cp.nc14;

%% Loop through all frames, at each frame, loop through all particles
% 1) sort out the particles that passes the quality control filter, 
% 2) extract the fluorescence from those particles, save into a structure.
% For now, we will filter off particles that have less than nFrameFilter
nFrameFilter = 3;

% result : a structure having (Frames) x (APbins)
sortedFluo = cell(length(Time), length(cp.APbinID)); % sortedFluo{i,j};

for i=1:length(Time)
    frameN = i; % frame index
    for j=1:length(cp.APbinID)-1
        APpos1 = cp.APbinID(j); % left boundary
        APpos2 = cp.APbinID(j+1); % right boundary
        % APpos1 <= particles.MeanAP < APpos2
        
        % empty array to save the fluo (temp)
        FluoTemp = [];
        % go through all particles to sort out the ones passes filter
        % (quality, frame, and APbin)
        for k=1:length(particles)
            % Note. we can use qc_flag field for nucleus_struct file
            if length(particles(k).Frame) > nFrameFilter &&... % nFrameFilter 
                    sum(particles(k).Frame == frameN)&&...     % particles(j) contains that specific frame (frameN) 
                    ((particles(k).MeanAP>=APpos1) && (particles(k).MeanAP<APpos2))% AP position check
               % frame index in this particle
               frameIndex = find(particles(k).Frame == frameN);
               fluo = particles(k).Fluo(frameIndex);
               FluoTemp = [FluoTemp, fluo];

            end
        end
        % Save that FluoTemp in the struct format at corresponding
        % Frame x APbin
        sortedFluo{i,j} = FluoTemp;
            
    end
end
        

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