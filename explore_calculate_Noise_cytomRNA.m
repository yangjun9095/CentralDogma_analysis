function explore_calculate_Noise_cytomRNA(DataType)
% Description
% I'm curious to see the noise (variability) at the level of  cytoplasmic mRNA
% which we can quantify with MS2 technique.
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

%% Plot the totalmRNA (with Error) for individual nuclei
% Define NC
NC = 13;
% Define the minimal length (this is our QC filter for now)
minLength = 5;
maxLength = nc14-nc13; % only for nc13 as there are particles that are tracked all the way to nc14
% These particles should be curated later by either QC script or manual
% curation.

% Loop through all particles and plot the totalmRNA(with Error) over the AP axis
figure(2)
hold on
for k=1:length(particles)
    if particles(k).nc == NC &&... % at specific NC
            length(particles(k).Frame)>=minLength &&... & QC for length of frames
            length(particles(k).Frame)< maxLength &&...
            max(particles(k).Frame) < nc14 &&...% QC for nc tracking (previous/next NCs)
            min(particles(k).Frame) > nc13

    errorbar(particles(k).MeanAP,particles(k).TotalmRNA,particles(k).TotalmRNAError,'ob') 
    end
end

xlim([0.15 0.6])
xlabel('AP axis (EL)')
ylabel('accumulated mRNA (AU)')
title('')

% Save the plot
% info : NC, minLength, and Prefix?
mkdir(['S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures',filesep,'TotalmRNA']);
figPath =  'S:\YangJoon\Dropbox\Garcia Lab\Figures\CentralDogmaFigures\TotalmRNA';
saveas(gcf, [figPath, filesep, 'cytomRNA_indPtcle_APaxis_NC13_minLength=5_Prefix=',Prefix,'.tif'])
saveas(gcf, [figPath, filesep, 'cytomRNA_indPtcle_APaxis_NC13_minLength=5_Prefix=',Prefix,'.pdf'])


%% Second, let's calculate the nucleus-to-nucleus variability of total mRNA 
% (accumulated mRNA or cytoplasmic mRNA)

% Let's check nc12, nc13, and nc14
% The logic here is to calculate nucleus-to-nucleus variability of total
% mRNA pulled from the compiledparticles. There are many caveats in here.
% For example, we have to be careful with 1) tracking, 2) false-positive spots,
% and 3) particles that are tracked from the previous cycle, etc.

% These quality control should be checked carefully later.

% Check the nucleus-to-nucleus variability in nc13 and nc14
% Let's sort out the particles in nc13, with some minimal length
% (frame), for each AP bin.
NC=13;
minLength = 5; % minimal frame length
maxLength = nc14-nc13; % max length for NC13
%maxLength = length(Time) -nc14; % max length for NC13

% a structure to save the sorted particles
% We're saving the totalmRNA/totalmRNAError calculated in the
% compileparitlces
%sortedTotalmRNA = cell(2,41);
%sortedTotalmRNAError = cell(2,41);

for i=1:length(cp.APbinID)-1
    APbin1 = cp.APbinID(i);
    APbin2 = cp.APbinID(i+1);
    totalmRNATemp = [];
    totalmRNAErrorTemp = [];
    
    for k=1:length(particles)
        if particles(k).nc == NC &&... % at specific NC
                length(particles(k).Frame)>=minLength &&...% Check the length
                length(particles(k).Frame)< maxLength &&...
                particles(k).MeanAP > APbin1 &&... % AP filter
                particles(k).MeanAP <= APbin2&&...        
                min(particles(k).Frame) > nc13&&...% QC for nc tracking (previous/next NCs)
                max(particles(k).Frame) < nc14 
                
            
            % append the vector with totalmRNA and totalmRNAError
            totalmRNATemp = [totalmRNATemp, particles(k).TotalmRNA];
            totalmRNAErrorTemp = [totalmRNAErrorTemp, particles(k).TotalmRNAError];           
        end
    end
    sortedTotalmRNA{NC-12,i} = totalmRNATemp;
    sortedTotalmRNAError{NC-12,i} = totalmRNAErrorTemp;
end

%% Now, calculate the nucleus-to-nucleus (particle-to-particle variability)
% Use the sortedTotalmRNA struct
% Go through AP bins, and go through each NC
Mean_mRNA = nan(2,41);
STD_mRNA = nan(2,41);
CV_sqrd_mRNA = nan(2,41);

for i=1:length(sortedTotalmRNA)
    for j=1:2 % NCs
        totalmRNATemp = sortedTotalmRNA{j,i};
        if ~isnan(nanmean(totalmRNATemp))
            mean_temp = nanmean(totalmRNATemp);
            std_temp = nanstd(totalmRNATemp,0,2);
            CV_temp = std_temp./mean_temp;
            
            CV_sqrd_mRNA(j,i) = CV_temp.^2;
            Mean_mRNA(j,i) = mean_temp;
            STD_mRNA(j,i) = std_temp;
        end
    end
end

%% Plot the noise (CV^2) along the AP axis

hold on
plot(0:0.025:1, CV_sqrd_mRNA(1,:))
plot(0:0.025:1, CV_sqrd_mRNA(2,:))

xlabel('')
ylabel('')
title('')
legend('nc13','nc14')

%% Plot the noise (CV^2) along with the mean
hold on
plot(Mean_mRNA(1,:), CV_sqrd_mRNA(1,:))%,'o')
plot(Mean_mRNA(2,:), CV_sqrd_mRNA(2,:),'o')

xlabel('mean cytoplasmic mRNA (AU)')
ylabel('$$CV^2$$')
title('')
legend('nc13','nc14')

%% Part2.
% Calculate the noise between particles at each time point.

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