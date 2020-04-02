% Variability vs AP vs Time Analysis
% This code is for the noise (variability) analysis at the level of
% transcription / translation.
% Yang Joon Kim, (12/1/2017~

% Load the dataset ( I start with two datasets :the MCP-GFP dataset for mRNA and
% Hb-nb-eGFP dataset for protein for now, since we don't know if we
% saturate MS2 loops with the MCP-mCherry.)
Data = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-03-02-Hb-P2P-MS2V5-MCP-GFP\InterredmRNAProtein.mat')
% Load the dataset

%% variability of RNAP Loading rate
RNAPLoadingRate3
% Use the SDVectorAP and NParticlesAP to get the Standard Error of Mean
% (SEM)

%% Variability of mRNA (accumulated) 
% Note that this result could be compared with Shawn Little's paper(Cell, 2013)
% To calculate the variability of mRNA (produced in single nucleus), I need
% to track single traces pretty accurately. mRNA1(t), mRNA2(t),...,etc.

% For now, I will pick some AP bins, with well-tracked MS2 spot traces
% Question : How should I think about the lineages? Do I need His-iRFP?
% Even then, do I need to assume that the mRNA is separated equally during
% mitosis, or exported from the nucleus during mitosis?

% Let's start from considering nc14 only.
%% (Test) Pick one AP bin, and check if the particle tracking is done well
% Note. Particles.mat and CompiledParticles.mat are different.
% I will start with CompiledParticles.
Particles = MS2.CompiledParticles;

% Definitely, I need to be careful for CheckParticleTracking.
% I should go back, and manually curate the missing spots.

Index =[]; % Get the indices of particles that satisfies conditions below.
for i=1:length(Particles)
    clf
    % First, let's start with nc 14
    if Particles(i).nc==14
        % Pick some AP bins, let's start with 0.2~0.3 AP bin
        if (Particles(i).MeanAP <0.225) && (Particles(i).MeanAP >0.2)
            if sum(Particles(i).FrameApproved)>5 % At least 5 frames should exist.
                plot(Particles(i).Frame,Particles(i).Fluo,'-o')
                Index = [Index i]
                pause
            end
        end
    end
end
%% Check the Particle's position
xpos =[];
ypos= [];
for i=1:length(Index)
    j=Index(i);
    xpos(i) = nanmean(Particles(j).xPos);
    ypos(i) = nanmean(Particles(j).yPos);
    
    plot(xpos,ypos,'o','MarkerSize',11)
    pause
end

%% TotalmRNA and TotalmRNAError @ end of nc 14
SortedParticles = Particles(Index);
TotalmRNA = [];
TotalmRNAError = [];

% Use the TotalmRNA calculated by CompileParticles
for i=1:length(SortedParticles)
    TotalmRNA = [TotalmRNA SortedParticles(i).TotalmRNA];
    TotalmRNAError = [TotalmRNAError SortedParticles(i).TotalmRNAError];
end


%% For each Particle, calculate the accumulated mRNA at each time point.

SortedParticles = Particles(Index);

% Define a cell to save all variables.
mRNA ={};
for i=1:length(SortedParticles)
    mRNA(i).Frame = SortedParticles(i).Frame;
    for j=2:length(SortedParticles(i).Frame)
        mRNA(i).TotalmRNA(j) = trapz(SortedParticles(i).Frame(1:j),SortedParticles(i).Fluo(1:j));
    end
end

%% plot to check the mRNA accumulation
% for i=1:length(SortedParticles)
%     plot(mRNA(i).Frame,mRNA(i).TotalmRNA,'-o')
%     pause
% end
%% Match the Accumulated mRNA with time points
% Generate the time frame series

%% Calculate the mean and SEM of TotalmRNA
% MeanTotalmRNA = mean(TotalmRNA);
% STDTotalmRNA = std(TotalmRNA);
% SEMTotalmRNA = STDTotalmRNA / sqrt(length(TotalmRNA));

%% variability of protein (inferred)

%% variability of protein (measured by NB signal)
% For this, I can use the SDVectorAP and NParticles in CompiledNuclei.mat
%% Compare the level of variabilities at different steps of Central Dogma