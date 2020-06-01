%Simulation for reaction-diffusion model for mRNA and Protein in fruit fly
%embryos.

%First, load the Data from Garcia et al., 2013, Curr.Bio paper
load('DataForPlots.mat');
%% Load RateWeightV2 -> the rate for different AP position (at each n.c.) data that we want
Rate = RateWeight;

APbin=[0:0.025:1];

hold on
for i=1:3
    plot(APbin,Rate(:,i))
end
    hold off
title('Mean Pol II loading rate over AP axis','Fontsize',20)
xlabel('AP position (x/L)','Fontsize',15)
ylabel('Pol II loading rate (au/min)','Fontsize',15)
%% mRNA simulation
clear M
clear T
%Let's start with considering nc14 only.
%First, I need to define variables and parameters.

%gamma_M : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
gammaM=log(2)/60; %Use Min. as time unit.
%r_m:Transcription rate for each AP bin, at nc14. from Hernan's 2013 Curr.Bio paper
rm=Rate(:,3);   %Normalized by # of RNAP
%Remove Nans
for i=1:41
    rm(isnan(rm))=0;
end
%Dm:diffusion constant of mRNA,
Dm=0.1*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=[1:41];
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
APbinpos = (AP-1)*dx;
km=Dm/((dx)^2); %1/Min. unit.

%Time window
T=[0:0.01:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.
dt=0.01;
%We need Transcription Time Window, along the AP axis, as in Fig.4B (1x41)
Twindow=[0,0,0,0,0,0,0,0,18,17,17.5,17,19,17,15.5,16,15.5,12.5,11,11.5,10,8.5,11,8.5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6];

%Fraction of Active Nuclei Fig.4.F
FractionActiveNuclei=MaxOnRatio(:,3);
FractionActiveNuclei(isnan(FractionActiveNuclei))=0;

%Effective rate of transcription along the AP axis. (1x41)
rmeff=rm.*FractionActiveNuclei;

%Make a zero-matrix
M=zeros(length(T),length(AP));
M(1,:)=0;  %mRNA(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%M(i,1)=M(i-1,1)+rm(1)*dt-gammaM*M(i-1,1)*dt+k*M(i-1,2)*dt-k*M(i-1,1)*dt;
%M(i,41)=...
%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Make a mRNA reaction-diffusion equation
for j=2:length(AP)-1 %for all APbins
    for i=2:length(T)  %for all Timepoints
        if i<Twindow(j)
            M(i,j)=M(i-1,j)+rm(j)*dt-gammaM*M(i-1,j)*dt+km*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km*M(i-1,j)*dt;
        else
            M(i,j)=M(i-1,j)-gammaM*M(i-1,j)*dt+km*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km*M(i-1,j)*dt;
        end
    end
end

%colorful plot for timecourse : To see how the shape changes along the
%time. 
 jStart=1;
 jEnd=length(T)+1;
 jTotal=jEnd-jStart;

 colormap(jet(256));
 cmap=colormap;
 Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);
    
hold on
for l=1:length(T)
    plot([0:0.025:1],M(l,:),'color',Color(l-jStart+1,:));
%     ylim([0,1000])
%     pause(0.01)
     drawnow
end
hold off


%% Make a Movie (for Transcription Rate vs AP axis, along time)
%This plots mRNAlevel with respect to each AP bin over time.
%Each color means each time point.

% %Color(gradation from blue to red) for each time point
% %colorful plot for timecourse : To see how the shape changes along the
% %time. 
%  jStart=1;
%  jEnd=length(T)+1;
%  jTotal=jEnd-jStart;
% 
%  colormap(jet(256));
%  cmap=colormap;
%  Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);
    
Z=zeros(1,41);
%Make a Movie
writerObj = VideoWriter('/Users/yangjoonkim/Documents/MATLAB/TranscriptionRate-nc14.avi'); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj); 

for l=1:10:length(T)
    if l<2000
        plot([0:0.025:1],Rate,'Color',[ 0.75 0.75 0 ],'LineWidth',10);       %'color',Color(l-jStart+1,:));
    else
        plot([0:0.025:1],Z,'Color',[ 0.75 0.75 0 ],'LineWidth',10)
    end
    title('Transcription Rate along the AP axis')
    xlabel('AP Position')
    ylabel('Transcription Rate (AU)')
    ylim([0,8000])
    set(gca,'fontsize',20)
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end


close(writerObj);

%% Make a Movie (for mRNA vs AP axis, along time)
%This plots mRNAlevel with respect to each AP bin over time.
%Each color means each time point.


%Make a Movie
writerObj = VideoWriter('/Users/yangjoonkim/Documents/MATLAB/mRNA-nc14.avi'); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj); 

for l=1:10:length(T)
    
    plot([0:0.025:1],M(l,:),'r','LineWidth',10)
    title('mRNA along the AP axis')
    xlabel('AP Position')
    ylabel('mRNA (AU)')
    ylim([0,125000])
    set(gca,'fontsize',20)
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end


close(writerObj);

%% Make a Movie (for Protein vs AP axis, along time)
%This plots Protein level with respect to each AP bin over time.
%Each color means each time point.


%Make a Movie
writerObj = VideoWriter('/Users/yangjoonkim/Documents/MATLAB/Protein-nc14.avi'); % Name it.
writerObj.FrameRate = 10; % How many frames per second.
open(writerObj); 

for l=1:10:length(T)
    
    plot([0:0.025:1],P(l,:),'b','LineWidth',10)
    title('Protein along the AP axis')
    xlabel('AP Position')
    ylabel('Protein (AU)')
    ylim([0,6500000])
    set(gca,'fontsize',20)
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
end


close(writerObj);
%% Protein simulation
clear P
%Let's start with considering nc14 only.
%First, I need to define variables and parameters.
%gamma_P, r_p, Dp(kp)
%gamma_P : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
gammaP=log(2)/50; %Use Min. as time unit. I can change this later
%r_p:Translation rate for each AP bin, at nc14. 
rp=2; %I can change it later. from Petkova et al., 2014 Curr. Bio (Gregor lab)

%Dp:diffusion constant of Protein,
Dp=0.1*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=[1:41];
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
APbinpos = (AP-1)*dx;
kp=Dp/(dx)^2; %1/Min. unit.

%Make a zero-matrix
P=zeros(length(T),length(AP));
P(1,:)=0;  %Protein(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%P(i,1)=P(i-1,1)+rp(1)*dt-gammaP*P(i-1,1)*dt+k*P(i-1,2)*dt-k*P(i-1,1)*dt;
%P(i,41)=...
%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Make a Protein reaction-diffusion equation
%I will assume that the translation occurs during the whole nc14.
for i=2:length(T)  %for all Timepoints
    for j=2:length(AP)-1 %for all APbins
        P(i,j)=P(i-1,j)+rp*M(i-1,j)*dt-gammaP*P(i-1,j)*dt+kp*dt*(P(i-1,j-1)+P(i-1,j+1))-2*kp*P(i-1,j)*dt;
    end
end

%% colorful plot for Protein : To see how the shape changes along the time. 

hold on
for k=1:length(T)
    plot([0:0.025:1],P(k,:));
end
hold off

%Let's see how the slope at the boundary change according to different Dm
%Dm=[0.1:0.1:10] -> Fig.4 or 5
%For example, for certain time points, compare the slopes for different Dm
%values.

%% Analysis (mRNA)
%I need to analyze the width of the boundary
%To do that, I would need to fit a linear line, using the least square
%approximation.

%Parameters
%1. I need to pick one time point for mRNA or Protein
%For mRNA, maybe when mRNA has its maximum? and for protein as well?

t=20 %min

%plot mRNA along the AP axis at time t.
mRNA = M(t/dt+1,:);
%MAX = max(mRNA);
%mRNA = mRNA/MAX; %for normalization



% hold on
%plot ([0:0.025:1],mRNA,'o') % Normalized mRNA along the AP axis
%plot ([0:0.025:1-0.025],slope,'o') %slope
% hold off


hold on
%Spline fitting
dxx=0.001
xx=0:dxx:1;
yy=spline(APbin,mRNA,xx);
plot(APbin,mRNA,'o',xx,yy)
ylim([-1 400])
for i=1:length(xx)-1
    slope(i)=(yy(i+1)-yy(i)) / 0.001;
end

Minslope=min(slope);
pos=find(slope==Minslope);
xpos=xx(pos);
ypos=yy(pos);
LineFit=Minslope*(xx-xpos)+ypos;
plot(xx,LineFit)

title('mRNA along AP axis','fontsize',25)
xlabel('AP axis','fontsize',20)
ylabel('total mRNA (AU)','fontsize',20)
legend('mRNA')

barh(max(yy),1);
barh(0,1);

hold off

%Get the xpos for 0.5*Maximum
%find(yy==0.5*max(mRNA))
for i=1:length(xx)
    if yy(i)>0.5*max(mRNA)
        if yy(i+1)<0.5*max(mRNA)
            xhalf=i*dxx;
        end
    end
end
xhalf;
%Get the width
%First, get the xmax, xmin
for i=1:length(xx)
    if LineFit(i)>max(yy)
        if LineFit(i+1)<max(yy)
            xmax=i*dxx;
        end
    end
    
end

for i=1:length(xx)
    if LineFit(i)>min(yy)
        if LineFit(i+1)<min(yy)
            xmin=i*dxx;
        end
    end
    
end

Width = xmin-xmax;

%% Analysis (Protein)
%I need to analyze the width of the boundary
%To do that, I would need to fit a linear line, using the least square
%approximation.

%Parameters
%1. I need to pick one time point for mRNA or Protein
%For mRNA, maybe when mRNA has its maximum? and for protein as well?

t=20 %min

%plot mRNA along the AP axis at time t.
Protein = P(t/dt+1,:);
%MAX = max(mRNA);
%mRNA = mRNA/MAX; %for normalization



% hold on
%plot ([0:0.025:1],mRNA,'o') % Normalized mRNA along the AP axis
%plot ([0:0.025:1-0.025],slope,'o') %slope
% hold off

hold on
%Spline fitting
xx=0:0.001:1;
yy=spline(APbin,Protein,xx);
plot(APbin,Protein,'o',xx,yy)
ylim([-1 23000])
for i=1:length(xx)-1
    slope(i)=(yy(i+1)-yy(i)) / 0.001;
end

Minslope=min(slope);
pos=find(slope==Minslope);
xpos=xx(pos);
ypos=yy(pos);
LineFit=Minslope*(xx-xpos)+ypos;
plot(xx,LineFit)

title('Protein along AP axis','fontsize',25)
xlabel('AP axis','fontsize',20)
ylabel('Protein concentration(AU)','fontsize',20)
legend('Protein')

barh(max(yy),1);
barh(0,1);

hold off