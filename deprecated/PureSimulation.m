%Simulation for CentralDogma with simple prediction

%First, we assume the transcription rate to be step-function.
%Define step-funciton like transcription activity
for i=1:41
    if i<18
        Rate(i)=7000;
    elseif i>18
        Rate(i)=0;
    else Rate(i)=3500;
    end
end
APbin=[0:0.025:1];

%plot(APbin,Rate,'k','LineWidth',2)
plot(APbin,zeros(1,41),'k','LineWidth',2)
ylim([0 7000])
title('Transcription Rate','Fontsize',35)
xlabel('AP position','Fontsize',25)
ylabel('Transcription Rate','Fontsize',25)
% xt = get(gca, 'XTick');
% set(gca, 'FontSize', 16)
% yt = get(gca, 'YTick');
% set(gca, 'FontSize', 16)

%% Second, simulation for mRNA
clear M
clear T
%Let's start with considering nc14 only.
%First, I need to define variables and parameters.

%gamma_M : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
gammaM=log(2)/60; %Use Min. as time unit.
%r_m:Transcription rate for each AP bin, at nc14. from Hernan's 2013 Curr.Bio paper
rm=Rate; %(:,3);   %Normalized by # of RNAP
% %Remove Nans
% for i=1:41
%     rm(isnan(rm))=0;
% end
%Dm:diffusion constant of mRNA,
Dm=0.2*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=[1:41];
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
APbinpos = (AP-1)*dx;
km=Dm/((dx)^2); %1/Min. unit.

%Time window
T=[0:0.01:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.
dt=0.01;
%We need Transcription Time Window, along the AP axis, as in Fig.4B (1x41)
%Twindow=[0,0,0,0,0,0,0,0,18,17,17.5,17,19,17,15.5,16,15.5,12.5,11,11.5,10,8.5,11,8.5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6];

% %Fraction of Active Nuclei Fig.4.F
% FractionActiveNuclei=MaxOnRatio(:,3);
% FractionActiveNuclei(isnan(FractionActiveNuclei))=0;

%Effective rate of transcription along the AP axis. (1x41)
%rmeff=rm.*FractionActiveNuclei;

%Make a zero-matrix
M=zeros(length(T),length(AP));
M(1,:)=0;  %mRNA(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).

for i=2:length(T)  %for all Timepoints
    if i<length(T)*1/3
        M(i,1)=M(i-1,1)+rm(1)*dt-gammaM*M(i-1,1)*dt+km*M(i-1,2)*dt-km*M(i-1,1)*dt;

    else
        M(i,1)=M(i-1,1)-gammaM*M(i-1,1)*dt+km*M(i-1,2)*dt-km*M(i-1,1)*dt;

    end
end

for i=2:length(T)  %for all Timepoints
    if i<length(T)*1/3
        M(i,41)=M(i-1,41)+rm(41)*dt-gammaM*M(i-1,41)*dt+km*M(i-1,41)*dt-km*M(i-1,41)*dt;
    else
        M(i,41)=M(i-1,41)-gammaM*M(i-1,41)*dt+km*M(i-1,41)*dt-km*M(i-1,41)*dt;

    end
end


%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Make a mRNA reaction-diffusion equation
for j=2:length(AP)-1 %for all APbins
    for i=2:length(T)  %for all Timepoints
        
        if i<length(T)*1/3
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
end
hold off

title('mRNA along the AP axis')
xlabel('AP axis')
ylabel('mRNA (AU)')

%% Third, Protein simulation 
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

for i=2:length(T)  %for all Timepoints
    P(i,1)=P(i-1,1)+rp(1)*dt-gammaP*P(i-1,1)*dt+k*P(i-1,2)*dt-k*P(i-1,1)*dt;
    P(i,41)=P(i-1,41)+rp(1)*dt-gammaP*P(i-1,41)*dt+k*P(i-1,40)*dt-k*P(i-1,41)*dt;
end

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

%% Ploting all on the same graph (needs scaling)
% first, I want to plot time course change of Transcription rate, mRNA,
% Protein
%I need to define some discrete time points for plotting.
t=[0,10,20,30,40,50,60];
tpoint=t/dt+1;

for i=1:length(t)
    figure(i)
    hold on
    plot(APbin,Rate,'k','LineWidth',2)
    plot(APbin,M(tpoint(i),:),'r','LineWidth',2)
    plot(APbin,P(tpoint(i),:)/200,'b','LineWidth',2)
    hold off
    title('Transcription Rate/mRNA/Protein over AP position')
    xlabel('AP position')
    ylabel('Transcription Rate/mRNA/Protein')
    
end

%% Plotting mRNA for specific time point
t=[0,15,30,45,60];
tpoint=t/dt+1;

for i=1:length(tpoint)
    figure(i)
    plot(APbin,M(tpoint(i),:),'r','LineWidth',2)
    
    title('Accumulated mRNA','Fontsize',35)
    xlabel('AP position','Fontsize',25)
    ylabel('mRNA(AU)','Fontsize',25)
    ylim([0 70000])
end

%% Plotting Protien for specific time point
t=[0,15,30,45,60];
tpoint=t/dt+1;

for i=1:length(tpoint)
    figure(i)
    plot(APbin,P(tpoint(i),:),'b','LineWidth',2)
    
    title('Protein','Fontsize',35)
    xlabel('AP position','Fontsize',25)
    ylabel('Protein(AU)','Fontsize',25)
    ylim([0 2500000])
end

