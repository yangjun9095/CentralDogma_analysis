%This code is for comparison of different MCP-mCherry lines
%For example, MCP-NLS-mCherry, MCP-NoNLS-mCherry, MCP-mCherry with
%different dosages (x2, x1.5)
% Note : For excitation of mCherry, I used 25uW, and for eGFP, I used 30uW
clear all
%Load the Datasets
DataNLS=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-01-02-Hb-P2P-MS2-NLSmCherry\CompiledParticles.mat')
DataNoNLS=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-01-04-Hb-P2P-MS2-mCherry\CompiledParticles.mat')
DataCTL=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-01-08-Hb-P2P-HisRFP-MCP-GFP\CompiledParticles.mat')

%Integrate spot fluorescence over time, to get accumulated mRNA

%% 1.NLS
NLSTime=DataNLS.ElapsedTime;
NLSFluo=DataNLS.MeanVectorAP;
NLSAPbinID=DataNLS.APbinID;
%plot(NLSTime,NLSSpot)
NLS13=DataNLS.nc13;
NLS14=DataNLS.nc14;

%Calculate the Accumulation of mRNA
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./60; %mRNA half-life = 60min
MRNA=NLSFluo;
MRNA(isnan(MRNA))=0;
Nt=length(NLSTime); %# of time points.
AccumulatedmRNA=zeros(Nt-4,41);
%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
%didn't consider diffusion in here.
t13=DataNLS.nc13;
t14=DataNLS.nc14;

tstart=t14;
tend=length(Nt);

for i=tstart:tend-4
    for j=1:41
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+(NLSTime(i)-NLSTime(i-1))*MRNA(i,j) - lambdam*AccumulatedmRNA(i-1,j)*(NLSTime(i)-NLSTime(i-1));
    end
end

%% Plot for NLS-MCP-mCherry

%Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(NLSTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
hold on
for i=2:Nt-4
    plot(0:0.025:1,AccumulatedmRNA(i,:),'color',Color(i-iStart+1,:))
end
hold off
title('Accumulated mRNA along the AP axis')
xlabel('AP axis')
ylabel('Accumulated mRNA (AU)')
%xlim([0.25 0.6])
%set(gca,'fontsize',30)

%% 2. No NLS
NoNLSTime=DataNoNLS.ElapsedTime;
NoNLSSpot=DataNoNLS.MeanVectorAP;

%plot(NLSTime,NLSSpot)
t13=DataNoNLS.nc13;
t14=DataNoNLS.nc14;

%Calculate the Accumulation of mRNA
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./60; %mRNA half-life = 60min
MRNA=NoNLSSpot;
MRNA(isnan(MRNA))=0;
Nt=length(NoNLSTime); %# of time points.
AccumulatedmRNA=zeros(Nt-4,41);
%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
%didn't consider diffusion in here.
tstart=t13
tend=t14
for i=tstart:tend-4
    for j=1:41
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+(NoNLSTime(i)-NoNLSTime(i-1))*MRNA(i,j) - lambdam*AccumulatedmRNA(i-1,j)*(NoNLSTime(i)-NoNLSTime(i-1));
    end
end

%% plot
%Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(NLSTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
hold on
for i=2:Nt-4
    plot(0:0.025:1,AccumulatedmRNA(i,:),'color',Color(i-iStart+1,:))
end
hold off
title('Accumulated mRNA along the AP axis')
xlabel('AP axis')
ylabel('Accumulated mRNA (AU)')
%xlim([0.25 0.6])
%set(gca,'fontsize',30)

%% 3. Control
CTLTime=DataCTL.ElapsedTime;
CTLSpot=DataCTL.MeanVectorAP;

%plot(NLSTime,NLSSpot)
t13=DataCTL.nc13;
t14=DataCTL.nc14;

%Calculate the Accumulation of mRNA
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./60; %mRNA half-life = 60min
MRNA=CTLSpot;
MRNA(isnan(MRNA))=0;
Nt=length(CTLTime); %# of time points.
AccumulatedmRNA=zeros(Nt-4,41);
%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
%didn't consider diffusion in here.
tstart=t13
tend=t14
for i=tstart:tend-4
    for j=1:41
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+(NoNLSTime(i)-NoNLSTime(i-1))*MRNA(i,j) - lambdam*AccumulatedmRNA(i-1,j)*(NoNLSTime(i)-NoNLSTime(i-1));
    end
end

%% plot
%Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(CTLTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
hold on
for i=2:Nt-4
    plot(0:0.025:1,AccumulatedmRNA(i,:),'color',Color(i-iStart+1,:))
end
hold off
title('Accumulated mRNA along the AP axis')
xlabel('AP axis')
ylabel('Accumulated mRNA (AU)')
%xlim([0.25 0.6])
%set(gca,'fontsize',30)

%% Comparison of AccumulatedmRNA
NLSmNRA=NLSmNRA/max(NLSmNRA);
NoNLSmNRA=NoNLSmNRA/max(NoNLSmNRA);
CTLmNRA=CTLmNRA/max(CTLmNRA);


hold on
plot(0:0.025:1,NLSmNRA,'r')
plot(0:0.025:1,NoNLSmNRA,'k')
plot(0:0.025:1,CTLmNRA,'g')

title('Comparison of AccumulatedmRNA')
xlabel('AP')
ylabel('Normalized AccumulatedmRNA(AU)')
legend('MCP-NLS-mCherry','MCP-NoNLS-mCherry','MCP-GFP')
hold off
saveas(gcf,'C:\Users\YangJoon\Desktop\Yang Joon\MCP Comparison-mRNA-nc13','tiff')

%% Offset Comparison (for different MCPs)

hold on
plot(DataNLS.ElapsedTime(DataNLS.nc13:end)-DataNLS.ElapsedTime(DataNLS.nc13),DataNLS.MeanOffsetVector(DataNLS.nc13:end),'r')
plot(DataNoNLS.ElapsedTime(DataNoNLS.nc13:end)-DataNoNLS.ElapsedTime(DataNoNLS.nc13),DataNoNLS.MeanOffsetVector(DataNoNLS.nc13:end),'b')
plot(DataCTL.ElapsedTime(DataCTL.nc13:end)-DataCTL.ElapsedTime(DataCTL.nc13),DataCTL.MeanOffsetVector(DataCTL.nc13:end),'g')
hold off

title('Offset for different MCPs')
xlabel('Time(min)')
ylabel('Mean Offset (AU)')
legend('MCP-NLS-mCherry(25uW)','MCP-NoNLS-mCherry(25uW)','MCP-GFP(30uW)')
saveas(gcf,'C:\Users\YangJoon\Desktop\Yang Joon\MCP Comparison','tiff')