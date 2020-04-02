%Analysis code for the CompiledParticles.mat and CompiledNuclei
%% Load the data set
clear all

mRNA = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-03-02-Hb-P2P-MS2V5-MCP-GFP\CompiledParticles.mat')
%Protein = load('CompiledNuclei.mat')

%BG=load('BG-CompiledNuclei')

% mRNA=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2016-11-23-Hb-nbGFP-MS2-NLSmCherry\CompiledParticles.mat');
% Protein=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2016-11-23-Hb-nbGFP-MS2-NLSmCherry\CompiledNuclei.mat');

% mRNA2=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2016-10-18-Hb-nbGFP-MS2-mCherry\CompiledParticles.mat');
% Protein2=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2016-10-18-Hb-nbGFP-MS2-mCherry\CompiledNuclei.mat');

%% Plot mRNA along the timecourse (for specific AP bin)

ap=12; % choose specific AP bin

%Use MeanVectorAP at APbinID = ap, and plot them for timecourse 
hold on
plot(mRNA.ElapsedTime,mRNA.MeanVectorAP(:,ap),'r')
%plot(mRNA.ElapsedTime(mRNA.nc13:end)-(mRNA.nc13-mRNA.nc13),mRNA.MeanVectorAll(mRNA.nc13:end),'b')
hold off

title('Transcription Rate over time')
xlabel('Time(min)')
ylabel('Transcription Rate')
legend('MCP-NLS-mcherry','MCP-NoNLS-mCherry')
%I need to sync the nc, and normalized them...
%For synchronization, use nc data for frame shifting
%How to normalize? Let's first use the intensity at nc13 as a standard

% nc12(1) = mRNA.nc12;
% nc12(2) = mRNA2.nc12;
% 
% nc13(1) = mRNA.nc13;
% nc13(2) = mRNA2.nc13;
% 
% nc14(1) = mRNA.nc14;
% nc14(2) = mRNA2.nc14;

%% Plot the Offset
hold on
%plot(mRNA.ElapsedTime,mRNA.MeanVectorAll,'b')
plot(mRNA.ElapsedTime,mRNA.MeanOffsetVector,'r')
hold off

title('Spot Fluorescence and Offset')
legend('MS2 spot fluorscence')%,'Offset')
xlabel('Time(min)')
ylabel('Fluorescence (AU)')
set(gca,'fontsize',30)

%% mRNA Accumulation from Spot fluorescence
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./60; %mRNA half-life = 60min
MRNA=mRNA.MeanVectorAP;%.*mRNA.OnRatioAP;
MRNA(isnan(MRNA))=0;
Nt=length(mRNA.ElapsedTime); %# of time points.
AccumulatedmRNA=zeros(Nt-4,41);
%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
%didn't consider diffusion in here.
for i=2:Nt-4
    for j=1:41
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+(mRNA.ElapsedTime(i)-mRNA.ElapsedTime(i-1))*MRNA(i,j) - lambdam*AccumulatedmRNA(i-1,j)*(mRNA.ElapsedTime(i)-mRNA.ElapsedTime(i-1));
    end
end

% plot
%Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(mRNA.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap ;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
hold on
for i=2:Nt-4
    plot(0:0.025:1,AccumulatedmRNA(i,:),'color',Color(i-iStart+1,:))
end
hold off
title('Accumulated mRNA along the AP axis')
xlabel('AP axis')
ylabel('Accumulated mRNA (AU)')
xlim([0.2 0.6])
set(gca,'fontsize',30)
colorbar

%% mRNA Smoothening (to get the boundary position and width)
%Basically, it's averaging AccumulatedmRNA (AP) over 5 adjacent AP bins
clear SmoothmRNA

%range of AP bins : 
for i=1:41  
    for j=1:length(AccumulatedmRNA)   %j for index of timepoint
        if i<13|i>17
            SmoothmRNA(j,i)=AccumulatedmRNA(j,i);
        else
        SmoothmRNA(j,i)=nanmean(AccumulatedmRNA(j,i-2:i+2));
        end
    end
end

%Plot the SmoothmRNA to check
%color code
jStart=1 %Protein.nc1
jEnd=length(AccumulatedmRNA);
jTotal=jEnd-jStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for k=1:length(AccumulatedmRNA)
    plot(0:0.025:1,SmoothmRNA(k,:),'color',Color(k-jStart+1,:));
    xlim([0.2 0.6]);
    %ylim([0 14000])
    %drawnow
    %pause(0.1)
    title('mRNA across AP')
    xlabel('AP axis')
    ylabel('mRNA (AU) ')
    set(gca,'fontsize',20)
    %saveas(gcf,['/Users/yangjoonkim/Documents/MATLAB/Boundary Problem/2016-11-23-Hb-nbGFP-MS2-NLSmCherry/AccumulatedmRNA/snapshots/AccumulatedmRNA',num2str(k)],'tiff')
end
hold off



%% Analyze the boundary position / width for AccumulatedmRNA
%In this dataset, I need to consider only 0.25~0.6 since 0-0.225 shows 0.
%Input : SmoothmRNA (time x APbin)
%For specific time point, SmoothmRNA(tpoint,:) -> mRNA at different time
%points.
APbin=[0:0.025:1];
tpoint=2;

% plot(0:0.025:1,SmoothmRNA(tpoint,:))
% xlim([0.25 0.6])
hold on
%Spline fitting
dxx=0.001;
xx=0.2:dxx:0.6;
yy=spline(APbin,SmoothmRNA(tpoint,:),xx);
yyy=pchip(APbin,SmoothmRNA(tpoint,:),xx);
yyyy=interp1(APbin,SmoothmRNA(tpoint,:),xx);
%figure(2)
hold on
plot(APbin,SmoothmRNA(tpoint,:),'o',xx,yy,'r')
plot(xx,yyy,'b')    %APbin,SmoothmRNA(tpoint,:),'o',
plot(xx,yyyy,'g')  %APbin,SmoothmRNA(tpoint,:),'o'
xlim([0.2 0.6])

legend('data','spline','pchip','interp')
%Find the x position where it has half of the Maximum intensity
for t=1:length(SmoothmRNA)
    MaxIntensity(t)=max(SmoothmRNA(t,:));
end

for i=1:length(xx)
    if yy(i)<0.5*MaxIntensity(tpoint)
        if yy(i-1)>0.5*MaxIntensity(tpoint)
            Xhalf=xx(i);
            HalfMaxY=yy(i);
            Index(tpoint)=i;
        end    
    end
end
Xhalf
%need to think about calculating the width.
line([Xhalf Xhalf],[0 HalfMaxY])
%ylim([0 19000])
title('Fitting Curve to the Data')
xlabel('AP')
ylabel('mRNA(AU)')
set(gca,'fontsize',30)
hold off
%% Plot for SmoothenedmRNA and boundary position for corresponding frames

SmoothTime=mRNA.ElapsedTime(1:length(SmoothmRNA));

iStart=2; %Start time
iEnd=length(mRNA.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

for tpoint=2:length(SmoothmRNA)
    
    APbin=[0:0.025:1];

    % plot(0:0.025:1,SmoothmRNA(tpoint,:))
    % xlim([0.25 0.6])
    subplot(1,4,1)
    hold on
    plot(APbin,SmoothmRNA(tpoint,:),'color',Color(tpoint-iStart+1,:))
    title('mRNA Accumulation')
    xlabel('AP')
    ylabel('mRNA Accumulation(AU)')
    xlim([0.2 0.6])
    ylim([0 30000])
    drawnow
    hold off
    
    %Spline fitting
    dxx=0.001;
    xx=0.2:dxx:0.6;
    yy=pchip(APbin,SmoothmRNA(tpoint,:),xx);
    
    %Re-define the pchip fit for width calculation
    xxx=0.2:dxx:0.6;
    yyy=pchip(APbin,SmoothmRNA(tpoint,:),xxx);
    
    
    %Find the x position where it has half of the Maximum intensity
    for t=1:length(SmoothmRNA)
        MaxIntensity(t)=max(SmoothmRNA(t,:));
    end

    for i=1:length(xx)
        if (yy(i)<0.5*MaxIntensity(tpoint)) & ((yy(i)==0));
            
        elseif (yy(i)<0.5*MaxIntensity(tpoint))
            if yy(i-1)>0.5*MaxIntensity(tpoint)|yy(i-1)==0.5*MaxIntensity(tpoint)
                Xhalf(tpoint)=xx(i);
                HalfMaxY(tpoint)=yy(i);
                Index(tpoint)=i;
            end    
        end
    end
    
    %Get the Width using xxx,yyy pchip fit
    %First, draw a tangent line at the midpoint.
    %Mid point
    MidX(tpoint)=Xhalf(tpoint);
    MidY(tpoint)=HalfMaxY(tpoint);
    
    %tangent line
    dX=0.001;
    X=0:dX:1;
    grad=diff(yy)./diff(xx);
    %grad(Index) : slope(gradient) at the midpoint.
    MidIndex=Index(tpoint);
    Y=grad(MidIndex)*(X-MidX(tpoint))+MidY(tpoint);
    
    Maximum=max(yyy);
    Minimum=min(yyy);
    
    for i=1:length(X)
        if Y(i)>Maximum
            if Y(i+1)<Maximum
                Left=i
            %else
                %print('sth went wrong')
            end
        end
    end
    
    for j=1:length(X)
        if Y(j)<Minimum
            if Y(j-1)>Minimum
                Right=j
            %else
                %print('sth went wrong')
            end
        end
    end
            
    Width(tpoint) = (Right-Left)*dX;
     
%     tpoint
%     Xhalf(tpoint)
    
    subplot(1,4,2)
    plot(APbin,SmoothmRNA(tpoint,:),'o','color',Color(tpoint-iStart+1,:))
    hold on
    plot(xx,yy,'color',Color(tpoint-iStart+1,:))
    xlim([0.2 0.6])
    ylim([0 30000])
    title('Boundary of mRNA')
    xlabel('AP')
    ylabel('mRNA Accumulation(AU)')
    %hold off
    
    %need to think about calculating the width.
    line([Xhalf(tpoint) Xhalf(tpoint)],[0 HalfMaxY(tpoint)])
    plot(X,Y,'k')
    bar(Maximum)
    bar(Minimum)
    line([Left*dX Right*dX],[Maximum+1 Maximum+1])
    
    hold off
    %drawnow
    %pause
    
    
    subplot(1,4,3)
    hold on
    plot(SmoothTime(tpoint),Xhalf(tpoint),'o','color',Color(tpoint-iStart+1,:))
    xlim([0 max(SmoothTime)+1])
    ylim([0.2 0.5])
    title({'Boundary Position','along Time'})
    xlabel('Time (min)')
    ylabel('mRNA Boundary')
    hold off
    %pause(0.5)
    
    subplot(1,4,4)
    hold on
    plot(SmoothTime(tpoint),Width(tpoint),'o','color',Color(tpoint-iStart+1,:))
    xlim([0 max(SmoothTime)+1])
    ylim([0 0.3])
    
    title({'Boundary Width','along Time'})
    xlabel('Time (min)')
    ylabel('mRNA Boundary Width')
    hold off
    
    %saveas(gcf,['/Users/yangjoonkim/Documents/MATLAB/Boundary Problem/2017-03-02-Hb-P2P-MS2V5-MCP-GFP/mRNA/AccumulatedmRNA',num2str(tpoint)],'tiff')

end

%v = VideoWriter('mRNA.avi');

%% function for calculating x position
hold on
for i=2:length(SmoothmRNA)
    Xhalf(i)=FindPos(SmoothmRNA,i);
end
hold off


%% Plot for X half along the time.
SmoothTime=mRNA.ElapsedTime(1:length(SmoothmRNA));

plotyy(SmoothTime,Xhalf,mRNA.ElapsedTime,mRNA.MeanVectorAll)
% bar(SmoothTime(nc11),1)
% bar(SmoothTime(nc12),1)
% bar(SmoothTime(nc13),1)
%bar(SmoothTime(nc14),1)
[hAx,hLine1,hLine2] = plotyy(SmoothTime,Xhalf,mRNA.ElapsedTime,mRNA.MeanVectorAll);

title('Boundary Position along the Time')
xlabel('Time (min)')

ylabel(hAx(1),'Boundary Position (AP)') % left y-axis
ylabel(hAx(2),'Spot Fluorescence','fontsize',30) % right y-axis


%ylim([0.39 0.47])
% title('Boundary Position along the Time')
% xlabel('Time (min)')
% ylabel('Boundary Position (AP)')
set(gca,'fontsize',30)
%Put a bar on the graph to indicate nuclear cycles.

%figure(2)
%plot(mRNA.ElapsedTime,mRNA.MeanVectorAll,'r')
%% Calibration for # of mRNA molecules (based on Garcia 2013, using hb P2P)
%Calibrate mRNA fluorescence to # of molecules
SmoothmRNA=SmoothmRNA/60; %I can change this more accurately later.

%% Calculate the amount of Protein by integrating the accumulatedmRNA (SmoothmRNA)

clear protein

%First, I need to define variables and parameters.
%gamma_P, r_p, Dp(kp)

%gamma_P : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
Phalf=50; %min, half-life of Protein, I can change this later.
gammaP=log(2)/Phalf; %Use Min. as time unit. I can change this later

%r_p:Translation rate for each AP bin, at nc14. 
rp=2; %(2 proteins/mRNA, min), I can change it later. from Petkova et al., 2014 Curr. Bio (Gregor lab)

%Dp:diffusion constant of Protein,
Dp=0*60; %[5:0.1:10]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=[1:41];
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
APbinpos = (AP-1)*dx;
kp=Dp/(dx)^2; %1/Min. unit.
Time = mRNA.ElapsedTime;
dt=mean(diff(Time));
%Make a zero-matrix
protein=zeros(length(SmoothmRNA),length(AP));
protein(1,:)=0;  %protein(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%P(i,1)=P(i-1,1)+rp(1)*dt-gammaP*P(i-1,1)*dt+k*P(i-1,2)*dt-k*P(i-1,1)*dt;
%P(i,41)=...
%But, right now, we don't have to worry too much, because we don't care
%about boundary now.


%Right now, I only have data from AP bin 0.225(10th)~0.525(24th), so only use this
%datasets

%Make a Protein reaction-diffusion equation
%I will assume that the translation occurs during the whole nc14.
for i=2:length(SmoothmRNA)  %for all Timepoints
    protein(i,9)=protein(i-1,9)+rp*SmoothmRNA(i-1,9)*dt-gammaP*protein(i-1,9)*dt; %+kp*protein(i-1,11)*dt-kp*protein(i-1,10)*dt;
    protein(i,20)=protein(i-1,20)+rp*SmoothmRNA(i-1,20)*dt-gammaP*protein(i-1,20)*dt; %+kp*protein(i-1,25)*dt-kp*protein(i-1,41)*dt;
    
    for j=10:19%length(AP)-1 %for all APbins
        protein(i,j)=protein(i-1,j)+rp*SmoothmRNA(i-1,j)*dt-gammaP*protein(i-1,j)*dt+kp*dt*(protein(i-1,j-1)+protein(i-1,j+1))-2*kp*protein(i-1,j)*dt;
    end
    
end

% plot for protein
jStart=1 %Protein.nc1
jEnd=length(SmoothmRNA);
jTotal=jEnd-jStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=1:length(protein)
    plot(0:0.025:1,protein(i,:),'color',Color(i-jStart+1,:))
    drawnow
    %pause(0.1)
    xlim([0.2 0.6])
    %ylim([0 4000])
end
hold off
title('Simulation of Protein concentration')
xlabel('Time(min)')
ylabel('Amount of Protein (AU)')
set(gca,'fontsize',30)

%% Plot for Protein (Prediction)boundary position for corresponding frames
clear PXhalf

SmoothTime=mRNA.ElapsedTime(1:length(SmoothmRNA));

iStart=2; %Start time
iEnd=length(SmoothmRNA); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
    
for tpoint=3:length(SmoothmRNA)
    
    APbin=[0:0.025:1];

    % plot(0:0.025:1,SmoothmRNA(tpoint,:))
    % xlim([0.25 0.6])
    subplot(1,4,1)
    hold on
    plot(APbin,protein(tpoint,:),'color',Color(tpoint,:))
    xlim([0.2 0.6])
    %ylim([0 10000])
    title('Protein over AP')
    xlabel('AP')
    ylabel('Protein Conc.(AU)')
    %drawnow
    
    %Spline fitting
    dxx=0.001;
    xx=0.3:dxx:0.6;
    yy=pchip(APbin,protein(tpoint,:),xx);
    
    xxx=0.2:dxx:0.6;
    yyy=pchip(APbin,protein(tpoint,:),xx);
    
    %Find the x position where it has half of the Maximum intensity
    for t=1:length(protein)
        MaxIntensity(t)=max(protein(t,:));
    end

    for i=1:length(xx)
        if (yy(i)<0.5*MaxIntensity(tpoint))&(yy(i)==0)
        elseif (yy(i)<0.5*MaxIntensity(tpoint))
            if yy(i-1)>0.5*MaxIntensity(tpoint)
                PXhalf(tpoint)=xx(i);
                HalfMaxY(tpoint)=yy(i);
            end    
        end
    end
    PXhalf(tpoint)
    %Get the Width using xxx,yyy pchip fit
    %First, draw a tangent line at the midpoint.
    %Mid point
    MidX(tpoint)=PXhalf(tpoint);
    MidY(tpoint)=HalfMaxY(tpoint);
    
    %tangent line
    dX=0.001;
    X=0:dX:1;
    grad=diff(yy)./diff(xx);
    %grad(Index) : slope(gradient) at the midpoint.
    MidIndex=Index(tpoint);
    Y=grad(MidIndex)*(X-MidX(tpoint))+MidY(tpoint);
    
    Maximum=max(yyy);
    Minimum=min(yyy);
    
    for i=1:length(X)
        if Y(i)>Maximum
            if Y(i+1)<Maximum
                Left=i
            %else
                %print('sth went wrong')
            end
        end
    end
    
    for j=1:length(X)
        if Y(j)<Minimum
            if Y(j-1)>Minimum
                Right=j
            %else
                %print('sth went wrong')
            end
        end
    end
            
    PWidth(tpoint) = (Right-Left)*dX;
    
    subplot(1,4,2)
    
    plot(APbin,protein(tpoint,:),'o','color',Color(tpoint,:))
    hold on
    plot(xx,yy,'color',Color(tpoint,:))
    xlim([0.2 0.6])
    title('Protein over AP')
    xlabel('AP')
    ylabel('Protein Conc.(AU)')
    %ylim([0 10000])
    %drawnow
    %need to think about calculating the width.
    line([PXhalf(tpoint) PXhalf(tpoint)],[0 HalfMaxY(tpoint)])
    %pause
    plot(X,Y,'k')
    bar(Maximum)
    bar(Minimum)
    line([Left*dX Right*dX],[Maximum+1 Maximum+1])
    
    hold off
    
    subplot(1,4,3)
    hold on
    plot(SmoothTime(tpoint),PXhalf(tpoint),'o','color',Color(tpoint-iStart+1,:))
    xlim([0 max(SmoothTime)+1])
    %ylim([0.42 0.46])
    title({'Boundary Position',' along Time'})
    xlabel('Time (min)')
    ylabel('Boundary Position(AP)')
    hold off
    drawnow
    %pause(0.5)
    
    subplot(1,4,4)
    hold on
    plot(SmoothTime(tpoint),PWidth(tpoint),'o','color',Color(tpoint-iStart+1,:))
    xlim([0 max(SmoothTime)+1])
    %ylim([0 0.3])
    
    title({'Boundary Width','along Time'})
    xlabel('Time (min)')
    ylabel('Boundary Width')
    hold off
    
    
    %saveas(gcf,['/Users/yangjoonkim/Documents/MATLAB/Boundary Problem/2017-03-14-Hb-P2P-MS2V5-NB-MCP-NLS-mCherry/Protein/Dp=0,6min/Protein',num2str(tpoint)],'tiff')
    
end



%% Protein simulation for different parameters
Dp=0:0.1:1;
Phalf=1:1:100;
dt=dt;
for i=1:length(Dp)
    for k=1:length(Phalf)
        PXshift(i,k)=ProteinSimulation(SmoothmRNA,Dp(i),Phalf(k),dt);
        
    end
end

        
%% function for calculating x position
hold on
for i=3:length(protein)
    PXhalf(i)=FindPos(protein,i);
end
hold off

%% Get the real protein fluorescence (concentration) data
Protein = load('CompiledNuclei.mat')
%% Protein(NB signal) 
%Extract the MeanVectorAP data as Protein at diff.time points, at different AP bins
NBprotein=Protein.MeanVectorAP;
NBprotein(isnan(NBprotein))=0;
NBTime=Protein.ElapsedTime;
nc14=Protein.nc14;
%hold on
for i=1:length(NBTime)
    plot (0:0.025:1,NBprotein(i,:))
    xlim([0.2 0.65])
    ylim([0 300])
    drawnow
    pause
end
hold off

%% Nanobody signal analysis - Time : how it increases
%first, I will plot Nanobody signal to time,
hold on
for i=1:41
    plot(Protein.ElapsedTime,Protein.MeanVectorAP(:,i))
    drawnow
    pause
end
hold off

%% NB signal at specific AP bin                                                                                                          
plot(Protein.ElapsedTime,NBprotein(:,10))
title('Nanobody signal along the Time(AP=0.25)')
xlabel('Time (min)')
ylabel('NB signal(AU)')
set(gca,'fontsize',30)


%% Nanobody signal analysis - AP (for different time points)
hold on
for i=1:length(Protein.ElapsedTime)
    plot(0:0.025:1,NBprotein(i,:))
    
    %text(1,max(NBprotein(i,:)),txt1)
    drawnow
    pause
end
    
 %% NB background subtraction
% % NB signal - background signal for corresponding frames. (This needs
% % longer imaging for background...since I shouldn't miss frames)
% BG=load('BG-CompiledNuclei.mat');
% BGFluo=BG.MeanVectorAP;
% BGTime=BG.ElapsedTime;
% BGnc13=BG.nc13;
% 
% % hold on
% % plot(BGTime,BGFluo(:,15))
% % plot(NBTime,NBprotein(:,15))
% % hold off
% 
% for i=3:length(BGTime)
%     NBprotein(i,:)=NBprotein(i,:)-BGFluo(i-2,:)
% end

%Subtracting background (preliminary version)
NBprotein=NBprotein-39 ;

NBprotein(NBprotein<0)=0;

%% NB signal Smoothening
% First, I need to smoothen the NB signal, and need to define the spatial
% averaging window.
Window=3;
% Second, I need to know which APbins I have. In this dataset( 2016-11-23,
% I have APbins from 10th to 16th.)
for i=1:length(NBprotein)
    for j=1:41
        if j<12|j>24
            SmoothNBprotein(i,j)=NBprotein(i,j);
        else 
            SmoothNBprotein(i,j)= nanmean(NBprotein(i,j-2:j+2));
        end
    end
end
%% Plot for Protein (NB data) boundary position for corresponding frames
%Right now, I only considered nc14, since I didn't figure out how to deal
%with Background for nc11,nc12,nc13

NBTime=Protein.ElapsedTime;
NBnc14=Protein.nc14;


iStart=2; %Start time
iEnd=length(NBTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
    
for tpoint=NBnc14+3:length(NBTime)
    
    APbin=[0:0.025:1];
    
    %Spline fitting
    dxx=0.001;
    xx=0.225:dxx:0.6;
    yy=pchip(APbin,SmoothNBprotein(tpoint,:),xx);
    
%     xxx=0.2:dxx:0.6;
%     yyy=pchip(APbin,SmoothNBprotein(tpoint,:),xxx);
    
    %Smoothing after interpolating (Change to this later)
%     for i=1:length(length(xxx))
%         if i<101|length(xxx)-100
%             Smoothyyy(tpoint,i)=yyy(i);
%         else
%             Smoothyyy(tpoint,i)=mean(yyy(i-100:i+100))
%         end
%     end

    % plot(0:0.025:1,SmoothmRNA(tpoint,:))
    % xlim([0.25 0.6])
    subplot(1,4,1)
    hold on
    plot(APbin,SmoothNBprotein(tpoint,:),'color',Color(tpoint,:))
    xlim([0.225 0.6])
    ylim([0 250])
    title('Protein over AP')
    xlabel('AP')
    ylabel('Protein Conc.(AU)')
    %drawnow
    
    
    
    %Find the x position where it has half of the Maximum intensity
    for t=1:length(NBTime)
        MaxIntensity(t)=max(SmoothNBprotein(t,:));
    end

    for i=2:length(xx)
        if (yy(i)<0.5*MaxIntensity(tpoint))&(yy(i)==0)
        elseif yy(i)<0.5*MaxIntensity(tpoint)
            if yy(i-1)>0.5*MaxIntensity(tpoint)
                NBXhalf(tpoint)=xx(i);
                HalfMaxY(tpoint)=yy(i);
                Index(tpoint)=i;
            end    
        end
    end
    NBXhalf(tpoint)
    
     %Get the Width using xxx,yyy pchip fit
    %First, draw a tangent line at the midpoint.
    %Mid point
    MidX(tpoint)=NBXhalf(tpoint);
    MidY(tpoint)=HalfMaxY(tpoint);
    
    %tangent line
    dX=0.001;
    X=0:dX:1;
    grad=diff(yy)./diff(xx);
    %grad(Index) : slope(gradient) at the midpoint.
    MidIndex=Index(tpoint);
    slope=mean(grad(MidIndex-10:MidIndex));
    Y=grad(MidIndex+1)*(X-MidX(tpoint))+MidY(tpoint);
    
    Maximum=max(yy);
    Minimum=min(yy);
    
    for i=1:length(X)
        if Y(i)>Maximum
            if Y(i+1)<Maximum
                Left=i
            %else
                %print('sth went wrong')
            end
        end
    end
    
    for j=1:length(X)
        if Y(j)<Minimum
            if Y(j-1)>Minimum
                Right=j
            %else
                %print('sth went wrong')
            end
        end
    end
            
    NBWidth(tpoint) = (Right-Left)*dX
    
    subplot(1,4,2)
    
    %plot(APbinID,NBprotein(tpoint,:),'o','color',Color(tpoint,:))
    plot(xx,yy,'color',Color(tpoint,:))
    
    hold on
    
    xlim([0.225 0.6])
    title('Protein over AP')
    xlabel('AP')
    ylabel('Protein Conc.(AU)')
    ylim([0 250])
    
    %drawnow
    %need to think about calculating the width.
    line([NBXhalf(tpoint) NBXhalf(tpoint)],[0 HalfMaxY(tpoint)])
    
    plot(X,Y,'k')
    bar(Maximum)
    bar(Minimum)
    line([Left*dX Right*dX],[Maximum+1 Maximum+1])

    hold off
    
    subplot(1,4,3)
    hold on
    plot(NBTime(tpoint)-NBTime(nc14),NBXhalf(tpoint),'o','color',Color(tpoint-iStart+1,:))
    xlim([0 max(NBTime)-NBTime(nc14)+1])
    ylim([0.38 0.5])
    title({'Protein Boundary','along Time'})
    xlabel('Time (min)')
    ylabel('NB protein Boundary')
    hold off
    drawnow
    %pause
    %pause(0.5)
    
    subplot(1,4,4) %Plot for Width
    hold on
    plot(NBTime(tpoint)-NBTime(nc14),NBWidth(tpoint),'o','color',Color(tpoint-iStart+1,:))
    xlim([0 max(NBTime)-NBTime(nc14)+1])
    ylim([0 0.4])
    
    title({'Boundary Width','along Time'})
    xlabel('Time (min)')
    ylabel('NB Protein Boundary Width')
    hold off
    
    %pause
    %saveas(gcf,['/Users/yangjoonkim/Documents/MATLAB/Boundary Problem/2016-11-23-Hb-nbGFP-MS2-NLSmCherry/NBProtein/NBProtein',num2str(tpoint)],'tiff')
    
end

%% Compare Smoothening vs Raw data
%At each time point, compare NBprotein vs SmoothNBprotein

for tpoint=NBnc14+5:length(NBTime)
    
    APbin=[0:0.025:1];

    % plot(0:0.025:1,SmoothmRNA(tpoint,:))
    % xlim([0.25 0.6])
    subplot(1,2,1)
    %hold on
    plot(APbin,NBprotein(tpoint,:),'color',Color(tpoint,:))
    xlim([0.25 0.6])
    ylim([0 250])
    title('Protein over AP')
    xlabel('AP')
    ylabel('Protein Conc.(AU)')
    drawnow
    
    subplot(1,2,2)
    plot(APbin,SmoothNBprotein(tpoint,:),'color',Color(tpoint,:))
    xlim([0.25 0.6])
    ylim([0 250])
    title('Smoothened Protein over AP')
    xlabel('AP')
    ylabel('Protein Conc.(AU)')
    %drawnow
end

%% Interpolation to the NB signal


%Right now, I only considered nc14, since I didn't figure out how to deal
%with Background for nc11,nc12,nc13

NBTime=Protein.ElapsedTime;
NBnc14=Protein.nc14;


iStart=2; %Start time
iEnd=length(NBTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
    
for tpoint=NBnc14+5:length(NBTime)
    
    APbin=[0:0.025:1];

    % plot(0:0.025:1,SmoothmRNA(tpoint,:))
    % xlim([0.25 0.6])
    subplot(1,4,1)
    hold on
    plot(APbin,NBprotein(tpoint,:),'color',Color(tpoint,:))
    xlim([0.25 0.6])
    ylim([0 250])
    title('Protein over AP')
    xlabel('AP')
    ylabel('Protein Conc.(AU)')
    %drawnow
    
    %Spline fitting
    dxx=0.001;
    xx=0.3:dxx:0.6;
    yy=interp1(APbin,NBprotein(tpoint,:),xx);
    
    xxx=0.25:dxx:0.6;
    yyy=interp1(APbin,NBprotein(tpoint,:),xxx);
    
    %Smoothening yyy (interpolated graph)
    for i=1:length(xxx)
        if i<5|i>(length(xxx)-4)
            Smoothyyy(i)=yyy(i);
        else
        Smoothyyy(i)=mean(yyy(i-4:i+4));
        end
    end
        
    
    %Find the x position where it has half of the Maximum intensity
    for t=1:length(NBprotein)
        MaxIntensity(t)=max(NBprotein(t,:));
    end

    for i=1:length(xx)
        if yy(i)<0.5*MaxIntensity(tpoint)
            if yy(i-1)>0.5*MaxIntensity(tpoint)
                NBXhalf(tpoint)=xx(i);
                HalfMaxY(tpoint)=yy(i);
            end    
        end
    end
    NBXhalf(tpoint)
    
     %Get the Width using xxx,yyy pchip fit
    %First, draw a tangent line at the midpoint.
    %Mid point
    MidX(tpoint)=NBXhalf(tpoint);
    MidY(tpoint)=HalfMaxY(tpoint);
    
    %tangent line (change yyy with Smoothyyy)
    dX=0.001;
    X=0:dX:1;
    grad=diff(yy)./diff(xx);
    %grad(Index) : slope(gradient) at the midpoint.
    MidIndex=Index(tpoint);
    Y=grad(MidIndex)*(X-MidX(tpoint))+MidY(tpoint);
    
    Maximum=max(Smoothyyy);
    Minimum=min(Smoothyyy);
    
    for i=1:length(X)
        if Y(i)>Maximum
            if Y(i+1)<Maximum
                Left=i
            %else
                %print('sth went wrong')
            end
        end
    end
    
    for j=1:length(X)
        if Y(j)<Minimum
            if Y(j-1)>Minimum
                Right=j
            %else
                %print('sth went wrong')
            end
        end
    end
            
    NBWidth(tpoint) = (Right-Left)*dX
    
    subplot(1,4,2)
    
    plot(APbin,NBprotein(tpoint,:),'o','color',Color(tpoint,:))
    hold on
    plot(xxx,Smoothyyy,'color',Color(tpoint,:))
    xlim([0.25 0.6])
    title('Protein over AP')
    xlabel('AP')
    ylabel('Protein Conc.(AU)')
    ylim([0 250])
    %drawnow
    %need to think about calculating the width.
    line([NBXhalf(tpoint) NBXhalf(tpoint)],[0 HalfMaxY(tpoint)])
    
    plot(X,Y,'k')
    bar(Maximum)
    bar(Minimum)
    line([Left*dX Right*dX],[Maximum+1 Maximum+1])
    
    
    hold off
    
    subplot(1,4,3)
    hold on
    plot(NBTime(tpoint),NBXhalf(tpoint),'o','color',Color(tpoint-iStart+1,:))
    xlim([0 max(NBTime)+1])
    ylim([0.38 0.46])
    title({'Protein Boundary','along Time'})
    xlabel('Time (min)')
    ylabel('Protein Boundary')
    hold off
    drawnow
    %pause
    %pause(0.5)
    
    subplot(1,4,4)
    hold on
    plot(SmoothTime(tpoint),NBWidth(tpoint),'o','color',Color(tpoint-iStart+1,:))
    xlim([0 max(SmoothTime)+1])
    %ylim([0.58 0.7])
    
    title({'Boundary Width','along Time'})
    xlabel('Time (min)')
    ylabel('mRNA Boundary Width')
    hold off
    
    
    %saveas(gcf,['/Users/yangjoonkim/Documents/MATLAB/Boundary Problem/2016-11-23-Hb-nbGFP-MS2-NLSmCherry/NBProtein-Width/NBProtein',num2str(tpoint)],'tiff')
    
end

%% Compare the Boundary shift for Simulation vs NB data
hold on
%1. NB data 

plot(NBTime(63:122),NBXhalf(63:122),'o','LineWidth',3)

%2. 
plot(SmoothTime(2:end),P1,'LineWidth',3)
plot(SmoothTime(2:end),P2,'LineWidth',3)
plot(SmoothTime(2:end),P3,'LineWidth',3)
plot(SmoothTime(2:end),P4,'LineWidth',3)
%plot(SmoothTime(2:end),P5,'LineWidth',3)
%plot(SmoothTime(2:end),P7,'LineWidth',3)
plot(SmoothTime(2:end),P9,'LineWidth',3)

ylim([0.35 0.5])
title('Boundary Position for different parameters')
xlabel('Time(min)')
ylabel('Boundary Position')
legend('NB data','0um^2/sec,60min','0.1um^2/sec,60min','1um^2/sec,60min','0um^2/sec,600min','0um^2/sec,3min','Location','northwest')
set(gca,'FontSize',30)
hold off


%% Plot for Protein(nanobody) boundary position for corresponding frames
%SmoothTime=mRNA.ElapsedTime(1:length(SmoothmRNA));

%NBprotein=NBprotein-40 ;

% for i=length(NBprotein)
%     for j=1:41
%         if NBprotein(i,j)<0
%             NBprotein(i,j)=0;
%         end
%     end
% end

iStart=2; %Start time
iEnd=length(NBprotein); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
    
%background subtraction

for tpoint=90:length(NBprotein)
    
    APbin=[0:0.025:1];

    % plot(0:0.025:1,SmoothmRNA(tpoint,:))
    % xlim([0.25 0.6])
    figure(1)
    hold on
    plot(APbin,protein(tpoint,:),'color',Color(tpoint,:))
    xlim([0.25 0.6])
    ylim([0 250])
    drawnow
   
    
    %Spline fitting
    dxx=0.001;
    xx=0.3:dxx:0.6;
    yy=spline(APbin,NBprotein(tpoint,:),xx);
    
    %Find the x position where it has half of the Maximum intensity
    for t=1:length(NBprotein)
        MaxIntensity(t)=max(NBprotein(t,:));
    end

    for i=1:length(xx)
        if yy(i)<0.5*MaxIntensity(tpoint)
            if yy(i-1)>0.5*MaxIntensity(tpoint)
                Xhalf=xx(i);
                HalfMaxY=yy(i);
            end    
        end
    end
    Xhalf
    
    figure(2)
    hold on
    plot(APbin,NBprotein(tpoint,:),'o','color',Color(tpoint,:))
    plot(xx,yy,'color',Color(tpoint,:))
    xlim([0.25 0.6])
    ylim([0 250])
    drawnow
    %need to think about calculating the width.
    line([Xhalf Xhalf],[0 HalfMaxY])
    pause
    hold off
    
%     hold on
%     figure(3)
%     for i=2:length(SmoothTime)
%         plot(SmoothTime(i),Xhalf,'o')
%     end
end




% %% background signal subtraction
% NBprotein=NBprotein-50 ;
% 
% for i=length(NBprotein)
%     for j=1:41
%         if NBprotein(i,j)<0
%             NBprotein(i,j)=0;
%         end
%     end
% end
%         
%% Protein (NB signal) smoothening (Averaging over 3 AP bins)
for i=2:length(ncTime)
    for j=1:41
        if j<3|j>39
            SmoothProtein(i,j)=NBprotein(i,j);
        else
        SmoothProtein(i,j)=nanmean(NBprotein(i,j-2:j+2));
        end
    end
end

hold on
for i=1:length(ncTime)
    plot (0:0.025:1,SmoothProtein(i,:))
    drawnow
end
hold off

%% Protein(NB signal)Boundary Calculation along the time 

for i=88:length(ncTime)
    NBpos(i)=FindPos(NBprotein,i);
end

plot(ncTime(nc14:end),NBpos(nc14:end))

%% Protein(NB signal)Boundary Calculation along the time 

for i=90 :length(ncTime)
    SmoothNBpos(i)=FindPos(SmoothProtein,i);
end

plot(ncTime,NBpos)

%% Plot for Boundary position for mRNA and Protein prediction
hold on
plot(Time(2:length(SmoothmRNA)),Xhalf(2:end),'b')
plot(Time(3:length(protein)),PXhalf(3:end),'r')
plot(Time,SmoothNBpos,'g')
ylim([0.4 0.5])

title('Boundary position along time')
xlabel('Time(min)')
ylabel('Boundary position for mRNA/Protein')
legend('mRNA','Protein','NB')
hold off
set(gca,'fontsize',30)

%% Interpolation of SmoothmRNA (for time)
%set the time gap
deltaT=0.01; %(min)
NewTime=0:deltaT:max(mRNA.ElapsedTime);
SmoothTime=round(SmoothTime,2);

%Make a matrix to save the interpolation data
InterpSmoothmRNA=zeros(length(NewTime),length(APbin));

%Interpolation for Time, at specific AP bin(2.5%)
for i=1:length(APbin)
    Interpolation=pchip(SmoothTime,SmoothmRNA(:,i),NewTime);
    InterpSmoothmRNA(:,i)=Interpolation;
end

%% 