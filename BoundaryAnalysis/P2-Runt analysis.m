%Analyzing Runt - hbP2 + 0,1,2,3 Runt binding sites 

%% Load the data set
clear all

mRNAr0=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r0/CompiledParticles.mat')

mRNAr1=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r1/CompiledParticles.mat')

mRNAr2=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r2/CompiledParticles.mat');

mRNAr3=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r3/CompiledParticles.mat');

%% Plot mRNA along the timecourse (for specific AP bin)

ap=10; % choose specific AP bin

%Use MeanVectorAP at APbinID = ap, and plot them for timecourse 
hold on
plot(mRNAr1.ElapsedTime,mRNAr1.MeanVectorAP(:,ap),'r')
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
plot(mRNAr0.ElapsedTime,mRNAr0.MeanVectorAll,'b')
%plot(mRNAr2.ElapsedTime,mRNAr2.MeanOffsetVector,'r')
hold off

title('Spot Fluorescence and Offset')
legend('MS2 spot fluorscence')%,'Offset')
xlabel('Time(min)')
ylabel('Fluorescence (AU)')
set(gca,'fontsize',30)

%% mRNA Accumulation from Spot fluorescence (r0)
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./180; %mRNA half-life = 60min
MRNA=mRNAr0.MeanVectorAP;
MRNA(isnan(MRNA))=0;
Nt=length(mRNAr0.ElapsedTime); %# of time points.
AccumulatedmRNA=zeros(Nt-4,41);
%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
%didn't consider diffusion in here.
for i=2:Nt-4
    for j=1:41
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+(mRNAr0.ElapsedTime(i)-mRNAr0.ElapsedTime(i-1))*MRNA(i,j) - lambdam*AccumulatedmRNA(i-1,j)*(mRNAr0.ElapsedTime(i)-mRNAr0.ElapsedTime(i-1));
    end
end

% plot
%Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(mRNAr0.ElapsedTime); %End time
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

%% mRNA Accumulation from Spot fluorescence (r1)
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./180; %mRNA half-life = 60min
MRNA=mRNAr1.MeanVectorAP;
MRNA(isnan(MRNA))=0;
Nt=length(mRNAr1.ElapsedTime); %# of time points.
AccumulatedmRNA=zeros(Nt-4,41);
%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
%didn't consider diffusion in here.
for i=2:Nt-4
    for j=1:41
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+(mRNAr1.ElapsedTime(i)-mRNAr1.ElapsedTime(i-1))*MRNA(i,j) - lambdam*AccumulatedmRNA(i-1,j)*(mRNAr1.ElapsedTime(i)-mRNAr1.ElapsedTime(i-1));
    end
end

% plot
%Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(mRNAr1.ElapsedTime); %End time
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

%% mRNA Accumulation from Spot fluorescence (r2)
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./180; %mRNA half-life = 60min
MRNA=mRNAr2.MeanVectorAP;
MRNA(isnan(MRNA))=0;
Nt=length(mRNAr2.ElapsedTime); %# of time points.
AccumulatedmRNA=zeros(Nt-4,41);
%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
%didn't consider diffusion in here.
for i=2:Nt-4
    for j=1:41
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+(mRNAr2.ElapsedTime(i)-mRNAr2.ElapsedTime(i-1))*MRNA(i,j) - lambdam*AccumulatedmRNA(i-1,j)*(mRNAr2.ElapsedTime(i)-mRNAr2.ElapsedTime(i-1));
    end
end

% plot
%Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(mRNAr1.ElapsedTime); %End time
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

%% mRNA Accumulation from Spot fluorescence (r3)
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./180; %mRNA half-life = 60min
MRNA=mRNAr3.MeanVectorAP;
MRNA(isnan(MRNA))=0;
Nt=length(mRNAr3.ElapsedTime); %# of time points.
AccumulatedmRNA=zeros(Nt-4,41);
%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
%didn't consider diffusion in here.
for i=2:Nt-4
    for j=1:41
        AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+(mRNAr3.ElapsedTime(i)-mRNAr3.ElapsedTime(i-1))*MRNA(i,j) - lambdam*AccumulatedmRNA(i-1,j)*(mRNAr3.ElapsedTime(i)-mRNAr3.ElapsedTime(i-1));
    end
end

% plot
%Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(mRNAr1.ElapsedTime); %End time
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
        if i<11|i>22
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
%% Plot for SmoothenedmRNA and boundary position for corresponding frames (X_half-maximum)

SmoothTime=mRNAr0.ElapsedTime(1:length(SmoothmRNA));

iStart=2; %Start time
iEnd=length(mRNAr0.ElapsedTime); %End time
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
    xx=0.225:dxx:0.6;
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

%% MeanVectorAP (RNAP loading rate over time)
for ap=10:27
    plot(mRNAr2.ElapsedTime,mRNAr2.MeanVectorAP(:,ap),'r')
    hold on
    plot(mRNAr3.ElapsedTime-10,mRNAr3.MeanVectorAP(:,ap),'b')
    hold off
    title('RNAP loading rate')
    xlabel('AP')
    ylabel('RNAP loading rate(AU)')
    pause
    
end

%% MeanVectorAP (over AP)
iStart=1; %Start time
iEnd=length(mRNAr2.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
    
figure(1)
hold on
for t=1:length(mRNAr3.ElapsedTime)
    plot(0:0.025:1,mRNAr2.MeanVectorAP(t,:),'color',Color(t-iStart+1,:));
end

%% MeanVectorAP (over AP)
iStart=1; %Start time
iEnd=length(mRNAr1.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
%figure(2)
hold on
for t=1:length(mRNAr1.ElapsedTime)
    plot(0:0.025:1,mRNAr1.MeanVectorAP(t,:),'color',Color(t-iStart+1,:));
    drawnow
end

%% Plot all
for k=1:40
        
    hold on
    plot(0:0.025:1,mRNAr0.MeanVectorAP(mRNAr0.nc14+k,:),'k')
    plot(0:0.025:1,mRNAr1.MeanVectorAP(mRNAr1.nc14+k,:),'b')
    plot(0:0.025:1,mRNAr2.MeanVectorAP(mRNAr2.nc14+k,:),'g')
    plot(0:0.025:1,mRNAr3.MeanVectorAP(mRNAr3.nc14+k,:),'r')
    
    title('mRNA rate over AP')
    xlabel('AP')
    ylabel('mRNA (AU)')
    legend('0','1','2','3')
    set(gca,'Fontsize',30)
    xlim([0.2 0.6])
    pause
end

%% Plot MeanVectorAP



hold off


%% Smoothening more for inflection points

%Basically, it's averaging AccumulatedmRNA (AP) over 5 adjacent AP bins
clear SmoothmRNA

%range of AP bins : 
for i=1:41  
    for j=1:length(AccumulatedmRNA)   %j for index of timepoint
        if i<12|i>22
            SmoothmRNA(j,i)=AccumulatedmRNA(j,i);
        else
        SmoothmRNA(j,i)=nanmean(AccumulatedmRNA(j,i-3:i+3));
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


%% Plot for SmoothenedmRNA and boundary position for corresponding frames (inflection point)
SmoothTime=mRNAr0.ElapsedTime(1:length(SmoothmRNA));

iStart=2; %Start time
iEnd=length(mRNAr0.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

for tpoint=30:length(SmoothmRNA)
    
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
    xx=0.3:dxx:0.6;
    yy=pchip(APbin,SmoothmRNA(tpoint,:),xx);
    
    
    %Find the x position where it has half of the Maximum intensity
    
    for i=1:length(xx)-1
        gradient(i)=(yy(i+1)-yy(i))/(xx(i+1)-xx(i));
    end
    
    for k=1:length(xx)-2
        doublegradient(k)=gradient(k+1)-gradient(k);
    end
    
    j=1;
    value=1;
    while value==1;
        if doublegradient(j)>0;
            Xhalf(tpoint)=xx(j);
            HalfMaxY(tpoint)=yy(j);
            value=2;
        else
            j=j+1;
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

%% MeanFit
%% Load the data set
clear all

mRNAr0=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r0/MeanFits.mat')

mRNAr1=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r1/MeanFits.mat')

mRNAr2=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r2/MeanFits.mat');

mRNAr3=load('/Users/yangjoonkim/Documents/MATLAB/hbP2-r0123/hbP2-r3/MeanFits.mat');

%% FitResults
r0=mRNAr0.FitResults;
r1=mRNAr1.FitResults;
r2=mRNAr2.FitResults;
r3=mRNAr3.FitResults;
%% Plotting Mean Rate of transcription at nc 13 along the AP aixs

for i=1:length(r0)
    if r0(i,2).Approved==1
        r00(i)=r0(i,2).RateFit;
    elseif i<22
        r00(i)=nan;
    else
        r00(i)=0;
    end
end

for i=1:length(r1)
    if r1(i,2).Approved==1
        r11(i)=r1(i,2).RateFit;
    elseif i<10
        r11(i)=nan;
    else
        r11(i)=0;
    end
end

for i=1:length(r2)
    if r2(i,2).Approved==1
        r22(i)=r2(i,2).RateFit;
    elseif i<10
        r22(i)=nan;
    else
        r22(i)=0;
    end
end

for i=1:length(r3)
    if r3(i,2).Approved==1
        r33(i)=r3(i,2).RateFit;
    elseif i<10
        r33(i)=nan;
    else
        r33(i)=0;
    end
end

%% Averaging over 3 AP bins
for j=2:length(r0)-1
    r000(j)=nanmean(r00(j-1:j+1));
end

for j=2:length(r1)-1
    r111(j)=nanmean(r11(j-1:j+1));
end

for j=2:length(r2)-1
    r222(j)=nanmean(r22(j-1:j+1));
end

for j=2:length(r3)-1
    r333(j)=nanmean(r33(j-1:j+1));
end

r000(41)=0;
r111(41)=0;
r222(41)=0;
r333(41)=0;
%% 
hold on
plot(0:0.025:1,r000,'k')
plot(0:0.025:1,r111,'b')
plot(0:0.025:1,r222,'g')
plot(0:0.025:1,r333,'r')

xlim([0 1])
ylim([0 250])
title('Rate of Transcription for different number of repressor binding sites')
xlabel('AP')
ylabel('Rate of Transcription (AU)')
legend('0','1','2','3')

set(gca,'Fontsize',30)
hold off