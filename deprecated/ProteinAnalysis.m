% New analysis code, example
%% Load the dataset
mRNA = load('CompiledParticles.mat')
Protein = load('CompiledNuclei.mat')

%% Protein Boundary Calculation along the time 
%Extract the MeanVectorAP data as Protein at diff.time points, at different AP bins
NBprotein=Protein.MeanVectorAP;
ncTime=Protein.ElapsedTime;
nc14=Protein.nc14;

for i=nc14:length(ncTime)
    Pos(i)=sum(NBprotein(i,:)>max(NBprotein(i,:))*0.5);
    Bpos(i)=Pos(i)*0.025;
end

plot(ncTime(nc14:end),Bpos(nc14:end))

%% First, plot for Protein vs AP axis over the time

jStart=1 %Protein.nc14;
jEnd=length(Protein.ElapsedTime)+1;
jTotal=jEnd-jStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=Protein.nc14:length(Protein.ElapsedTime)
    plot(Protein.APbinID,Protein.MeanVectorAP(i,:),'color',Color(i-jStart+1,:))  %need to consider background level
end
hold off
%ylim([0 230]) 
xlim([0.2 0.6])
title('Protein Concentration over AP position')
xlabel('AP axis')
ylabel('Protein Concentration(AU)')
set(gca,'fontsize',60)
%% Smoothening the protein gradient(for AP position) for finding boundary position and slope(width)
%Basically, it's averaging nuclear fluorescence (AP) over 5 adjacent AP bins

for i=1:41
    for j=1:length(Protein.ElapsedTime)
        if i<3|i>39
            SmoothProtein(j,i)=Protein.MeanVectorAP(j,i);
        else
        SmoothProtein(j,i)=nanmean(Protein.MeanVectorAP(j,i-2:i+2));
        end
    end
end

%% Plot the SmoothProtein for AP axis (along the time) (either MeanVectorAP(T(i),:) or SmoothProtein(T(i),:)
T=[30,44,58,72,87]; %Time point index for ncTime (20min,30min,40min,50min,60min)

hold on
for i=1:length(T)
    plot(Protein.APbinID,Protein.MeanVectorAP(T(i),:),'color',Color(T(i)-jStart+1-10,:))  %need to consider background level
end
hold off
%ylim([0 230]) 
xlim([0.2 0.6])
title('Protein Concentration over AP position')
xlabel('AP axis')
ylabel('Protein Concentration(AU)')
set(gca,'fontsize',60)


%% Plot 
%Boundary position and width for different time points
iStart=2; %Start time
iEnd=length(mRNA.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

Timepoint=[30,44,58,73,87]; %(20min,30min,40min,50min,60min)

for tpoint=1:length(Timepoint)
    
        %Spline fitting
    dxx=0.001;
    xx=0.3:dxx:0.6;
    yy=pchip(APbin,SmoothProtein(Timepoint(tpoint),:),xx);
    
    figure(1)
    hold on
    plot(APbin,SmoothProtein(Timepoint(tpoint),:),'o','color',Color(Timepoint(tpoint)-iStart+1,:))
   
    plot(xx,yy,'color',Color(Timepoint(tpoint)-iStart+1,:))
    xlim([0.2 0.6])
    %ylim([0 30000])
    title('Boundary of mRNA')
    xlabel('AP')
    ylabel('mRNA Accumulation(AU)')
    
    hold off
    %drawnow
    %pause
    set(gca,'fontsize',30)
    
    figure(2)
    hold on
    plot(ncTime(Timepoint(tpoint)),Xhalf(Timepoint(tpoint)),'o','color',Color(Timepoint(tpoint)-iStart+1,:))
    xlim([0 max(ncTime)+1])
    %ylim([0.37 0.42])
    title({'Boundary Position','along Time'})
    xlabel('Time (min)')
    ylabel('mRNA Boundary')
    hold off
    %pause(0.5)
    set(gca,'fontsize',30)
    
    figure(3)
    hold on
    plot(ncTime(Timepoint(tpoint)),Width(Timepoint(tpoint)),'o','color',Color(Timepoint(tpoint)-iStart+1,:))
    xlim([0 max(ncTime)+1])
    ylim([0 0.3])
    
    title({'Boundary Width','along Time'})
    xlabel('Time (min)')
    ylabel('mRNA Boundary Width')
    hold off
    set(gca,'fontsize',30)
end

%% Protein plot (for different time points)
iStart=2; %Start time
iEnd=length(mRNA.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

Timepoint=[30,44,58,73,87];

for tpoint=1:length(Timepoint)
    
        %Spline fitting
    dxx=0.001;
    xx=0.3:dxx:0.6;
    yy=pchip(APbin,protein(Timepoint(tpoint),:),xx);
    
    figure(1)
    hold on
    plot(APbin,protein(Timepoint(tpoint),:),'o','color',Color(Timepoint(tpoint)-iStart+1,:))
   
    plot(xx,yy,'color',Color(Timepoint(tpoint)-iStart+1,:))
    xlim([0.2 0.6])
    %ylim([0 4000])
    title('Boundary of Protein')
    xlabel('AP')
    ylabel('Protein(AU)')
    
    hold off
    %drawnow
    %pause
    set(gca,'fontsize',30)
    
    figure(2)
    hold on
    plot(SmoothTime(Timepoint(tpoint)),PXhalf(Timepoint(tpoint)),'o','color',Color(Timepoint(tpoint)-iStart+1,:))
    xlim([0 max(SmoothTime)+1])
    %ylim([0.37 0.42])
    title({'Boundary Position','along Time'})
    xlabel('Time (min)')
    ylabel('Protein Boundary')
    hold off
    %pause(0.5)
    set(gca,'fontsize',30)
    
    figure(3)
    hold on
    plot(SmoothTime(Timepoint(tpoint)),PWidth(Timepoint(tpoint)),'o','color',Color(Timepoint(tpoint)-iStart+1,:))
    xlim([0 max(SmoothTime)+1])
    ylim([0 0.3])
    
    title({'Boundary Width','along Time'})
    xlabel('Time (min)')
    ylabel('Protein Boundary Width')
    hold off
    set(gca,'fontsize',30)
end

