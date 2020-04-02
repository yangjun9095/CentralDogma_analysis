%Plot for Figure.3.
%Boundary position and width for different time points
iStart=2; %Start time
iEnd=length(mRNA.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

Timepoint=[30,44,58,73,87];


for tpoint=1:length(Timepoint)
    
    MaxIntensity(Timepoint(tpoint))=max(SmoothmRNA(Timepoint(tpoint),:));
        %Spline fitting
    dxx=0.001;
    xx=0.3:dxx:0.6;
    yy=pchip(APbin,SmoothmRNA(Timepoint(tpoint),:),xx);
    
    figure(1)
    hold on
    plot(APbin,SmoothmRNA(Timepoint(tpoint),:),'o','color',Color(Timepoint(tpoint)-iStart+1-10,:))
   
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
    plot(ncTime(Timepoint(tpoint)),NBXhalf(Timepoint(tpoint)),'o','color',Color(Timepoint(tpoint)-iStart+1-10,:))
    xlim([0 max(SmoothTime)+1])
    ylim([0.37 0.42])
    title({'Boundary Position','along Time'})
    xlabel('Time (min)')
    ylabel('mRNA Boundary')
    hold off
    %pause(0.5)
    set(gca,'fontsize',30)
    
    figure(3)
    hold on
    plot(SmoothTime(Timepoint(tpoint)),Width(Timepoint(tpoint)),'o','color',Color(Timepoint(tpoint)-iStart+1-10,:))
    xlim([0 max(SmoothTime)+1])
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
    xx=0.2:dxx:0.6;
    yy=pchip(APbin,protein(Timepoint(tpoint),:),xx);
    
    figure(1)
    hold on
    plot(APbin,protein(Timepoint(tpoint),:),'o','color',Color(Timepoint(tpoint)-iStart+1,:))
   
    plot(xx,yy,'color',Color(Timepoint(tpoint)-iStart+1,:))
    xlim([0.2 0.6])
    ylim([0 4000])
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
    ylim([0.37 0.42])
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

%% Nanobody signal plot
Protein=load('/Users/yangjoonkim/Documents/MATLAB/Boundary Problem/2017-03-14-Hb-P2P-MS2V5-NB-MCP-NLS-mCherry/CompiledNuclei')

%% Boundary Width (stripe+P2 vs P2P)
hold on
plot(StripeNBTime(Stripenc14+1:end)-StripeNBTime(Stripenc14),StripeNBWidth(Stripenc14:end),'r','LineWidth',5)
plot(P2PNBTime(P2Pnc14+1:end)-P2PNBTime(P2Pnc14),P2PNBWidth(P2Pnc14:end),'b','LineWidth',5)
title('Boundary Width over Time')
xlabel('Time after nc14')
ylabel('Boundary Width')
legend('P2+stripe','P2P')
set(gca,'Fontsize',30)


%% 
hold on
plot(StripeNBTime(Stripenc14+5:end)-StripeNBTime(Stripenc14),StripeNBXhalf(Stripenc14+4:end),'r','LineWidth',5)
plot(P2PNBTime(P2Pnc14+5:end)-P2PNBTime(P2Pnc14),P2PNBXhalf(P2Pnc14+4:end),'b','LineWidth',5)
title('Boundary Position over Time')
xlabel('Time after nc14')
ylabel('Boundary Position')
legend('P2+stripe','P2P')
set(gca,'Fontsize',30)
