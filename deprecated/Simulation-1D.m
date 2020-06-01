%1D Simulation Generator (mRNA, Protein)


%% Transcriptional Rate
%Define step-funciton like transcription activity
%Start from thinking embryo as 1D boxes, and divide 10um, such that the
%width is around diameter of one nucleus. (500um as embryo length)

for i=1:50
    if i<26
        Rate(i)=20;  %RNA Pol II loading rate, [molecules/min], Garcia et al. 2013, Curr.Bio. 
    else
        Rate(i)=0;
    end
end
APbin=[0.01:0.02:0.99];

plot(APbin,Rate,'k','LineWidth',5)
%plot(APbin,zeros(1,41),'k','LineWidth',2)
ylim([0 25])
title('Transcriptional Rate')
xlabel('AP position')
ylabel('Transcription Rate (mRNA molecules/min)')
set(gca, 'FontSize', 30)

%% mRNA simulation-1D
%Using the CalculatemRNA.m function to calculate the mRNA for various Dm and Tm
%Rate is given as a step-function
clear mRNA
% Dm=[0,0.1,1,5,10];
% Tm=[1,10,50,100,1000];

Dm=[0:0.1:0.9,1:1:10];
Tm=[1:9,10:10:50,100:100:1000];

mRNA=cell(length(Dm),length(Tm));

for i=1:length(Dm)
    for j=1:length(Tm)
        M=CalculatemRNA(Rate,Dm(i),Tm(j),60);
        mRNA(i,j)={M};
    end
end

%% Plot to check different mRNA profiles
jStart=1;
jEnd=length(Tm)+1;
jTotal=jEnd-jStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=1:length(Tm)
    for j=1:length(Dm)
        M=cell2mat(mRNA(j,i));
        
        plot(0.01:0.02:0.99,M(1:500:6001,:),'color',Color(j-jStart+1,:))
        pause
    end
end

%% Boundary position and width as a function of Dm and Tm
clear Bposition
clear BWidth

% Dm=[0,0.1,1,5,10];
% Tm=[1,10,50,100,1000];

Dm=[0:0.1:0.9,1:1:10];
Tm=[1:9,10:10:50,100:100:1000];

%mRNA=cell(length(Dm),length(Tm));

for i=1:length(Dm)
    for j=1:length(Tm)
        M=cell2mat(mRNA(i,j));
        [Xhalf,Width] = GetBoundary(M);
        Bposition(i,j)={Xhalf};
        BWidth(i,j)={Width};
    end
end

%% 3D plot for the boundary position

clear Boundary

% Dm=[0,0.1,1,5,10];
% Tm=[1,10,50,100,1000];

Dm=[0:0.1:0.9,1:1:10];
Tm=[1:9,10:10:50,100:100:1000];

for i=1:length(Dm)
    for j=1:length(Tm)
        B=cell2mat(Bposition(i,j));
        for k=1:6 %number of time points
            Boundary(i,j,k)=B(k*1000+1);
        end
    end
end

%Colorful surfacees for time
% jStart=1;
% jEnd=6+1;
% jTotal=jEnd-jStart;
% 
% colormap(jet(256));
% cmap=colormap;
% Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);
% %Color(j-jStart+1,:)

hold on
%for l=1:6
l=6;
    C=Boundary(:,:,l); %For coloring according to the Z-value
    hSurface=surfc(Tm,Dm,Boundary(:,:,l),C)
    caxis([0.5 0.8])
    axis([0 1000 0 10 0.5 0.8])  %Z range : 0.5~0.8
    set(gca, 'XScale', 'log', 'YScale', 'log')
    %set(hSurface, 'FaceColor',Color(l,:), 'FaceAlpha',1, 'EdgeAlpha', 0);
    
    title('Boundary Position from Simulation')
    xlabel('Tm (min)')
    ylabel('Dm (\mu m^2/sec)')
    zlabel('Boundary Position (AP)')
    %legend('10 min','20 min','30 min','40 min','50 min','60 min')
    
    set(gca, 'Fontsize', 30)
    ax=gca
    ax.XGrid = 'on'
    ax.YGrid = 'on'
    ax.ZGrid = 'on'
    
    %pause
%end

%% 3D plot for the boundary width
clear BoundaryWidth

% Dm=[0,0.1,1,5,10];
% Tm=[1,10,50,100,1000];

Dm=[0:0.1:0.9,1:1:10];
Tm=[1:9,10:10:50,100:100:1000];

for i=1:length(Dm)
    for j=1:length(Tm)
        B=cell2mat(BWidth(i,j));
        for k=1:6
            BoundaryWidth(i,j,k)=B(k*1000+1);
        end
    end
end

hold on
% for l=1:6
l=6;
    C=BoundaryWidth(:,:,l); %For coloring according to the Z-value
    hSurface=surf(Tm,Dm,BoundaryWidth(:,:,l),C)
    caxis([0 1])
    axis([0 1000 0 10 0 1])
    set(gca, 'XScale', 'log', 'YScale', 'log')

    %set(hSurface, 'FaceColor',Color(l,:), 'FaceAlpha',1, 'EdgeAlpha', 0);
    
    title('mRNA Boundary Width')
    xlabel('Tm (min)')
    ylabel('Dm (\mu m^2/sec)')
    zlabel('Boundary Width (AP)')
    %legend('10 min','20 min','30 min','40 min','50 min','60 min')
    
    set(gca, 'Fontsize', 30)
    ax=gca
    ax.XGrid = 'on'
    ax.YGrid = 'on'
    ax.ZGrid = 'on'
    %pause
% end

%% mRNA input (Dm=0,Tm=60min) - sanity check for boundary
clear Xhalf
clear Width

%M=CalculatemRNA(Rate,0,60,60);
M=cell2mat(mRNA(11,11));
[Xhalf,Width,yyy,Left,Right,Yhalf,yy]=GetBoundary(M);


xxx=0:0.001:1;
%Plot mRNA profile at different time points, with boundary position and
%width
figure(2)
hold on
for t=1:6
    plot(0.01:0.02:0.99,M(1000*t+1,:),'o')
    plot(xxx,yyy(1000*t+1,:))
    plot(xxx,yy(1000*t+1,:))
    line([Xhalf(1000*t+1) Xhalf(1000*t+1)],[0 Yhalf(1000*t+1)])
    line([Left(1000*t+1)*0.001 Right(1000*t+1)*0.001],[max(M(1000*t+1,:))+1 max(M(1000*t+1,:))+1])
    ylim([-1 1000])
    pause
end

%% 
hold on
%Plot Boundary position along time
plot((1001:1000:6001)*0.01,Xhalf(1001:1000:6001))
%Plot for Width
plot((1001:1000:6001)*0.01,Width(1001:1000:6001)) 

%% Protein simulation (1D)

%First, I need certain mRNA profile to start with.
%Let's start with Dm=0,Tm=60min

M=CalculatemRNA(Rate,5,60,60);

clear Protein
% Dm=[0,0.1,1,5,10];
% Tm=[1,10,50,100,1000];

Dp=[0:0.1:0.9,1:1:10];
Tp=[1:9,10:10:50,100:100:1000];


Protein=cell(length(Dp),length(Tp));

for i=1:length(Dp)
    for j=1:length(Tp)
        P=CalculateProtein(M,Dp(i),Tp(j),60); %Assuming that the translation is keep occuring 
        Protein(i,j)={P};
    end
end

%% Plot to check different Protein profiles
jStart=1;
jEnd=length(Tp)+1;
jTotal=jEnd-jStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=1:length(Tp)
    for j=1:length(Dp)
        P=cell2mat(Protein(j,i));
        
        plot(0.01:0.02:0.99,P(1:500:6001,:),'color',Color(j-jStart+1,:))
        pause
    end
end

%% Boundary position and width as a function of Dp and Tp
clear PBposition
clear PBWidth

Dp=[0:0.1:0.9,1:1:10];
Tp=[1:9,10:10:50,100:100:1000];

%mRNA=cell(length(Dm),length(Tm));

for i=1:length(Dp)
    for j=1:length(Tp)
        P=cell2mat(Protein(i,j));
        [Xhalf,Width] = GetBoundary(P);
        PBposition(i,j)={Xhalf};
        PBWidth(i,j)={Width};
    end
end

%% 3D plot for the boundary position (Protein)
Dp=[0:0.1:0.9,1:1:10];
Tp=[1:9,10:10:50,100:100:1000];

for i=1:length(Dp)
    for j=1:length(Tp)
        B=cell2mat(PBposition(i,j));
        for k=1:6 %number of time points
            Boundary(i,j,k)=B(k*1000+1);
        end
    end
end

%Colorful surfacees for time
jStart=1;
jEnd=6+1;
jTotal=jEnd-jStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);
%Color(j-jStart+1,:)

hold on
%for l=1:6
l=6;
    C=Boundary(:,:,l);
    hSurface=surf(Tp,Dp,Boundary(:,:,l),C)
    caxis([0.5 0.8])
    axis([0 1000 0 10 0.5 0.8])
    set(gca, 'XScale', 'log', 'YScale', 'log')
    %set(hSurface, 'FaceColor',[0,0,1], 'FaceAlpha',0.8, 'EdgeAlpha', 0);
    
    title('Protein Boundary Position')
    xlabel('Tp (min)')
    ylabel('Dp (\mu m^2/sec)')
    zlabel('Boundary Position (AP)')
    %legend('10 min','20 min','30 min','40 min','50 min','60 min')
    
    set(gca, 'Fontsize', 30)
    ax=gca
    ax.XGrid = 'on'
    ax.YGrid = 'on'
    ax.ZGrid = 'on'
    %pause
%end

%% 3D plot for the boundary width (Protein)

Dp=[0:0.1:0.9,1:1:10];
Tp=[1:9,10:10:50,100:100:1000];


for i=1:length(Dp)
    for j=1:length(Tp)
        B=cell2mat(PBWidth(i,j));
        for k=1:6
            BoundaryWidth(i,j,k)=B(k*1000+1);
        end
    end
end

%hold on
%for l=1:6
l=6;

    C=BoundaryWidth(:,:,l); %For coloring according to the Z-value
    
    hSurface=surf(Tp,Dp,BoundaryWidth(:,:,l),C)
    caxis([0 1])
    axis([0 1000 0 10 0 1])
    set(gca, 'XScale', 'log', 'YScale', 'log')
    
    %set(hSurface, 'FaceColor',[0,0,1], 'FaceAlpha',0.8, 'EdgeAlpha', 0);
    
    title('Boundary Width from Simulation')
    xlabel('Tp (min)')
    ylabel('Dp (\mu m^2/sec)')
    zlabel('Protein Boundary Width (AP)')
    legend('10 min','20 min','30 min','40 min','50 min','60 min')
    
    set(gca, 'Fontsize', 30)
    ax=gca
    ax.XGrid = 'on'
    ax.YGrid = 'on'
    ax.ZGrid = 'on'
    %pause
%end


%% Comparison of mRNA and protein (Boundary Position)

%This code will compare the boundary position of mRNA and protein at time
%60min (at the end of nc 14)

%Bposition vs PBposition

%% Plot for mRNA and Protein at the same time
hold on
for i=1:500:length(M)
    %hold on
    subplot(1,2,1)
    plot(0.01:0.02:0.99,M(i,:))
    
    subplot(1,2,2)
    plot(0.01:0.02:0.99,P(i,:))
    title('Title')
    xlabel('x axis')
    ylabel('y axis')
    legend('legend')
    hold off
    pause
    
    
end
%% Protein increase over time

for j=1:length(APbin)
    hold on
    subplot(1,2,1)
    plot((1:6001)*0.01,M(:,j))
    hold off
    
    hold on
    subplot(1,2,2)
    plot((1:6001)*0.01,P(:,j))
    hold off
    
    pause
end

