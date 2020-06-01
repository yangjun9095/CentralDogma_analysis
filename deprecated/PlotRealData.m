%Work on the real data that I took.
%First, I need to load the data as mRNA and Protein (I need to fix this
%later)
mRNA = load('CompiledParticles.mat')
Protein = load('CompiledNuclei.mat')

%% First, plot for Protein vs AP axis over the time

jStart=Protein.nc14;
jEnd=length(Protein.ElapsedTime)+1;
jTotal=jEnd-jStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=Protein.nc14:length(Protein.ElapsedTime)
    plot(Protein.APbinID,Protein.MeanVectorAP(i,:)-38,'color',Color(i-jStart+1,:))
end
hold off
ylim([0 230]) 
xlim([0.2 0.6])
title('Protein Concentration over AP position')
xlabel('AP axis')
ylabel('Protein Concentration(AU)')
set(gca,'fontsize',60)
%% Smoothening the protein gradient(for AP position) for finding boundary position and slope(width)
for i=1:41
    for j=1:length(Protein.ElapsedTime)
        if i<3|i>39
            SmoothProtein(j,i)=Protein.MeanVectorAP(j,i);
        else
        SmoothProtein(j,i)=mean(Protein.MeanVectorAP(j,i-2:i+2));
        end
    end
end

%% Analysis (Protein)

for i=90:length(Protein.ElapsedTime)     %Protein.nc14:
    Boundarypos(i)=FindProteinPosition(SmoothProtein(i,:),i);
    
end
hold on
plot(Protein.ElapsedTime(90:end)-Protein.ElapsedTime(Protein.nc14),Boundarypos(90:end),'o')
%plot()
hold off
xlim([0 60]);
ylim([0 1]);
title('Boundary Position for Protein')
xlabel('Time(min)')
ylabel('Boundary Position (AP axis)')

%% Second, we have to plot polII loading rate
%Use mRNA (CompiledParticles.MeanVectorAP)

%for colorful plots
jStart=1;
jEnd=length(mRNA.ElapsedTime)+1;
jTotal=jEnd-jStart;

colormap(jet(256));
cmap=colormap;
Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=1:length(mRNA.ElapsedTime)
    plot(mRNA.ElapsedTime(i),mRNA.MeanVectorAll(1,i),'o-','color',Color(i-jStart+1,:))
end
hold off

title('Active PolII over AP position','Fontsize',40)
xlabel('AP axis','Fontsize',30)
ylabel('Active PolII(AU)','Fontsize',30)