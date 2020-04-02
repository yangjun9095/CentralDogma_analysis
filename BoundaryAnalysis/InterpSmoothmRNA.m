%Interpolation of AccumulatedmRNA
%% mRNA Accumulation from Spot fluorescence
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~4 frames in ElapsedTime)
lambdam=log(2)./60; %mRNA half-life = 60min
MRNA=mRNA.MeanVectorAP;
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
xlim([0.2 0.6])
set(gca,'fontsize',30)

%% Interpolation 
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

