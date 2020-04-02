%CentralDogma : Prediction for mRNA and protein from given transcriptional
%rate profile.
%(written by Yang Joon Kim)

%Assumptions:
%1) I am predicting the number of mRNA transcripts or protein molecules in
%a nucleus. Thus, I can use the MeanVectorAP (averaged over only active
%nuclei). Later, to calculate the total number of mRNA / protein over AP, I
%need to consider the fraction of active nuclei.

%First, we assume the transcription rate to be a step-function (Perfect
%On/Off fashion)

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

plot(APbin,Rate,'k','LineWidth',2)
%plot(APbin,zeros(1,41),'k','LineWidth',2)
ylim([0 25])
title('Transcription Rate')
xlabel('AP position')
ylabel('Transcription Rate (mRNA molecules/min)')
set(gca, 'FontSize', 30)


%% Second, simulation for mRNA (use Dm, gammaM)
clear M
clear T
%Let's start with considering nc14 only.
%First, I need to define variables and parameters.

Tm=60; %half-life of mRNA, measured by injecting alpha-amanitin
%gamma_M : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
gammaM=log(2)/Tm; %Use Min. as time unit.(I can change this later, for nc14)
%r_m:Transcription rate for each AP bin, at nc14. from Hernan's 2013 Curr.Bio paper
rm=Rate; %(:,3);   %Normalized by # of RNAP

%Dm:diffusion constant of mRNA,
Dm=0*60; %0.1*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=1:50;
dx=500/50; %delta x = 500um/50(bins) = 10um/bin (diameter of one nucleus)
APbinpos = (AP-1)*dx+5; %Center of nuclei, defined as position of that Bin.
km=Dm/((dx)^2); %1/Min. unit.

%Time window
dt=0.01;
T=[0:dt:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.

%We need Transcription Time Window, along the AP axis, as in Fig.4B (1x41)
%Twindow=[0,0,0,0,0,0,0,0,18,17,17.5,17,19,17,15.5,16,15.5,12.5,11,11.5,10,8.5,11,8.5,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6];

% %Fraction of Active Nuclei Fig.4.F
% FractionActiveNuclei=MaxOnRatio(:,3);
% FractionActiveNuclei(isnan(FractionActiveNuclei))=0;

%Effective rate of transcription along the AP axis.
%rmeff=rm.*FractionActiveNuclei;

%Make a zero-matrix
M=zeros(length(T),length(AP));
M(1,:)=0;  %mRNA(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).

%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Make a mRNA reaction-diffusion equation
for i=2:length(T) %for all Timepoints
    for j=1:length(AP)  %for all AP bins
        
        if i<length(T)*1/3
            if j==1
                M(i,1)=M(i-1,1)+rm(j)*dt-gammaM*M(i-1,1)*dt+km*M(i-1,2)*dt-km*M(i-1,1)*dt;
            elseif j==50
                M(i,50)=M(i-1,50)+rm(j)*dt-gammaM*M(i-1,50)*dt+km*M(i-1,49)*dt-km*M(i-1,50)*dt;
            else
            M(i,j)=M(i-1,j)+rm(j)*dt-gammaM*M(i-1,j)*dt+km*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km*M(i-1,j)*dt;
            end
        else
            if j==1
                M(i,1)=M(i-1,1)-gammaM*M(i-1,1)*dt+km*M(i-1,2)*dt-km*M(i-1,1)*dt;
            elseif j==50
                M(i,50)=M(i-1,50)-gammaM*M(i-1,50)*dt+km*M(i-1,49)*dt-km*M(i-1,50)*dt;
            else M(i,j)=M(i-1,j)-gammaM*M(i-1,j)*dt+km*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km*M(i-1,j)*dt;
            end
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
for l=1:500:length(T)
    plot([0.01:0.02:0.99],M(l,:),'color',Color(l-jStart+1,:));
%     ylim([0,1000])
%     pause(0.01)
end
hold off

title('mRNA along the AP axis')
xlabel('AP axis')
ylabel('mRNA (molecules)')
set(gca, 'FontSize', 30)

%% For various Dm (mRNA Diffusion constants)
clear M
clear T
clear mRNA
%Let's start with considering nc14 only.
%First, I need to define variables and parameters.

Tm=60; %half-life of mRNA, measured by injecting alpha-amanitin
%gamma_M : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
gammaM=log(2)/Tm; %Use Min. as time unit.(I can change this later, for nc14)
%r_m:Transcription rate for each AP bin, at nc14. from Hernan's 2013 Curr.Bio paper
rm=Rate; %(:,3);   %Normalized by # of RNAP

%Dm:diffusion constant of mRNA,
Dm=[0,0.1,1,5,10]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=1:50;
dx=500/50; %delta x = 500um/50(bins) = 10um/bin (diameter of one nucleus)
APbinpos = (AP-1)*dx+5; %Center of nuclei, defined as position of that Bin.
km=Dm/((dx)^2); %1/Min. unit.

%Time window
dt=0.01;
T=[0:dt:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.

mRNA=cell(1,length(km));

%Make a mRNA reaction-diffusion equation
for k=1:length(km)
    clear M
    M=zeros(length(T),length(AP));
    M(1,:)=0;  %mRNA(t=0)=0 at all APbins. Initial Condition
    for i=2:length(T) %for all Timepoints
        for j=1:length(AP)  %for all AP bins
        
            if i<length(T)*1/3
                if j==1
                    M(i,1)=M(i-1,1)+rm(j)*dt-gammaM*M(i-1,1)*dt+km(k)*M(i-1,2)*dt-km(k)*M(i-1,1)*dt;
                elseif j==50
                    M(i,50)=M(i-1,50)+rm(j)*dt-gammaM*M(i-1,50)*dt+km(k)*M(i-1,49)*dt-km(k)*M(i-1,50)*dt;
                else
                M(i,j)=M(i-1,j)+rm(j)*dt-gammaM*M(i-1,j)*dt+km(k)*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km(k)*M(i-1,j)*dt;
                end
            else
                if j==1
                    M(i,1)=M(i-1,1)-gammaM*M(i-1,1)*dt+km(k)*M(i-1,2)*dt-km(k)*M(i-1,1)*dt;
                elseif j==50
                    M(i,50)=M(i-1,50)-gammaM*M(i-1,50)*dt+km(k)*M(i-1,49)*dt-km(k)*M(i-1,50)*dt;
                else M(i,j)=M(i-1,j)-gammaM*M(i-1,j)*dt+km(k)*dt*(M(i-1,j-1)+M(i-1,j+1))-2*km(k)*M(i-1,j)*dt;
                end
            end
        end
    end
    mRNA(k)={M};
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
for k=1:length(km)
    M=cell2mat(mRNA(k));
    for l=1:500:length(T)
        plot([0.01:0.02:0.99],M(l,:),'color',Color(l-jStart+1,:));
    %     ylim([0,1000])
    %     pause(0.01)
    
    end
    pause
end
hold off

title('mRNA along the AP axis')
xlabel('AP axis')
ylabel('mRNA (molecules)')
set(gca, 'FontSize', 30)

%% Using the CalculatemRNA.m function to calculate the mRNA for various Dm and Tm
%Rate is given as step-function
Dm=[0,0.1,1,5,10];
Tm=[1,10,50,100,1000];

mRNA=cell(length(Dm),length(Tm));

for i=1:length(Dm)
    for j=1:length(Tm)
        M=CalculatemRNA(Rate,Dm(i),Tm(j),20);
        mRNA(i,j)={M};
    end
end

%% Plot to check different mRNA profiles
% jStart=1;
% jEnd=length(Tm)+1;
% jTotal=jEnd-jStart;
% 
% colormap(jet(256));
% cmap=colormap;
% Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for i=1:length(Tm)
    for j=1:length(Dm)
        M=cell2mat(mRNA(j,i));
        
        plot(0.01:0.02:0.99,M(1:500:6001,:),'color',Color(j-jStart+1,:))
        pause
    end
end

%% Boundary position and width as a function of Dm and Tm
Dm=[0,0.1,1,5,10];
Tm=[1,10,50,100,1000];

%mRNA=cell(length(Dm),length(Tm));

for i=1:length(Dm)
    for j=1:length(Tm)
        M=cell2mat(mRNA(i,j));
        [Xhalf,Width] = GetBoundary(M);
        Bposition(i,j)={Xhalf};
        BWidth(i,j)={Width};
    end
end


%% 3D plot (X:Dm,Y:Tm,Z: Xhalf or Width)
Dm=[0,0.1,1,5,10];
Tm=[1,10,50,100,1000];

Z1=Bposition;
Z2=BWidth;



%% To capture the boundary position and width at one specific time point

%I need to interpolate the x axis (AP), and I will use pchip for
%interpolation. Alternatives : interp1, 

M=cell2mat(mRNA(5));

%hold on
for T=2:6001
    MRNA(T,:)=M(T,:);
    %plot([0.01:0.02:0.99],MRNA(T,:))
    
    clear Left
    clear Right
    
    clear yyy
    %At one time point, T, get the boundary position and width
    
    %First, I need to do interpolation
    xxx=0:0.001:1;
    
    yyy(T,:)=pchip(APbin,MRNA(T,:),xxx);
    
    for i=1:length(xxx)
        if (yyy(T,i)<0.5*max(MRNA(T,:)))|(yyy(T,i)==0.5*max(MRNA(T,:)))
            if yyy(T,i-1)>0.5*max(MRNA(T,:))
                Xhalf(T)=xxx(i);
                Yhalf(T)=yyy(T,i);
                Slope(T)=(yyy(T,i+1)-yyy(T,i-1))/(xxx(i+1)-xxx(i-1));
            end
        else
            Xhalf(T)=0;
            Yhalf(T)=0;
            Slope(T)=0;
        end
        
        %Tangential line at the boundary position
        yy(T,:)=Slope(T)*(xxx-Xhalf(T))+Yhalf(T);
    end
    
    for j=1:length(xxx)
        if (yy(T,:)<max(MRNA(T,:)))+(yy(T,:)==max(MRNA(T,:))) %In case the tangent line is always smaller than the mRNA profile
            Left=0;
        elseif yy(T,j)>max(MRNA(T,:))
            if yy(T,j+1)<max(MRNA(T,:))
                Left=j;
            end
        end
        
        if (yy(T,:)>min(MRNA(T,:)))+(yy(T,:)==min(MRNA(T,:)))
            Right=0;
        elseif yy(T,j)<min(MRNA(T,:))
            if yy(T,j-1)>min(MRNA(T,:))
                Right=j;
            end
        end
        
    end
    
    Width(T)=(Right-Left)*0.001;
    Sharp(T)=abs(Slope(T));
    
end
%hold off
%% Plot to check Xhalf and Width
hold on
plot(0:0.01:60,Xhalf)
plot(0:0.01:60,Width)
%plot(0:0.01:60,Width2)
%% Plot to check curve fitting
hold on
plot(xxx,yy(5757,:))
plot(xxx,yyy(5757,:))

%% Check if the boundary position and width calculation is correct.
%I will check if the boundary position and width that I am calculating is
%correct, by plotting them at different time points.

% parameters :
%1) mRNA profile: plot(APbin,M(T,:))
%2) fitted line : plot(APbin,yy)

for tpoint=2:6001
    plot(APbin,M(tpoint,:),'o','color',Color(tpoint,:))
    hold on
    plot(xxx,yyy(tpoint,:),'color',Color(tpoint,:))
    %xlim([0.2 0.6])
    title('mRNA over AP')
    xlabel('AP')
    ylabel('mRNA (molecules)')
    %ylim([0 10000])
    %drawnow
    %need to think about calculating the width.
    line([Xhalf(tpoint) Xhalf(tpoint)],[0 Yhalf(tpoint)])
    %pause
    %plot(X,Y,'k')
    %bar(Maximum)
    %bar(Minimum)
    %line([Left*dX Right*dX],[Maximum+1 Maximum+1])
    pause
    hold off
end
%% Boundary Position (half-maximum)
hold on
plot(0:0.01:60,Xhalf1,'r')
%pause
plot(0:0.01:60,Xhalf2,'k')
%pause
plot(0:0.01:60,Xhalf3,'b')
%pause
plot(0:0.01:60,Xhalf4,'g')
%pause
plot(0:0.01:60,Xhalf5)

ylim([0.4 0.8])
set(gca,'Fontsize',30)
title('Boundary Position along Time')
xlabel('Time (min)')
ylabel('Boundary Position (AP)')
legend('0 \mu m^2/sec','0.1\mu m^2/sec','1\mu m^2/sec','5\mu m^2/sec','10\mu m^2/sec')
%% Boundary Width
hold on
plot(0:0.01:60,Width1,'r')
pause
plot(0:0.01:60,Width2,'k')
pause
plot(0:0.01:60,Width3,'b')
pause
plot(0:0.01:60,Width4,'g')
pause
plot(0:0.01:60,Width5)

set(gca,'Fontsize',30)
title('Boundary Width along Time')
xlabel('Time (min)')
ylabel('Boundary Width(AP)')
legend('0 \mu m^2/sec','0.1\mu m^2/sec','1\mu m^2/sec','5\mu m^2/sec','10\mu m^2/sec')
%% Sharpness
hold on
plot(0.01:0.01:60,Sharp1(2:end),'r')
pause
plot(0.01:0.01:60,Sharp2(2:end),'k')
pause
plot(0.01:0.01:60,Sharp3(2:end),'b')
pause
plot(0.01:0.01:60,Sharp4(2:end),'g')
pause
plot(0.01:0.01:60,Sharp5(2:end))

set(gca,'Fontsize',30)
title('Boundary Sharpness along Time')
xlabel('Time (min)')
ylabel('Boundary Sharpness(AP)')

%% Third, Protein simulation 
clear P
%Let's start with considering nc14 only.

%I am calculating the protein profile from mRNA profile (predicted from
%previous code)
M=cell2mat(mRNA(1));
%First, I need to define variables and parameters.
%gamma_P, r_p, Dp(kp)
%gamma_P : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
Tp=5; %min
gammaP=log(2)/Tp; %Use Min. as time unit. I can change this later
%r_p:Translation rate for each AP bin, at nc14. 
rp=2; %(2 proteins/mRNA, min), I can change it later. from Petkova et al., 2014 Curr. Bio (Gregor lab)

%Dp:diffusion constant of Protein,
Dp=[0,0.1,1,5,10]*60; %[5:0.1:10]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=[1:50];
dx=500/50; %delta x = 500um/50(bins) = 10um/bin
APbinpos = (AP-1)*dx+5;
kp=Dp/(dx)^2; %1/Min. unit.


%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%P(i,1)=P(i-1,1)+rp(1)*dt-gammaP*P(i-1,1)*dt+k*P(i-1,2)*dt-k*P(i-1,1)*dt;
%P(i,41)=...
%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Time window
dt=0.01;
T=[0:dt:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.

Protein=cell(1,length(kp));

%Make a Protein reaction-diffusion equation
%I will assume that the translation occurs during the whole nc14.
for k=1:length(kp)
    %Make a zero-matrix
    P=zeros(length(T),length(AP));
    P(1,:)=0;  %Protein(t=0)=0 at all APbins. Initial Condition
    
    for i=2:length(T)  %for all Timepoints
        P(i,1)=P(i-1,1)+rp*M(i-1,1)*dt-gammaP*P(i-1,1)*dt+kp(k)*P(i-1,2)*dt-kp(k)*P(i-1,1)*dt;
        P(i,50)=P(i-1,50)+rp*M(i-1,50)*dt-gammaP*P(i-1,50)*dt+kp(k)*P(i-1,49)*dt-kp(k)*P(i-1,50)*dt;
    
        for j=2:length(AP)-1 %for all APbins
            P(i,j)=P(i-1,j)+rp*M(i-1,j)*dt-gammaP*P(i-1,j)*dt+kp(k)*dt*(P(i-1,j-1)+P(i-1,j+1))-2*kp(k)*P(i-1,j)*dt;
        end
    
    end
    
    Protein(k)={P};
end

%% Simulation - Protein profile for different half-life
clear P
clear Protein
%Let's start with considering nc14 only.

%I am calculating the protein profile from mRNA profile (predicted from
%previous code)
M=cell2mat(mRNA(1));
%First, I need to define variables and parameters.
%gamma_P, r_p, Dp(kp)
%gamma_P : mRNA degradation rate, log(2)/half-life of mRNA (~60min)
Tp=[0.1,1,10,100,1000]; %min, half-life of protein
gammaP=log(2)./Tp; %Use Min. as time unit. I can change this later
%r_p:Translation rate for each AP bin, at nc14. 
rp=2; %(2 proteins/mRNA, min), I can change it later. from Petkova et al., 2014 Curr. Bio (Gregor lab)

%Dp:diffusion constant of Protein,
Dp=5*60; %[5:0.1:10]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016

AP=[1:50];
dx=500/50; %delta x = 500um/50(bins) = 10um/bin
APbinpos = (AP-1)*dx+5;
kp=Dp/(dx)^2; %1/Min. unit.


%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%P(i,1)=P(i-1,1)+rp(1)*dt-gammaP*P(i-1,1)*dt+k*P(i-1,2)*dt-k*P(i-1,1)*dt;
%P(i,41)=...
%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Time window
dt=0.01;
T=[0:dt:60]; %Time matrix, spanning from 0min to 60min, with 0.01 min gap.

Protein=cell(1,length(gammaP));

%Make a Protein reaction-diffusion equation
%I will assume that the translation occurs during the whole nc14.
for l=1:length(gammaP)
    %Make a zero-matrix
    P=zeros(length(T),length(AP));
    P(1,:)=0;  %Protein(t=0)=0 at all APbins. Initial Condition
    
    for i=2:length(T)  %for all Timepoints
        P(i,1)=P(i-1,1)+rp*M(i-1,1)*dt-gammaP(l)*P(i-1,1)*dt+kp*P(i-1,2)*dt-kp*P(i-1,1)*dt;
        P(i,50)=P(i-1,50)+rp*M(i-1,50)*dt-gammaP(l)*P(i-1,50)*dt+kp*P(i-1,49)*dt-kp*P(i-1,50)*dt;
    
        for j=2:length(AP)-1 %for all APbins
            P(i,j)=P(i-1,j)+rp*M(i-1,j)*dt-gammaP(l)*P(i-1,j)*dt+kp*dt*(P(i-1,j-1)+P(i-1,j+1))-2*kp*P(i-1,j)*dt;
        end
    
    end
    
    Protein(l)={P};
end

%% Plot for check if simulation makes sense
%colorful plot for timecourse : To see how the shape changes along the
%time. 
 jStart=1;
 jEnd=length(T)+1;
 jTotal=jEnd-jStart;

 colormap(jet(256));
 cmap=colormap;
 Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);
    
hold on
for k=1:length(gammaP)
    P=cell2mat(Protein(k));
    for l=1:500:length(T)
        plot([0.01:0.02:0.99],P(l,:),'color',Color(l-jStart+1,:));
    %     ylim([0,1000])
    %     pause(0.01)
    
    end
    pause
end
hold off

title('Protein along the AP axis')
xlabel('AP axis')
ylabel('Protein (molecules)')
set(gca, 'FontSize', 30)


%% Capture the boundary position and width for Protein
%I need to interpolate the x axis (AP), and I will use pchip for
%interpolation. Alternatives : interp1, 
clear P
clear PXhalf
clear PWidth
P=cell2mat(Protein(5)); %select protein profile for specific Dp or Tp

%hold on
for T=3:6001
    %P(T,:)=P(T,:);
    %plot([0.01:0.02:0.99],MRNA(T,:))
    
    clear Left
    clear Right
    
    clear yyy
    %At one time point, T, get the boundary position and width
    
    %First, I need to do interpolation
    xxx=0:0.001:1;
    
    yyy(T,:)=pchip(APbin,P(T,:),xxx);
    
    for i=1:length(xxx)
        if (yyy(T,i)<0.5*max(P(T,:)))|(yyy(T,i)==0.5*max(P(T,:)))
            if yyy(T,i-1)>0.5*max(P(T,:))
                PXhalf(T)=xxx(i);
                PYhalf(T)=yyy(T,i);
                PSlope(T)=(yyy(T,i+1)-yyy(T,i-1))/(xxx(i+1)-xxx(i-1));
            end
        else
            PXhalf(T)=0;
            PYhalf(T)=0;
            PSlope(T)=0;
        end
        
        %Tangential line at the boundary position
        yy(T,:)=PSlope(T)*(xxx-PXhalf(T))+PYhalf(T);
    end
    
    for j=1:length(xxx)
        %if (yy(T,:)<max(P(T,:)))+(yy(T,:)==P(MRNA(T,:))) %In case the tangent line is always smaller than the mRNA profile
        %    Left=0;
        if yy(T,j)>max(P(T,:))
            if yy(T,j+1)<max(P(T,:))
                Left=j;
            end
        end
        
        %if (yy(T,:)>min(P(T,:)))+(yy(T,:)==min(P(T,:)))
        %    Right=0;
        if yy(T,j)<min(P(T,:))
            if yy(T,j-1)>min(P(T,:))
                Right=j;
            end
        end
        
    end
    
    PWidth(T)=(Right-Left)*0.001;
    PSharp(T)=abs(Slope(T));
    
end
%hold off

%% Check for Protein boundary position and width
%% Plot to check Xhalf and Width
hold on
plot(0:0.01:60,PXhalf)
plot(0:0.01:60,PWidth)
%plot(0:0.01:60,Width2)
%% Plot to check curve fitting
hold on
plot(xxx,yy(6001,:))
plot(xxx,yyy(6001,:))

%% Check if the boundary position and width calculation is correct.
%I will check if the boundary position and width that I am calculating is
%correct, by plotting them at different time points.

% parameters :
%1) mRNA profile: plot(APbin,M(T,:))
%2) fitted line : plot(APbin,yy)

for tpoint=2:6001
    plot(APbin,P(tpoint,:),'o','color',Color(tpoint,:))
    hold on
    plot(xxx,yyy(tpoint,:),'color',Color(tpoint,:))
    %xlim([0.2 0.6])
    title('mRNA over AP')
    xlabel('AP')
    ylabel('mRNA (molecules)')
    %ylim([0 10000])
    %drawnow
    %need to think about calculating the width.
    line([(PXhalf(tpoint)-0.5*PWidth(tpoint)) (PXhalf(tpoint)+0.5*PWidth(tpoint))],[PYhalf(tpoint) PYhalf(tpoint)])
    %pause
    %plot(X,Y,'k')
    %bar(Maximum)
    %bar(Minimum)
    %line([Left*dX Right*dX],[Maximum+1 Maximum+1])
    pause
    hold off
end
%% Protein Boundary Position (half-maximum)
hold on
plot(0:0.01:60,PXhalf1,'r')
%pause
plot(0:0.01:60,PXhalf2,'k')
%pause
plot(0:0.01:60,PXhalf3,'b')
%pause
plot(0:0.01:60,PXhalf4,'g')
%pause
plot(0:0.01:60,PXhalf5)

ylim([0.4 0.8])
set(gca,'Fontsize',30)
title('Boundary Position along Time')
xlabel('Time (min)')
ylabel('Boundary Position (AP)')
legend('0 \mu m^2/sec','0.1\mu m^2/sec','1\mu m^2/sec','5\mu m^2/sec','10\mu m^2/sec')

%% Protein Boundary Position (half-maximum) - half-life of 5min
hold on
plot(0:0.01:60,PXhalf11,'r')
%pause
plot(0:0.01:60,PXhalf12,'k')
%pause
plot(0:0.01:60,PXhalf13,'b')
%pause
plot(0:0.01:60,PXhalf14,'g')
%pause
plot(0:0.01:60,PXhalf15)

ylim([0.4 0.8])
set(gca,'Fontsize',30)
title('Boundary Position along Time')
xlabel('Time (min)')
ylabel('Boundary Position (AP)')
legend('0 \mu m^2/sec','0.1\mu m^2/sec','1\mu m^2/sec','5\mu m^2/sec','10\mu m^2/sec')


%% Protein Boundary Position (half-maximum) - different life time
hold on
plot(0:0.01:60,PXhalf6)
%pause
plot(0:0.01:60,PXhalf7)
%pause
plot(0:0.01:60,PXhalf8)
%pause
plot(0:0.01:60,PXhalf9)
%pause
plot(0:0.01:60,PXhalf10)

ylim([0.4 0.8])
set(gca,'Fontsize',30)
title('Boundary Position along Time')
xlabel('Time (min)')
ylabel('Boundary Position (AP)')
legend('0.1 min','1 min','10 min','100 min','1000 min')
%% Protein Boundary Width - Dp
hold on
plot(0:0.01:60,PWidth1,'r')
pause
plot(0:0.01:60,PWidth2,'k')
pause
plot(0:0.01:60,PWidth3,'b')
pause
plot(0:0.01:60,PWidth4,'g')
pause
plot(0:0.01:60,PWidth5)

ylim([0 0.9])
set(gca,'Fontsize',30)
title('Boundary Width along Time')
xlabel('Time (min)')
ylabel('Boundary Width(AP)')
legend('0 \mu m^2/sec','0.1\mu m^2/sec','1\mu m^2/sec','5\mu m^2/sec','10\mu m^2/sec')


%% Protein Boundary Width - Dp (for protein half-life of 5 min)
hold on
plot(0:0.01:60,PWidth11,'r')
pause
plot(0:0.01:60,PWidth12,'k')
pause
plot(0:0.01:60,PWidth13,'b')
pause
plot(0:0.01:60,PWidth14,'g')
pause
plot(0:0.01:60,PWidth15)

ylim([0 0.9])
set(gca,'Fontsize',30)
title('Boundary Width along Time')
xlabel('Time (min)')
ylabel('Boundary Width(AP)')
legend('0 \mu m^2/sec','0.1\mu m^2/sec','1\mu m^2/sec','5\mu m^2/sec','10\mu m^2/sec')

%% Protein Boundary Width - half-life

hold on
plot(0:0.01:60,PWidth6)
pause
plot(0:0.01:60,PWidth7)
pause
plot(0:0.01:60,PWidth8)
pause
plot(0:0.01:60,PWidth9)
pause
plot(0:0.01:60,PWidth10)

ylim([0 0.9])
set(gca,'Fontsize',30)
title('Boundary Width along Time')
xlabel('Time (min)')
ylabel('Boundary Width(AP)')
legend('0.1 min','1 min','10 min','100 min','1000 min')
%% Plot for protein
hold on
for k=1:length(kp)
    P=cell2mat(Protein(k));
    for l=1:500:length(T)
        plot([0.01:0.02:0.99],P(l,:),'color',Color(l-jStart+1,:));
    %     ylim([0,1000])
    %     pause(0.01)
        title('Protein along the AP axis')
    xlabel('AP axis')
    ylabel('Protein (molecules)')
    set(gca, 'FontSize', 30)

    end
    pause
end
hold off



%% colorful plot for Protein : To see how the shape changes along the time. 
Tindex=[1:100:6001];

%colorful plot for timecourse : To see how the shape changes along the
%time. 
 jStart=1;
 jEnd=length(T)+1;
 jTotal=jEnd-jStart;

 colormap(jet(256));
 cmap=colormap;
 Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for k=Tindex
    plot([0.01:0.02:0.99],P(k,:),'color',Color(k-jStart+1,:));
end
hold off
%xlim([0.2 0.6])
title('Protein over AP axis(simulation)')
xlabel('AP axis')
ylabel('Protein Concentration(molecules)')
set(gca,'fontsize',30)
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
    %ylim([0 70000])
end

%% Plotting Protien for specific time point
t=[0,15,30,45,60];
tpoint=t/dt+1;

for i=1:length(tpoint)
    figure(i)
    plot(APbin,P(tpoint(i),:),'b','LineWidth',2)
    
    title('Protein','Fontsize',40)
    xlabel('AP position','Fontsize',30)
    ylabel('Protein(AU)','Fontsize',30)
    %ylim([0 2500000])
end

%% colorful plot for Protein : To see how the shape changes along the time. 
Tindex=[1:50:6001];

%colorful plot for timecourse : To see how the shape changes along the
%time. 
 jStart=1;
 jEnd=length(T)+1;
 jTotal=jEnd-jStart;

 colormap(jet(256));
 cmap=colormap;
 Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);

hold on
for k=Tindex
    plot([0:0.025:1],M(k,:),'color',Color(k-jStart+1,:));
end
hold off
%xlim([0.2 0.6])
title('Accumulated mRNA over AP axis(simulation)')
xlabel('AP axis')
ylabel('Accumulated mrNA(AU)')
xlim([0.25 0.6])
%Let's see how the slope at the boundary change according to different Dm
%Dm=[0.1:0.1:10] -> Fig.4 or 5
%For example, for certain time points, compare the slopes for different Dm
%values.
