function Pipeline

%% Load the data set

mRNA=load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1\CompiledParticles.mat');
Protein=load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1\CompiledNuclei.mat');

%% Plot protein vs AP
%This plots Protein(Fluorescence) level with respect to AP axis over time.
%Each color means each time point.

%Color(gradation from blue to red) for each time point
iStart=1; %Start time
iEnd=length(Protein.ElapsedTime); %End time
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
%Plot Mean protein concentration(au) for each AP bin along time
hold on
for i=iStart:iEnd
    plot(Protein.APbinID,Protein.MeanVectorAP(i,:),'color',Color(i-iStart+1,:));
    %ylim([0,250])  %without this, the graph figure axis would chage (for making movies)
end
hold off
title('Protein Concentration along the AP axis')
xlabel('AP Position')
ylabel('Protein (MeanVectorAP) (au)')


%% Plot mRNA production(accumulation) (spot flurescence) vs AP
%This plots mRNA production rate(MS2 spot Fluorescence with respect to AP axis over time.
%Each color means each time point.
mRNApercell(1,1:41)=0;
%Color(gradation from blue to red) for each time point
iStart=1;
iEnd=length(mRNA.ElapsedTime);
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

%Plot (mean)mRNA rate for each AP bin along time
hold on
for i=iStart:iEnd
    for j=1:41
        mRNApercell(i,j)=nansum(mRNA.MeanVectorAP(1:i,j));
    plot(mRNA.APbinID,mRNApercell(i,:),'color',Color(i-iStart+1,:))
   
%     xlim([0,1])
     %ylim([0,1000])
    end
end
hold off
xlabel('AP Position')
ylabel('mRNA production rate(MeanVectorAP)')


%% mRNA & Protein vs AP position (for each time frame) - Combine above two codes
%Plot "Protein" and "mRNA production rate" with AP axis (over time) at the
%same time. Each color means each time points.
iStart=1;
iEnd=length(mRNA.ElapsedTime);
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
hold on
for i=iStart:iEnd
    plot(Protein.APbinID,Protein.MeanVectorAP(i,:),'color',Color(i-iStart+1,:))
    plot(mRNA.APbinID,mRNA.MeanVectorAP(i,:),'color',Color(i-iStart+1,:))
 
    xlabel('AP Position')
    ylabel('MeanFluoProtein')
end
hold off
%% mRNA(rate) or Protein vs time(Frame) [for each AP]

%Color(gradation from blue to red) for each AP bin
%for iStart&iEnd, we can put some specific range based on data's AP bins
iStart=1;
iEnd=length(mRNA.APbinID);
colormap(jet(256));
cmap=colormap;
Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

%Plot (used "time for x-axis" from mRNA.ElapsedTime because the dimension seems to be
%different for Protein.ElapsedTime)
%By changing from either mRNA rate or Protein, we can get mRNA rate or
%Protein dynamics graph along the time. (for different AP bins)
hold on
for i=16%iStart:iEnd
    plot(mRNA.ElapsedTime(60:end)-25,3*Protein.MeanVectorAP(60:end,i)-200,'g')%'color',Color(i-iStart+1,:))  %Protein(each AP bin) vs time
    plot(mRNA.ElapsedTime(60:end)-25,mRNA.MeanVectorAP(60:end,i),'r')%'color',Color(i-iStart+1,:))     %mRNA(rate, at each AP bin) vs time
end

%We have to change title or label according to mRNA rate or Protein.
title('\fontsize{20} Protein concentration / mRNA transcription rate over time ')
xlabel('\fontsize{18} Time(min)')
ylabel('\fontsize{18} Protein concentration / mRNA transcription rate(au)')
axis([0,50,0,700])
legend('\fontsize{18} Protein concentration','\fontsize{18} mRNA transcription rate')%'\fontsize{12}3','\fontsize{12}4',...
    %'\fontsize{12}5','\fontsize{12}6','\fontsize{12}7','\fontsize{12}8','location','northwest')
hold off

%% mRNA Accumulation over time (MeanVectorAll, for average over all AP)
%First, start with mRNA trapz for MeanVectorAPAll (Mean mRNA rate for
%all spatial window)
%This mRNA Accumulation calculation is under the assumption that there is
%no mRNA degradation.
%I assume Nan values as 0 for calculation.

%Define AccumulatedmRNA and set its initial value.
AccumulatedmRNAall(1)=0;
if isnan(mRNA.MeanVectorAll(1))
    mRNA.MeanVectorAll(1)=0;
end
%If the MeanVectorAll value for mRNA rate is Nan, I assume it to be 0.
for i=2:length(mRNA.ElapsedTime)
    if isnan(mRNA.MeanVectorAll(i))
       mRNA.MeanVectorAll(i)=0;
       AccumulatedmRNAall(i)=AccumulatedmRNAall(i-1);
    else
    AccumulatedmRNAall(i)=trapz(mRNA.ElapsedTime(1:i),mRNA.MeanVectorAll(1:i));
    end
end

%Plot mRNA accumulation along time (plus, plot Protein)
hold on
plot(mRNA.ElapsedTime,AccumulatedmRNAall,'r')   %Red means accumulated mRNA
plot(mRNA.ElapsedTime,100*Protein.MeanVectorAll,'g')    %Green means protein, x50 is for scaling for plot.
hold off

%cross-correlation for whole embryo Mean mRNA and Mean protein
%To see if there is cross-correlation between mRNA and protein in whole
%embryo level.
clear Ntotal
clear Time
clear NegTime
clear NewTime
%It should be cross-correlation between mRNA(accumulated) and rate of
%Protein production

%Protein Nan -> 0
SmoothDerivProtein(isnan(SmoothDerivProtein))=0;
N(isnan(N))=0;
%Use AccumulatedmRNA for the amount of mRNA
%Time difference (tau)
Time = Protein.ElapsedTime;
Time = Time(:,2:length(SmoothDerivProtein));
Ntotal=length(Time)+1;
NegTime = fliplr(Time)*(-1); % reverse the time series to make negative time series
NewTime = [NegTime,0, Time]; % This is time series spanning from -T:0:T

%Make a matrix of index for Time Delay (Tau)
Tau = [-1*Ntotal+1:1:Ntotal-1];

    %We will separate the time series to 3 categories, -,0,+ 
    
    for i=1:length(Tau) %for different Tau(as an index), get G(tau)
        G(i,AP)=0;
        
        if Tau(i)<0
            for n=1:Ntotal-abs(Tau(i))
                G(i,AP)=G(i,AP)+Protein.MeanVectorAll(n)*AccumulatedmRNAall(n+abs(Tau(i)));
            end
                G(i,AP)=G(i,AP) / (Ntotal-abs(Tau(i)));
            
        elseif Tau(i)>0
            for n=1:Ntotal-Tau(i)
                G(i,AP)=G(i,AP)+Protein.MeanVectorAll(n+abs(Tau(i)))*AccumulatedmRNAall(n);  %cross-correlation with mRNA decay
                %G(tau,AP)=G(tau,AP)+SmoothDerivProtein(n,AP)*AccumulatedmRNA(n+tau,AP);
            end
                G(i,AP)=G(i,AP) / (Ntotal-Tau(i));
                
        else 
            for n=1:Ntotal
                G(i,AP)=G(i,AP)+Protein.MeanVectorAll(n)*AccumulatedmRNAall(n);
            end
                G(i,AP)=G(i,AP) / Ntotal ;
        end
    end
    %Normalization of Cross-correlation
    %G(:,AP)=G(:,AP)/std(N(:,AP))/std(SmoothDerivProtein(:,AP)); %normalization
%     F=sum(SmoothDerivProtein(:,AP).*SmoothDerivProtein(:,AP)) / length(SmoothDerivProtein(:,AP));
%     H=sum(N(:,AP,1).*N(:,AP,1)) / length(N(:,AP,1));
%     G(:,AP)=G(:,AP) / (F*H);

%For colorful index for AP index
iStart=15%1;
iEnd=20;%length(mRNA.APbinID);
    colormap(jet(128));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*127)+1,:); 
%Plot along the AP axis
hold on
figure(1)
%plot()
plot(NewTime,G(:),'.')


%% mRNA Accumulation / Protein vs Time (for each APbin)
clear AccumulatedmRNA
%Future plan : I think we can get fold-change for AP axis from this data.

%Color for different AP bins
iStart=1;
iEnd=length(mRNA.APbinID);
AccumulatedmRNA(1,iStart:iEnd)=0;
    colormap(jet(128));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*127)+1,:);
%If MeanVectorAP is Nan, we assume it 0 for integration.
    for k=iStart:iEnd
    if isnan(mRNA.MeanVectorAP(1,k))
        mRNA.MeanVectorAP(1,k)=0;
    end
    end

for i=iStart:iEnd  %for color
%If MeanVectorAP is Nan, we assume it 0 for integration.
%After changing Nans to zeros, use trapz to integrate mRNA rate to get mRNA
%accumulation.
    for j=2:length(mRNA.ElapsedTime)
        if isnan(mRNA.MeanVectorAP(j,i))
        mRNA.MeanVectorAP(j,i)=0;
        AccumulatedmRNA(j,i)=AccumulatedmRNA(j-1,i);
        else
        AccumulatedmRNA(j,i)=trapz(mRNA.ElapsedTime(1:j),mRNA.MeanVectorAP(1:j,i));
        end
    
    end
    
    %Plot AccumulatedmRNA, Protein for n.c 13~end (We can change time frame)
    hold on
    for l=iStart:iEnd
    plot(mRNA.ElapsedTime(mRNA.nc13:end),AccumulatedmRNA(mRNA.nc13:end,i),'color',Color(i-iStart+1,:));  %mRNA accumulation plot
    plot(mRNA.ElapsedTime(mRNA.nc13:end),70*Protein.MeanVectorAP(mRNA.nc13:end,i),'color',Color(i-iStart+1,:));  %Protein plot, for scaling multiply 70 or sth else.
    end
    xlabel('ElapsedTime')
    ylabel('total mRNA / Protein')
hold off
end

%% mRNA accumulated at each AP bin, for each n.c. (Hernan's Curr.Bio paper Fig.3-A)
%I need whole cycle data for n.c12, 13, 14
clear Int12
clear Int13
clear Int14
%For nc12, 13, 14, calculate the accumulated mRNA at the end of each n.c.
for i=1:length(mRNA.APbinID)
    Int12(i)=trapz(mRNA.ElapsedTime(1:mRNA.nc13),mRNA.MeanVectorAP(1:mRNA.nc13,i));
    Int13(i)=trapz(mRNA.ElapsedTime(1:mRNA.nc14),mRNA.MeanVectorAP(1:mRNA.nc14,i));
    Int14(i)=trapz(mRNA.ElapsedTime(1:end),mRNA.MeanVectorAP(1:end,i));    
end
hold on
plot(mRNA.APbinID(14:20),Int12(14:20),'b')  %(14:20)
plot(mRNA.APbinID(14:20),Int13(14:20),'g')
plot(mRNA.APbinID(14:20),Int14(14:20),'r')

hold off
title('\fontsize{18} mRNA production along the AP axis')
xlabel('\fontsize{15} AP axis')
ylabel('\fontsize{15} mRNA produced per cell')
legend('\fontsize{12}nc12','\fontsize{12}nc13','\fontsize{12}nc14')

%% Protein Smoothening - Time Window
clear SmoothProtein
clear W
%We need time window that we would average the values, W
W=5; %averaging window (1Frame~0.37 sec)
SmoothProtein(1)=0;  %Initial Condition
%Color
iStart=1;
iEnd=length(mRNA.APbinID);
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
%For specific AP, i
hold on
for i=iStart:iEnd
%Average over Time Window, "W"  (it will decrease the # of rows as W-1 )
for j=1:length(mRNA.ElapsedTime)-W+1
    SmoothProtein(j,i)=sum(Protein.MeanVectorAP(j:j+W-1,i))/W; %I used sum()/W rather than nanmean because it would give weird protein derivative.
end
%Plot from the n.c13 to the end. (We can change the time frame)
plot(mRNA.ElapsedTime(mRNA.nc13:end-W+1),SmoothProtein(mRNA.nc13:end,i),'color',Color(i-iStart+1,:));
end
hold off
title('\fontsize{18}Protein concentration over time (Smoothened)')
xlabel('\fontsize{15}Time(min)')
ylabel('\fontsize{15}Protein concentration(au)')
%% Protein Derivative (using previous SmoothProtein) for each AP
clear DerivProtein
clear DeltaP
clear DeltaT
%For colorful index for AP index
iStart=1;
iEnd=length(mRNA.APbinID);
    colormap(jet(128));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*127)+1,:);
    
% deltaP, deltaT 
DeltaP=diff(SmoothProtein);
DeltaT=diff(mRNA.ElapsedTime);

hold on
%Calculate the derivative of the Protein w.r.t. time
for i=iStart:iEnd
    for j=1:length(SmoothProtein)-1
        DerivProtein(j,i)=DeltaP(j,i)/DeltaT(j); 
    end
    plot(mRNA.ElapsedTime(mRNA.nc14:end-W),DerivProtein(mRNA.nc14:end,i),'color',Color(i-iStart+1,:))
end
hold off
title('\fontsize{18}Protein derivative over time')
xlabel('\fontsize{15}Time(min)')
ylabel('\fontsize{15}Protein derivative(au)')
%% Smoothening the Protein Derivative (+protein degradation or bleaching)
clear SmoothDerivProtein
clear W1
%We need time window that we would average the values, W1
W1=10;
SmoothDerivProtein(1,1:41)=0;
%Color
iStart=1;
iEnd=length(mRNA.APbinID);
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:); 
%For specific AP, i
hold on
for i=iStart:iEnd
%Average over Time Window, "W1"
for j=1:length(DerivProtein)-W1+1
    SmoothDerivProtein(j,i)=sum(DerivProtein(j:j+W1-1,i))/W1;
end
plot(mRNA.ElapsedTime(mRNA.nc14:end-W-W1+1),SmoothDerivProtein(mRNA.nc14:end,i),'color',Color(i-iStart+1,:));
end
hold off
title('\fontsize{18}Protein derivative over time(Smoothened)')
xlabel('\fontsize{15}Time(min)')
ylabel('\fontsize{15}Protein derivative(au)')
%% Protein Derivative vs mRNA (Accumulated, w/o degradation)

%We should run the SmoothRrotein and AccumulatedmRNA before we run this.
%We should set the length of two matrixes same, so that we can plot
%DerivProtein and AccumulatedmRNA together.

%Color
iStart=1;
iEnd=length(mRNA.APbinID);
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
hold on
for i=iStart:iEnd
    plot(AccumulatedmRNA(1:length(SmoothDerivProtein),i),SmoothDerivProtein(:,i),'color',Color(i-iStart+1,:));
end

hold off
title('\fontsize{18}Protein derivative vs total mRNA')
xlabel('\fontsize{15}mRNA Accumulation(w/o degradation)')
ylabel('\fontsize{15}Protein derivative(au)')

%% mRNA accumulation vs time, for each AP bin. (with Degradation, half life of mRNA)
clear Time
clear DeltaT
clear T1
clear T2
clear N
clear lambdam
% I will add Protein derivative O.D.E based on mRNA accumulation simulation
% I would need O.D.E to calculate this.
% Half-life of mRNA is depicted as lambdam
% We have to define Time, DeltaT, N(total mRNA, at specific AP position),
% so that we can change time range (n.c) as we want.i.e.
% ex) Time=mRNA.ElapsedTime(mRNA.nc13:mRNA:nc14) for nc13
Time=mRNA.ElapsedTime;
DeltaT=diff(Time); %for integration(with dt)
%Integration Time range
T1=1;
T2=length(Time);
%Define decay constant using half-life of mRNA
lambdam=log(2)./[0:10:300];  %0~300 : half-life of mRNA (min)
N(1,length(mRNA.APbinID),length(lambdam))=0; %Assuming N(1)=0, initial condition
% Coloring for different AP bins.
iStart=14%1;
iEnd=21;%length(mRNA.APbinID);
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:); 
%Plot new figure for each lambdam
for k=1:length(lambdam)
    figure(k)
hold on
    for i=iStart:iEnd
        for j=T1:T2-1
            %mRNA(t+dt)=mRNA(t)+mRNArate*dt - lambdam*mRNA(t)*dt
            N(j+1,i,k)= N(j,i,k)+ mean(mRNA.MeanVectorAP(j:j+1,i))*DeltaT(j) - lambdam(k)*N(j,i,k)*DeltaT(j);
        end
    plot(Time,N(:,i,k),'color',Color(i-iStart+1,:))
    end
title('\fontsize{18}total mRNA (w/ degradation)')
xlabel('\fontsize{15}Time(min)')
ylabel('\fontsize{15}mRNA Accumulation(au)')
legend('\fontsize{12} 1','\fontsize{12}2','\fontsize{12}3','\fontsize{12}4',...
    '\fontsize{12}5','\fontsize{12}6','\fontsize{12}7','\fontsize{12}8','location','northwest')
hold off
end

%% Time Delay between transcription and translation (T)
clear DelayedSmoothDerivProtein
%There is time delay between transcription(mRNA) and translation(Protein)
%ProteinDerivative(t)=k*mRNA(t-T), T= delay time
T=[0:5:100];   %T=0~40 (frames, ~0.37min/frame)
DelayedSmoothDerivProtein(1:length(SmoothDerivProtein),1:41,1:length(T))=0;
%color
iStart=1;
iEnd=length(mRNA.APbinID);  %We need to change these values to get diff. colors.
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
    
%for l=1:length(T)
     for k=1:length(lambdam)
     figure(k)
    %figure(l)
    hold on
    for i=iStart:iEnd
        %We have to choose delay time from T vector, and put the index as
        %j=T(index)+1:length(SmoothDerivProtein)
        for j=T(5)+1:length(SmoothDerivProtein)
            DelayedSmoothDerivProtein(j,i,5)=SmoothDerivProtein(j-T(5),i);
        end
    plot(N(1:length(DelayedSmoothDerivProtein),i,k),DelayedSmoothDerivProtein(:,i,5),'color',Color(i-iStart+1,:));%We should set k th value for the N(j,i,k), which indicates mRNA half-life
    end
%     end
title('\fontsize{18}protein derivative vs total mRNA (w/ degradation)')
xlabel('\fontsize{15}mRNA(t)(au)')
ylabel('\fontsize{15}protein derivative(t+T) (au)')
legend('\fontsize{12} 1','\fontsize{12}2','\fontsize{12}3','\fontsize{12}4',...
    '\fontsize{12}5','\fontsize{12}6','\fontsize{12}7','\fontsize{12}8','location','northwest')
    hold off
end
%Should be editted more.....

%% Protein Derivative(SmoothDerivProtein / with Degradation, lambdap) vs mRNA accumulation (N,integrtion of real data w/ degradation)
%dp(t+T)/dt + lambdap*p(t+T)=k*mRNA(t);
clear lambdap
clear Lefthand
clear Righthand
%We should run the DelayedSmoothDerivProtein and N before we run this.
%We have to pick one AP bin, which is best among one data set. -> needs to
%be editted so that it can plot multiple AP bins.

APbin=14;   %set one AP bin  / we can make for loop for AP bins.
lambdap=log(2)./[0:10:300];  %half-life of protein. 0~100 (min)
T=[0:5:25];   %Delay time for mRNA->Protein , T=0~10 (min) we have to
% multiply 0.37 to frames.

%Color for each time points (j)
jStart=1;
jEnd=length(SmoothDerivProtein); %Time points
    colormap(jet(256));
    cmap=colormap;
    Color=cmap(round(((jStart:jEnd)-jStart)/(jEnd-jStart)*255)+1,:);
    
%Because there are too many variables, I just choose one APbin, and one
%mRNA half-life as one index in N(j,k,i) previously defined.
%So, now we can change delay time, protein half-life, and see how dp/dt vs
%mRNA changes along the time. I think we can define multi-dimensional
%variable which contains mRNA,protein half-life, delay time, AP bins, n.c.
%and by changing plot index, we can compare those factors' effect.

for l=1:length(T)   %for different delay time
    for k=1:length(lambdap)  %for different protein-half-life
        for j=jStart:jEnd-T(l)  %for time points (Color)
                Lefthand(j,k,l)= SmoothDerivProtein(j+T(l),APbin)+lambdap(k)*SmoothProtein(j+T(l),APbin); %dp/dt(t+T) + lambdap*p(t+T) 
                %We can set multiple APbin using for loop
                Righthand(j,k,l)=N(j,APbin,6);  %N(j,k,i) ->we should set i th value for lambdam manually %6th value is 60min
        end
    end
end
%Plot
for n=1:length(T)  %delay time (0~10min)
    figure(n)
    hold on
    % for m=1:length(lambdap)
        for i=mRNA.nc13:length(Lefthand)  %i : time!, we can put mRNA.nc13, nc14, for example.
            plot(Righthand(i,20,n),Lefthand(i,20,n),'Marker','o','color',Color(i-iStart+1,:))
        end
    % end
    title('\fontsize{18}protein derivative vs total mRNA (w/ degradation)')
    xlabel('\fontsize{15}mRNA(t)(au)')
    ylabel('\fontsize{15}protein derivative(t+T) (au)')
    hold off
    
end

end