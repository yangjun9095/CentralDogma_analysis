%% Calculation of mRNA and Protein from the Dataset (with calibration)
function [AccumulatedmRNA,Protein,RNAPLoadingRate3,InferredProtein, NBProteinFluo] = ...
                CorrelatemRNAProtein_function(Prefix1,Prefix2,Dm,Dp,Tm,Tp,rp)

% This code is for the calculation of the mRNA and protein from the Imaging
% data with the proper calibration, and parameters, such as production
% rate, degradation rate, and diffusion (But I will include this later)

% First, I need to calibrate the MS2 signals from the MCP-mCherry (No NLS,
% NLS, and tdMCP) to the number of mRNA produced per given time.
% Note. I need to check if MCP-mCherry construct is saturating the MS2
% loops, by comparing with higher dosages. If it saturates, then I can use
% the smFISH data for calibration. 
% I think I need to do the smFISH for the VK33 line P2P-MS2.V5-Hunchback constructs.
% Use the OneNote that Hernan wrote about calibration for now. All I need to do is
% taking a couple of datasets of P2P.V1-MS2 x MCP-mCherry, 
% and the number of mRNA from HG's calculation.

% Input arguments :
% Prefix1 : MS2 dataset
% Prefix2 : NB dataset
% Dm, Dp : Diffusion coefficient of mRNA and protein respectively
% (um^2/sec)
% Tm, Tp : half-life of mRNA and protein respectively (min)
% Note. Units for space and time in this script is um and minutes.

% Use the dataset of hb P2P-MS2.V5-Hunchback (no nanobody) x MCP-GFP for mRNA profile.
% For now, I will use Prefix1= '2017-03-02-Hb-P2P-MS2V5-MCP-GFP'
% Load the dataset

MS2 = load(['E:\YangJoon\LivemRNA\Data\DynamicsResults\',Prefix1,'\CompiledParticles.mat']);

%% Calibrate the MS2 signals to number of mRNAs
% I'm basing on Hernan's paper (Garcia, 2013), using smFISH to calibrate
% the P2P.MS2.V1 signals
% hold on 
% title('Mean Spot Fluorescence over Time')
% xlabel('Time (min)')
% ylabel('Mean Spot Fluorescence (AU)')
% for AP=1:41
%     plot(MS2.ElapsedTime,MS2.MeanVectorAP(:,AP))
%     pause
% end

%%  Sketchy calibration by just comparing fluorescence (AU) with Hernan's
% 2013 paper. The AU maximum is ~1200, and the corresponding # of active
% RNAP loaded is ~100. 
% So, I will divide the fluorescence intensity (of my construct) with 12. 
% This calibration should be done more rigorously later.
RNAPLoadingRate = MS2.MeanVectorAP/12;

%Plot the RNAP loading rate along time (over AP)
% for i=1:length(RNAPLoadingRate3(51:100:end,:))
%     plot(0:0.025:1,RNAPLoadingRate3(51+100*(i-1),:))
%     pause
% end
% Here, I made assumptions 1) hb P2P-MS2-lacZ in 38F1 expression is same as 
% hb P2P-MS2.V5-Hb in VK33. 
%% First, plot the Averaged Spot fluo trace
ap = 12;
%Use MeanVectorAP at APbinID = ap, and plot them for timecourse 
hold on
plot(MS2.ElapsedTime,MS2.MeanVectorAP(:,ap),'g','LineWidth',5)
%plot(mRNA.ElapsedTime(mRNA.nc13:end)-(mRNA.nc13-mRNA.nc13),mRNA.MeanVectorAll(mRNA.nc13:end),'b')
hold off

title('MS2 spot fluorescence')
xlabel('Time(min)')
ylabel('Spot Fluorescence (AU)')
legend('MCP-GFP')


%% Plot the Offset
% hold on
% %plot(mRNA.ElapsedTime,mRNA.MeanVectorAll,'b')
% plot(MS2.ElapsedTime,MS2.MeanOffsetVector,'g')
% hold off
% 
% title('Spot Fluorescence and Offset')
% legend('MS2 spot fluorscence')%,'Offset')
% xlabel('Time(min)')
% ylabel('Fluorescence (AU)')
% set(gca,'fontsize',30)

% %% Interpolation of Spot fluorescence (Spatial)
% % This is to make the APbin size smaller, for proper numerical calculation
% % of Reaction-Diffusion equation.
% % Define a new x-axis
% xx = 0:1.25:500;
% RNAPLoadingRate2 = zeros(length(RNAPLoadingRate),length(xx));
% for i=1:length(RNAPLoadingRate)
%     RNAPLoadingRate2(i,:) = pchip(0:12.5:500,RNAPLoadingRate(i,:),xx);
% end 
% 
% %% Check for spatial interpolation
% % Plot the AP profile of RNAP loading rate at each time point
% for j=1:length(RNAPLoadingRate)
%     clf
%     hold on
%     plot(0:0.025:1,RNAPLoadingRate(j,:)) % Data
%     plot(0:0.0025:1,RNAPLoadingRate2(j,:)) % Interpolation
%     pause
% end
% 
% % The result seems to be a good match.

%% Only think about nc13 and nc14
% Time = MS2.ElapsedTime;
% nc12Frame = 6; % the frame when nc12 ends
% Time = Time(6:end) - Time(6); % Cut out before the nc13, and sync it to 0 min
% RNAPLoadingRate = RNAPLoadingRate(6:end,:);
%% Interpolation of Spot Fluorescence (Temporal)
% This is to make the time bin smaller, for proper numerical calculation
% of Reaction-Diffusion equation.
% Define a new time-axis
NBTime = MS2.ElapsedTime;
NewTime =0:0.01:max(NBTime) ;% now the increment is 0.01 (min), previously, it was around 0.68 (min)
RNAPLoadingRate(isnan(RNAPLoadingRate))=0;

RNAPLoadingRate3 = zeros(length(NewTime),41);
for i=1:length(RNAPLoadingRate(1,:))
    RNAPLoadingRate3(:,i) = pchip(NBTime,RNAPLoadingRate(:,i),NewTime);
end

% ignore nc12,
RNAPLoadingRate3(1:600,:)=0;

%% Check the temporal interpolation by plot
% AP1 = 9; % APbin start
% AP2 = 21; % APbin end
% for j=AP1:AP2
%     clf
%     hold on
%     plot(NBTime,RNAPLoadingRate(:,j),'o') % Data
%     plot(NewTime,RNAPLoadingRate3(:,j)) % Interpolation
%     pause
% end

% %% mRNA Accumulation from Spot fluorescence (using SDD model, from the MS2 spot data)
% clear AccumulatedmRNA
% % First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
% %Hb gene length : 3232 bp,
% %r_elongation : 1.54kb/min
% %Et = 2min(~4 frames in ElapsedTime)
% 
% % Big assumption here is that I have only 0.2~0.6 AP bins, so I assume that
% % 0.0~0.2 has has almost same mRNA production as 0.2 AP bin, and 0.6~1.0
% % has almost no mRNA production. ( I should acquire data for these AP bins
% % later)
% 
% % Parameters (um for length, and minutes for time)
% gammaM = 60; %mRNA half-life = 60min, lambdam : degradation rate.
% lambdam=log(2)./gammaM; 
% Dm=1*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016
% dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
% km=Dm/((dx)^2); %1/min. unit.
% 
% RNAPLoadingRate=MS2.MeanVectorAP / 12; %This is sketchy calibration of MS2 signal to # of active RNAP .*mRNA.OnRatioAP;
% RNAPLoadingRate(isnan(RNAPLoadingRate))=0;
% Nt=length(MS2.ElapsedTime); %# of time points.
% AccumulatedmRNA=zeros(Nt-4,41);
% 
% %Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
% AP1 = 9; % APbin start
% AP2 = 21; % APbin end
% for i=4:Nt
%     dt = (MS2.ElapsedTime(i)-MS2.ElapsedTime(i-1));
%     % For AP bins in the middle
%     for j=AP1+1:41
%         if j~=41
%             AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+... % mRNA from the previous time point
%                 dt*RNAPLoadingRate(i,j)-... % mRNA synthesis
%                 lambdam*AccumulatedmRNA(i-1,j)*dt+...% mRNA degradation
%                 km*dt*(AccumulatedmRNA(i-1,j-1)+AccumulatedmRNA(i-1,j+1))-2*km*AccumulatedmRNA(i-1,j)*dt;
%         else
%             AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+...% mRNA from the previous time point
%                 dt*RNAPLoadingRate(i,j)-... % mRNA synthesis
%                 lambdam*AccumulatedmRNA(i-1,j)*dt+...% mRNA degradation
%                 km*dt*(AccumulatedmRNA(i-1,j-1))-km*dt*AccumulatedmRNA(i-1,j);
%         end
%     end
%     % For the first AP bin, we assume that the diffusion to and from the left is same.
%         AccumulatedmRNA(i,AP1) = AccumulatedmRNA(i-1,AP1)+... % mRNA from the previous time point
%             dt*RNAPLoadingRate(i,AP1)-... % mRNA synthesis
%             lambdam*AccumulatedmRNA(i-1,AP1)*dt+...% mRNA degradation
%             km*dt*(AccumulatedmRNA(i-1,j+1)-AccumulatedmRNA(i-1,j));
%     % For the last AP bin, we assume that the diffusion from the right is
%     % zero (so, I included AP2 into for for loop above)
%     
%     
% end
% 
% % plot
% % Color(gradation from blue to red) for each time point
% iStart=2; %Start time
% iEnd=length(MS2.ElapsedTime); %End time
%     colormap(jet(256));
%     cmap=colormap ;
%     Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
% hold on
% for i=2:Nt-4
%     plot(0:0.025:1,AccumulatedmRNA(i,:),'color',Color(i-iStart+1,:))
% end
% hold off
% title('Accumulated mRNA along the AP axis')
% xlabel('AP axis')
% ylabel('Accumulated mRNA (Number of mRNA)')
% xlim([0.2 0.6])
% set(gca,'fontsize',30)
% colorbar

%% mRNA Accumulation from Spot fluorescence (From the Data, temporally interpolated)
clear AccumulatedmRNA
% First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
%Hb gene length : 3232 bp,
%r_elongation : 1.54kb/min
%Et = 2min(~3 frames in ElapsedTime)

% Big assumption here is that I have only 0.2~0.6 AP bins, so I assume that
% 0.0~0.2 has has almost same mRNA production as 0.2 AP bin, and 0.6~1.0
% has almost no mRNA production. ( I should acquire data for these AP bins
% later)

% Parameters (um for length, and minutes for time)
%Tm = 60; %mRNA half-life = 60min, lambdam : degradation rate.
GammaM=log(2)./Tm; 
Dm=Dm*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
km=Dm/((dx)^2); %1/min. unit.

%RNAPLoadingRate=MS2.MeanVectorAP / 12; %This is sketchy calibration of MS2 signal to # of active RNAP .*mRNA.OnRatioAP;
RNAPLoadingRate3(isnan(RNAPLoadingRate3))=0;
Nt=length(NewTime); %# of time points.
AccumulatedmRNA=zeros(Nt,length(RNAPLoadingRate3(1,:)));

%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
AP1 = 9; % APbin start
AP2 = 21; % APbin end


for i=600:Nt % Integration starts from nc 13
    dt = (NewTime(i)-NewTime(i-1));
    % For AP bins in the middle
    for j=AP1+1:41
        if j~=41
            AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+... % mRNA from the previous time point
                dt*RNAPLoadingRate3(i-100,j)-... % mRNA synthesis, Et=2min, thus dNp/dt ~Fluo(t-Et/2)=Fluo(t-1min)
                GammaM*AccumulatedmRNA(i-1,j)*dt+...% mRNA degradation
                km*dt*(AccumulatedmRNA(i-1,j-1)+AccumulatedmRNA(i-1,j+1))-2*km*AccumulatedmRNA(i-1,j)*dt; %mRNA diffusion
        else % the last APbin (41th)
            AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+...% mRNA from the previous time point
                dt*RNAPLoadingRate3(i-100,j)-... % mRNA synthesis, Et=2min, thus dNp/dt ~Fluo(t-Et/2)=Fluo(t-1min)
                GammaM*AccumulatedmRNA(i-1,j)*dt+...% mRNA degradation
                km*dt*(AccumulatedmRNA(i-1,j-1))-km*dt*AccumulatedmRNA(i-1,j);
        end
    end
    % For the first AP bin(AP1), we assume that the diffusion to and from the
    % left is same.
        AccumulatedmRNA(i,AP1) = AccumulatedmRNA(i-1,AP1)+... % mRNA from the previous time point
            dt*RNAPLoadingRate3(i-100,AP1)-... % mRNA synthesis
            GammaM*AccumulatedmRNA(i-1,AP1)*dt+...% mRNA degradation
            km*dt*(AccumulatedmRNA(i-1,AP1+1)-AccumulatedmRNA(i-1,AP1));
    % For the last AP bin, we assume that the diffusion from the right is
    % zero (so, I included AP2 into for for loop above)
    
    
end

% plot
% Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(NewTime); %End time
    colormap(jet(256));
    cmap=colormap ;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

figure(1)
hold on
for i=51:50:Nt
    plot(0:0.025:1,AccumulatedmRNA(i,:),'color',Color(i-iStart+1,:))
end
hold off
title('Accumulated mRNA along the AP axis')
xlabel('AP axis')
ylabel('Accumulated mRNA (Number of mRNA)')
xlim([0.2 0.6])
set(gca,'fontsize',30)
colorbar
%% 
% %% mRNA Accumulation from Spot fluorescence (From the Data, spatially, temporally interpolated)
% clear AccumulatedmRNA
% % First, I need to use Hernan and Jacques' 2014 PNAS paper's method.
% %Hb gene length : 3232 bp,
% %r_elongation : 1.54kb/min
% %Et = 2min(~4 frames in ElapsedTime)
% 
% % Big assumption here is that I have only 0.2~0.6 AP bins, so I assume that
% % 0.0~0.2 has has almost same mRNA production as 0.2 AP bin, and 0.6~1.0
% % has almost no mRNA production. ( I should acquire data for these AP bins
% % later)
% 
% % Parameters (um for length, and minutes for time)
% lambdam=log(2)./60; %mRNA half-life = 60min, lambdam : degradation rate.
% Dm=0.1*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016
% dx=500/400; %delta x = 500um/40(bins) = 12.5um/bin
% km=Dm/((dx)^2); %1/min. unit.
% 
% %RNAPLoadingRate=MS2.MeanVectorAP / 12; %This is sketchy calibration of MS2 signal to # of active RNAP .*mRNA.OnRatioAP;
% RNAPLoadingRate2(isnan(RNAPLoadingRate2))=0;
% Nt=length(MS2.ElapsedTime); %# of time points.
% AccumulatedmRNA=zeros(Nt-4,length(RNAPLoadingRate2));
% 
% %Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
% AP1 = 9; % APbin start
% InterpAP1 = (AP1-1)*10+1;
% AP2 = 21; % APbin end
% InterpAP2 = (AP2-1)*10+1;
% 
% for i=2:Nt-4
%     dt = (MS2.ElapsedTime(i)-MS2.ElapsedTime(i-1));
%     % For AP bins in the middle
%     for j=InterpAP1+1:InterpAP2
%         AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+... % mRNA from the previous time point
%             dt*RNAPLoadingRate2(i,j)-... % mRNA synthesis
%             lambdam*AccumulatedmRNA(i-1,j)*dt+...% mRNA degradation
%             km*dt*(AccumulatedmRNA(i-1,j-1)+AccumulatedmRNA(i-1,j+1))-2*km*AccumulatedmRNA(i-1,j)*dt;
%     end
%     % For the first AP bin, we assume that the diffusion to and from the
%     % left is same.
%         AccumulatedmRNA(i,InterpAP1) = AccumulatedmRNA(i-1,InterpAP1)+... % mRNA from the previous time point
%             dt*RNAPLoadingRate2(i,InterpAP1)-... % mRNA synthesis
%             lambdam*AccumulatedmRNA(i-1,InterpAP1)*dt+...% mRNA degradation
%             km*dt*(AccumulatedmRNA(i-1,InterpAP1+1)-AccumulatedmRNA(i-1,InterpAP1));
%     % For the last AP bin, we assume that the diffusion from the right is
%     % zero (so, I included AP2 into for for loop above)
%     
%     
% end
% 
% % plot
% % Color(gradation from blue to red) for each time point
% iStart=2; %Start time
% iEnd=length(MS2.ElapsedTime); %End time
%     colormap(jet(256));
%     cmap=colormap ;
%     Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);
% hold on
% for i=2:Nt-4
%     plot(0:0.0025:1,AccumulatedmRNA(i,:),'color',Color(i-iStart+1,:))
% end
% hold off
% title('Accumulated mRNA along the AP axis')
% xlabel('AP axis')
% ylabel('Accumulated mRNA (Number of mRNA)')
% xlim([0.2 0.6])
% set(gca,'fontsize',30)
% colorbar

%% Calculation of protein from the AccumulatedmRNA
% Make sure to check the dimension of AccumulatedmRNA (depending on
% which interpolation I did)

% Note that I didn't consider the delay between the transcription and
% translation. Transcripts should be shuttled to the cytoplasm, then
% translated by the ribosome. Thus, there could be delay, which I'm not
% sure how to address now.
clear Protein
%First, I need to define variables and parameters.
%gamma_P, r_p, Dp(kp)
%gamma_P : protein degradation rate, log(2)/half-life of protein (not
%known). We can use Bcd lifetime in here for now. (~50 min)
gammaP=log(2)/Tp; %Use Min. as time unit. I can change this later
%r_p:Translation rate for each AP bin, at nc14. 
%rp=2; %protein molecules / mRNA / min. 
% I can change it later. Now this value is from Petkova et al., 2014 Curr. Bio (Gregor lab)

%Dp:diffusion coefficient of Protein,
Dp=Dp*60; %[0.1:0.1:0.4]*60; %um^2/min. estimated protein diffusion in embryo. From 7um^2/sec for Bcd by FCS measurements (Abu-Arish)
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
kp=Dp/(dx)^2; %1/min. unit.

NewTime;
dt = NewTime(2)-NewTime(1);
%Make a zero-matrix
Protein=zeros(length(NewTime),41);
Protein(1,:)=0;  %Protein(t=0)=0 at all APbins. Initial Condition

%Additional Note.
%Think about the first and the last AP bins, as separately, 
%because they have diffusion to only one side(left or right).
%P(i,1)=P(i-1,1)+rp(1)*dt-gammaP*P(i-1,1)*dt+k*P(i-1,2)*dt-k*P(i-1,1)*dt;
%P(i,41)=...
%But, right now, we don't have to worry too much, because we don't care
%about boundary now.

%Make a Protein reaction-diffusion equation
%I will assume that the translation occurs during the whole time( nc13 and
%nc14)
AP=1:41;
% I will use the same argument as I calculated the mRNA,
% that 
for i=2:length(NewTime)  %for all Timepoints
    for j=AP1+1:length(AP) %for all APbins
        if j~=41
            Protein(i,j)=Protein(i-1,j)+...
                rp*AccumulatedmRNA(i-1,j)*dt-...
                gammaP*Protein(i-1,j)*dt+...
                kp*dt*(Protein(i-1,j-1)+Protein(i-1,j+1))-2*kp*Protein(i-1,j)*dt;
        else
            Protein(i,j)=Protein(i-1,j)+...
                rp*AccumulatedmRNA(i-1,j)*dt-...
                gammaP*Protein(i-1,j)*dt+...
                kp*dt*Protein(i-1,j-1)-kp*Protein(i-1,j)*dt;
        end
    end
    % AP1 bin, assuming the influx from the left side and outflux toward
    % the AP1 bin is the same
                Protein(i,AP1)=Protein(i-1,AP1)+...
                rp*AccumulatedmRNA(i-1,AP1)*dt-...
                gammaP*Protein(i-1,AP1)*dt+...
                kp*dt*Protein(i-1,AP1+1)-kp*dt*Protein(i-1,AP1);
end

%% color plot for Protein : To see how the shape changes along the time. 
% hold on
% for k=51:100:length(NewTime)
%     plot([0:0.025:1],Protein(k,:));
% end
% hold off

%% save the result variables
% save(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-03-02-Hb-P2P-MS2V5-MCP-GFP\InterredmRNAProtein.mat'],...
%     'RNAPLoadingRate','RNAPLoadingRate3','Dm','Dp','dt','dx','gammaM','gammaP','rp',...
%     'AP1','AP2','AccumulatedmRNA','Protein')

%% Measurement of Protein signal using the Nanobody-eGFP signal
% Prefix2 = '2017-03-14-Hb-P2P-MS2V5-NB-MCP-NLS-mCherry' (For now)
% Load the dataset
NBData = load(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\',Prefix2,'\CompiledNuclei.mat'])

% Now, think about the protein

% The fluorescence of Hb-NB bound by eGFP inside the nucleus
NBTime = NBData.ElapsedTime;
NBFluoAll = NBData.MeanVectorAll;
NBFluo = NBData.MeanVectorAP;
% Background fluorescence, for now, average nuclear fluo over posterior
% AP bins.
BGFluoAll = (nanmean(NBData.MeanVectorAP(:,23:end),2))'; 

NBProteinFluoAll = NBFluoAll - BGFluoAll;
for i=1:length(NBFluo)
    NBProteinFluo(i,:) = NBFluo(i,:)-BGFluoAll(i);
end
% Plot the Average mRNA and Protein over time
A = 50; % just a constant for scaling of protein

% figure(2)
% hold on
% %plot(NewTime,AccumulatedmRNA,'LineWidth',3)
% plot(NBTime,A*NBProteinFluoAll,'LineWidth',3)
% 
% title('Averaged protein level over time')
% xlabel('Time (min)')
% ylabel('Fluorescence (AU)')
% legend('Protein (NB signal)')

%% Compare the inferred Protein with the measured Protein
% Compare Protein (inferred) and NBFluo (nanobody-eGFP signal)
% 1) need to pick AP bin 
AP = 15;
% 2) synchronize the two datasets
% Nanobody data nc
nc13_NB = NBData.nc13;
nc14_NB = NBData.nc14;
% MS2 data nc (inferred protein)
nc13_Protein = MS2.nc13;
nc14_Protein = MS2.nc14;

hold on
plot(NewTime,AccumulatedmRNA)
plot(NBTime,10*NBProteinFluo)

% From these two datasets, there is a 3 min-gap between NB dataset and MS2
% dataset. NB dataset's 24 min (beginning of nc 14) corresponds to MS2
% dataset's 27 min
Tgap = 3; % min
%% Compare the NB-protein data and inferred protein data
clear time
MS2TimeIndex = zeros(length(NBTime),1);
for i=36:length(NBTime)
    time(i) = NBTime(i)+Tgap; %(3 min gap, this should be changed for different datasets)
    time(i) = round(time(i),2); % round off at 10^(-2) to get the corresponding time point
    %[row,col] = find(NewTime==time(i));
    %MS2TimeIndex(i)=row;
    MS2TimeIndex(i)=int16(time(i)*100);
    InferredProtein(i,:) = Protein(MS2TimeIndex(i),:);
end

%% Comparing the Inferred protein vs measured protein
clf
% color code for the time
iStart=33; %Start time
iEnd=length(NBTime); %End time
colormap(jet(256));
cmap=colormap ;
Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

hold on
for i=33:15:length(NBTime)

    for j=1:41
        plot(InferredProtein(i,j),NBProteinFluo(i,j),'o','color',Color(i-iStart+1,:))
    end
        
        title(['Protein : Prediction vs Measurement ','@',num2str(NBTime(i)-NBTime(nc14_NB)),'min'])
        xlabel('Inferred Protein (AU)')
        ylabel('Measured Protein (AU')
        set(gca,'fontsize',15)
        pause
end
end