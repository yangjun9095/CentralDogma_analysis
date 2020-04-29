% Plot for accumulated mRNA/ inferred protein / measured protein
function Figure_CentralDogma_infermRNAProtein (Prefix1)
%% Load datasets
Prefix1 = '2016-10-18-Hb-nbGFP-MS2-mCherry';
%Prefix1 = '2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1';
%Prefix1 = '2018-08-09-hbP2P-MS2V5-2xMCP-mCherry-vasa-eGFP1-3';
MS2 = load(['E:\YangJoon\LivemRNA\Data\DynamicsResults\',Prefix1,'\CompiledParticles.mat']);

%%  Sketchy calibration by just comparing fluorescence (AU) with Hernan's
% 2013 paper. The AU maximum is ~1200, and the corresponding # of active
% RNAP loaded is ~100. 
% So, I will divide the fluorescence intensity (of my construct) with 12. 
% This calibration should be done more rigorously later.
RNAPLoadingRate = MS2.MeanVectorAP(MS2.nc14:end,:)/12;
SDLoadingRate = MS2.SDVectorAP(MS2.nc14:end,:)/12;
NParticles = MS2.NParticlesAP(MS2.nc14:end,:);

NBTime = MS2.ElapsedTime(MS2.nc14:end);

%% Plot RNAP Loading rate
% plot
Nt=length(NBTime); %# of time points.
% Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(NBTime); %End time (MS2.nc14:end)
    colormap(jet(256));
    cmap=colormap ;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

figure(1)
hold on
for i=2:Nt
    errorbar(0:0.025:1,RNAPLoadingRate(i,:),SDLoadingRate(i,:)./sqrt(NParticles(i,:)),'color',Color(i-iStart+1,:))
    %pause
end
hold off
title('RNAP loading rate over AP')
xlabel('AP axis')
ylabel('RNAP loading rate (number of RNAPs/ min)')
xlim([0.2 0.6])
set(gca,'fontsize',30)
colorbar
%% Calculate the accumulated mRNA
RNAPLoadingRate(isnan(RNAPLoadingRate))=0;
SDLoadingRate(isnan(SDLoadingRate))=0;
% Parameters (um for length, and minutes for time)
Tm = 60; %mRNA half-life = 60min, lambdam : degradation rate.
GammaM=log(2)./Tm; 
Dm = 0;
Dm=Dm*60; %[0.1:0.1:0.4]*60; %um^2/min. mRNA diffusion in Mouse, from BioNumbers.Yan et al., Dynamics of Translation of Single mRNA Molecules In Vivo. Cell. 2016
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
km=Dm/((dx)^2); %1/min. unit.
    
Nt=length(NBTime); %# of time points.
AccumulatedmRNA=zeros(Nt,length(RNAPLoadingRate(1,:)));
ErrorAccumulatedmRNA = zeros(Nt,length(SDLoadingRate(1,:)));

%Integrate the # of mRNA using the equation from Hernan and Jacques' paper,
AP1 = 9; % APbin start
AP2 = 21; % APbin end

%%
for i=3:Nt % Integration starts from nc 14
    dt = (NBTime(i)-NBTime(i-1));
    % For AP bins in the middle
    for j=AP1+1:41
        if j~=41
            AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+... % mRNA from the previous time point
                dt*RNAPLoadingRate(i-2,j)-... % mRNA synthesis, Et=2min, thus dNp/dt ~Fluo(t-Et/2)=Fluo(t-1min)
                GammaM*AccumulatedmRNA(i-1,j)*dt+...% mRNA degradation
                km*dt*(AccumulatedmRNA(i-1,j-1)+AccumulatedmRNA(i-1,j+1))-2*km*AccumulatedmRNA(i-1,j)*dt; %mRNA diffusion
            ErrorAccumulatedmRNA(i,j) = sqrt((ErrorAccumulatedmRNA(i-1,j).^2 + SDLoadingRate(i-2,j).^2) * dt);
        else % the last APbin (41th)
            AccumulatedmRNA(i,j) = AccumulatedmRNA(i-1,j)+...% mRNA from the previous time point
                dt*RNAPLoadingRate(i-2,j)-... % mRNA synthesis, Et=2min, thus dNp/dt ~Fluo(t-Et/2)=Fluo(t-1min)
                GammaM*AccumulatedmRNA(i-1,j)*dt+...% mRNA degradation
                km*dt*(AccumulatedmRNA(i-1,j-1))-km*dt*AccumulatedmRNA(i-1,j);
            ErrorAccumulatedmRNA(i,j) = sqrt((ErrorAccumulatedmRNA(i-1,j).^2 + SDLoadingRate(i-2,j).^2) * dt);
        end
    end
    % For the first AP bin(AP1), we assume that the diffusion to and from the
    % left is same.
        AccumulatedmRNA(i,AP1) = AccumulatedmRNA(i-1,AP1)+... % mRNA from the previous time point
            dt*RNAPLoadingRate(i-2,AP1)-... % mRNA synthesis
            GammaM*AccumulatedmRNA(i-1,AP1)*dt+...% mRNA degradation
            km*dt*(AccumulatedmRNA(i-1,AP1+1)-AccumulatedmRNA(i-1,AP1));
        ErrorAccumulatedmRNA(i,AP1) = sqrt((ErrorAccumulatedmRNA(i-1,AP1).^2 + SDLoadingRate(i-2,AP1).^2) * dt);
    % For the last AP bin, we assume that the diffusion from the right is
    % zero (so, I included AP2 into for for loop above)
    
    
end

% plot
% Color(gradation from blue to red) for each time point
iStart=2; %Start time
iEnd=length(NBTime); %End time (MS2.nc14:end)
    colormap(jet(256));
    cmap=colormap ;
    Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

figure(1)
hold on
for i=3:Nt
    errorbar(0:0.025:1,AccumulatedmRNA(i,:),ErrorAccumulatedmRNA(i,:),'color',Color(i-iStart+1,:))
end
hold off
title('Accumulated mRNA along the AP axis')
xlabel('AP axis')
ylabel('Accumulated mRNA (Number of mRNA)')
xlim([0.2 0.6])
set(gca,'fontsize',30)
colorbar

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
Tp = 50;
gammaP=log(2)/Tp; %Use Min. as time unit. I can change this later
%r_p:Translation rate for each AP bin, at nc14. 
rp=2; %protein molecules / mRNA / min. 
% I can change it later. Now this value is from Petkova et al., 2014 Curr. Bio (Gregor lab)

%Dp:diffusion coefficient of Protein,
Dp = 1;
Dp=Dp*60; %[0.1:0.1:0.4]*60; %um^2/min. estimated protein diffusion in embryo. From 7um^2/sec for Bcd by FCS measurements (Abu-Arish)
dx=500/40; %delta x = 500um/40(bins) = 12.5um/bin
kp=Dp/(dx)^2; %1/min. unit.

%NBTime;
dt = NBTime(2)-NBTime(1);
%Make a zero-matrix
Protein=zeros(length(NBTime),41);
Protein(1,:)=0;  %Protein(t=0)=0 at all APbins. Initial Condition
ErrorProtein = zeros(length(NBTime),41);
ErrorProtein(1,:) = 0;

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
for i=2:length(NBTime)  %for all Timepoints
    for j=AP1+1:length(AP) %for all APbins
        if j~=41
            Protein(i,j)=Protein(i-1,j)+...
                rp*AccumulatedmRNA(i-1,j)*dt-...
                gammaP*Protein(i-1,j)*dt+...
                kp*dt*(Protein(i-1,j-1)+Protein(i-1,j+1))-2*kp*Protein(i-1,j)*dt;
            ErrorProtein(i,j) = sqrt(ErrorProtein(i-1,j).^2 + rp*ErrorAccumulatedmRNA(i-1,j).^2*dt);
        else
            Protein(i,j)=Protein(i-1,j)+...
                rp*AccumulatedmRNA(i-1,j)*dt-...
                gammaP*Protein(i-1,j)*dt+...
                kp*dt*Protein(i-1,j-1)-kp*Protein(i-1,j)*dt;
            ErrorProtein(i,j) = sqrt(ErrorProtein(i-1,j).^2 + rp*ErrorAccumulatedmRNA(i-1,j).^2*dt);
        end
    end
    % AP1 bin, assuming the influx from the left side and outflux toward
    % the AP1 bin is the same
                Protein(i,AP1)=Protein(i-1,AP1)+...
                rp*AccumulatedmRNA(i-1,AP1)*dt-...
                gammaP*Protein(i-1,AP1)*dt+...
                kp*dt*Protein(i-1,AP1+1)-kp*dt*Protein(i-1,AP1);
            
                ErrorProtein(i,AP1) = sqrt(ErrorProtein(i-1,AP1).^2 +...
                                    rp*ErrorAccumulatedmRNA(i-1,AP1).^2*dt);
end

figure(2)
hold on
for i=3:Nt
    errorbar(0:0.025:1,Protein(i,:),ErrorProtein(i,:),'color',Color(i-iStart+1,:))
end
hold off
title('Predicted protein along the AP axis')
xlabel('AP axis')
ylabel('Predicted Protein (Number of protein molecules)')
xlim([0.2 0.6])
set(gca,'fontsize',30)
colorbar

%% Measurement of Protein signal using the Nanobody-eGFP signal
% Prefix2 = '2017-03-14-Hb-P2P-MS2V5-NB-MCP-NLS-mCherry' (For now)
% Load the dataset
Prefix2 = Prefix1;
%Prefix2 = '2018-08-08-hbP2P-MS2V5-NB-2xMCP-mCherry-vasa-eGFP1';
%NBData = load('E:\YangJoon\LivemRNA\Data\DynamicsResults\2017-12-31-Hb-P2P-MS2V5-NB-MCP-mCherry-vasa-eGFP\CompiledNuclei.mat');

NBData = load(['E:\YangJoon\LivemRNA\Data\DynamicsResults\',Prefix2,'\CompiledNuclei.mat']);
% Now, think about the protein

% The fluorescence of Hb-NB bound by eGFP inside the nucleus
NBTime = NBData.ElapsedTime(NBData.nc14:end);
%NBFluoAll = NBData.MeanVectorAll(NBData.nc14:end);
NBFluo = NBData.MeanVectorAP(NBData.nc14:end,:);
NBFluoError = NBData.SDVectorAP(NBData.nc14:end,:);
% Background fluorescence, for now, average nuclear fluo over posterior
% AP bins.
%PosteriorAPbin = 18;
%BGFluoAll = (nanmean(NBData.MeanVectorAP(NBData.nc14:end,PosteriorAPbin:end),2))'; 
%BGFluoAll(isnan(BGFluoAll))=0;
%NBProteinFluoAll = NBFluoAll - BGFluoAll;

% for i=1:length(NBTime)
%     NBProteinFluo(i,:) = NBFluo(i,:)-BGFluoAll(i);
% end

% Plot the Average mRNA and Protein over time
A = 50; % just a constant for scaling of protein

figure(3)
hold on
%plot(NewTime,AccumulatedmRNA,'LineWidth',3)
for i=3:Nt
    errorbar(0:0.025:1,NBFluo(i,:),NBFluoError(i,:),'color',Color(i-iStart+1,:))
end
xlim([0.2 0.6])
set(gca,'fontsize',30)
title('Measured protein over AP')
xlabel('Time (min)')
ylabel('Measured protein concentration (AU)')
legend('Protein (LlamaTag signal)')

%% 
hold on
for i=1:Nt-1
    plot(0:0.025:1,NBData.MeanVectorAP(NBData.nc14+i,:))
    pause
    
end