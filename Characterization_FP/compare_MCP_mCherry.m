%% Compare MCP-mCherry
% This code is to figure out whether MCP-NoNLS-mCherry is at the
% saturation concentration (for MS2 loops).
% The idea is that the fluorescence accumulation over nc 13 at one AP bin
% should be the same for MCP-NoNLS-mCherry, and MCP-NLS-mCherry

%% Load the Datasets
%clear all

DataNLS=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-01-02-Hb-P2P-MS2-NLSmCherry\CompiledParticles.mat')
DataNoNLS=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-01-04-Hb-P2P-MS2-mCherry\CompiledParticles.mat')
DataNoNLS2=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-10-16-P2P-MCP-NoNLS-mCherry\CompiledParticles.mat')
%DataCTL=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-01-08-Hb-P2P-HisRFP-MCP-GFP\CompiledParticles.mat')
DataNoNLSDouble1=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-09-10-P2P-MCP-NoNLS-mCherry-doubledosage\CompiledParticles.mat')
%DataNoNLSDouble2=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-09-14-P2P-MCP-NoNLS-mCherry-doubledosage2\CompiledParticles.mat')
DataNoNLSDouble2=load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-12-13-P2P-MCP-NoNLS-mCherry-doubledosage\CompiledParticles.mat')
DataNoNLS_15 = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-09-10-P2P-MCP-NoNLS-mCherry-threeovertwodosage\CompiledParticles.mat')
DataNoNLS_15_2 = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-09-15-P2P-MCP-NoNLS-mCherry-threeovertwodosage\CompiledParticles.mat')
DataNoNLS_15_3 = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-09-17-P2P-MCP-NoNLS-mCherry-threeovertwodosage\CompiledParticles.mat')
%% mRNA accumulation
% Next, I need to accumulate the MS2 spot fluorescence at specific AP bin
% over the nc 13. Note that here I ignore the degradation.

%% 1.MCP-NLS-mCherry
NLSTime=DataNLS.ElapsedTime;
NLSFluo=DataNLS.MeanVectorAP;
NLSAPbinID=DataNLS.APbinID;
%plot(NLSTime,NLSSpot)


%Calculate the Accumulation of mRNA

% Choose a specific AP bin,
% AP = 14;  % AP/0.025 + 1 = position in AP

NLSt13=DataNLS.nc13;
NLSt14=DataNLS.nc14;

for AP=1:41
    clear MS2SpotFluo
    MS2SpotFluo=NLSFluo(:,AP);
    MS2SpotFluo(isnan(MS2SpotFluo))=0;
    
    Sum_NLS_mCherry_nc13(AP) = trapz(NLSTime(NLSt13:NLSt14),MS2SpotFluo(NLSt13:NLSt14));
    Sum_NLS_mCherry_nc14(AP) = trapz(NLSTime(NLSt14:end),MS2SpotFluo(NLSt14:end));
end


%% 2.MCP-NoNLS-mCherry
NoNLSTime=DataNoNLS.ElapsedTime;
NoNLSFluo=DataNoNLS.MeanVectorAP;
NoNLSAPbinID=DataNoNLS.APbinID;
%plot(NLSTime,NLSSpot)


%Calculate the Accumulation of mRNA
%clear MS2SpotFluo

% Choose a specific AP bin,
%AP = 14;  % AP/0.025 + 1 = position in AP

NoNLSt13=DataNoNLS.nc13;
NoNLSt14=DataNoNLS.nc14;

% For all AP bins,
for AP=1:41
    clear MS2SpotFluo
    MS2SpotFluo=NoNLSFluo(:,AP);
    MS2SpotFluo(isnan(MS2SpotFluo))=0;
    
    Sum_NoNLS_mCherry_nc13(AP) = trapz(NoNLSTime(NoNLSt13:NoNLSt14),MS2SpotFluo(NoNLSt13:NoNLSt14));
    Sum_NoNLS_mCherry_nc14(AP) = trapz(NoNLSTime(NoNLSt14:end),MS2SpotFluo(NoNLSt14:end));
end

%% 3. MCP-eGFP
CTLTime=DataCTL.ElapsedTime;
CTLFluo=DataCTL.MeanVectorAP;
CTLAPbinID=DataCTL.APbinID;
%plot(NLSTime,NLSSpot)


%Calculate the Accumulation of mRNA
%clear MS2SpotFluo

% Choose a specific AP bin,
%AP = 14;  % AP/0.025 + 1 = position in AP

CTLt13=DataCTL.nc13;
CTLt14=DataCTL.nc14;

clear MS2SpotFluo
% For all AP bins,
for AP=1:41
    clear MS2SpotFluo
    MS2SpotFluo=CTLFluo(:,AP);
    MS2SpotFluo(isnan(MS2SpotFluo))=0;
    
    Sum_MCP_eGFP_nc13(AP) = trapz(CTLTime(CTLt13:CTLt14),MS2SpotFluo(CTLt13:CTLt14));
    Sum_MCP_eGFP_nc14(AP) = trapz(CTLTime(CTLt14:end),MS2SpotFluo(CTLt14:end));
end

%% 4. MCP-NoNLS-mCherry (double dosage)
NoNLSDoubleTime=DataNoNLSDouble2.ElapsedTime;
NoNLSDoubleFluo=DataNoNLSDouble2.MeanVectorAP;
NoNLSDoubleAPbinID=DataNoNLSDouble2.APbinID;
%plot(NLSTime,NLSSpot)


%Calculate the Accumulation of mRNA
%clear MS2SpotFluo

% Choose a specific AP bin,
%AP = 14;  % AP/0.025 + 1 = position in AP

NoNLSDoublet13=DataNoNLSDouble2.nc13;
NoNLSDoublet14=DataNoNLSDouble2.nc14;

% For all AP bins,
for AP=1:41
    clear MS2SpotFluo
    MS2SpotFluo=NoNLSDoubleFluo(:,AP);
    MS2SpotFluo(isnan(MS2SpotFluo))=0;
    
    Sum_NoNLS_mCherry_Double_nc13(AP) = trapz(NoNLSDoubleTime(NoNLSDoublet13:NoNLSDoublet14),MS2SpotFluo(NoNLSDoublet13:NoNLSDoublet14));
    Sum_NoNLS_mCherry_Double_nc14(AP) = trapz(NoNLSDoubleTime(NoNLSDoublet14:end),MS2SpotFluo(NoNLSDoublet14:end));
end

%% 5. MCP-NoNLS-mCherry (1.5 dosage, MCP-NoNLS-mCherry(3)/CyO ; MCP-NoNLS-mCherry(6))
NoNLSTime_15=DataNoNLS_15.ElapsedTime;
NoNLSFluo_15=DataNoNLS_15.MeanVectorAP;
NoNLSAPbinID_15=DataNoNLS_15.APbinID;
%plot(NLSTime,NLSSpot)


%Calculate the Accumulation of mRNA
%clear MS2SpotFluo

% Choose a specific AP bin,
%AP = 14;  % AP/0.025 + 1 = position in AP

NoNLSt13_15=DataNoNLS_15.nc13;
NoNLSt14_15=DataNoNLS_15.nc14;

% For all AP bins,
for AP=1:41
    clear MS2SpotFluo
    MS2SpotFluo=NoNLSFluo_15(:,AP);
    MS2SpotFluo(isnan(MS2SpotFluo))=0;
    
    Sum_NoNLS_mCherry_15_nc13(AP) = trapz(NoNLSTime_15(NoNLSt13_15:NoNLSt14_15),MS2SpotFluo(NoNLSt13_15:NoNLSt14_15));
    Sum_NoNLS_mCherry_15_nc14(AP) = trapz(NoNLSTime_15(NoNLSt14_15:end),MS2SpotFluo(NoNLSt14_15:end));
end

%% plot to check

% Define the constant that calibrates the MCP-mCherry to MCP-eGFP
% A = 1; % Scaling factor for MCP-NoNLS-mCherry
% B = 1; % Scaling factor for MCP-NLS-mCherry
APbins = 0:0.025:1;

figure(1)
hold on
% plot(APbins(8:24),Sum_NoNLS_mCherry_nc13(8:24),'b','LineWidth',5)
% plot(APbins(11:27),Sum_NLS_mCherry_nc13(11:27),'r','LineWidth',5)
% plot(APbins(8:25),Sum_NoNLS_mCherry_Double_nc13(8:25),'g','LineWidth',5)
% plot(APbins(12:28),Sum_NoNLS_mCherry_15_nc13(12:28),'k','LineWidth',5)
% Plot only the APbins that doesn't have Nans.
plot(APbins(8:24),Sum_NoNLS_mCherry_nc13(8:24),'b','LineWidth',5)
plot(APbins(10:26),Sum_NLS_mCherry_nc13(11:27),'r','LineWidth',5)
plot(APbins(8:25),Sum_NoNLS_mCherry_Double_nc13(8:25),'g','LineWidth',5)
plot(APbins(12:28),Sum_NoNLS_mCherry_15_nc13(12:28),'k','LineWidth',5)

title('Comparison of MCPs (Fluorescence integration over nc13)')
xlabel('AP')
ylabel('Fluorescence (AU)')
legend('MCP-NoNLS-mCherry','MCP-NLS-mCherry','MCP-NoNLS-mCherry(double)','MCP-NoNLS-mCherry(x1.5)')

figure(2)
hold on
plot(APbins(8:24),Sum_NoNLS_mCherry_nc14(8:24),'b','LineWidth',5)
plot(APbins(11:27),Sum_NLS_mCherry_nc14(11:27),'r','LineWidth',5)
plot(APbins(8:25),Sum_NoNLS_mCherry_Double_nc14(8:25),'g','LineWidth',5)
plot(APbins(12:28),Sum_NoNLS_mCherry_15_nc14(12:28),'k','LineWidth',5)


title('Comparison of MCPs (Fluorescence integration over nc14)')
xlabel('AP')
ylabel('Fluorescence (AU)')
legend('MCP-NoNLS-mCherry','MCP-NLS-mCherry','MCP-NoNLS-mCherry(double)','MCP-NoNLS-mCherry(x1.5)')



% Think about the error bars -> Use the mean value of SDMeanVector?
% Is there any way that we can extract the information of transcription
% with the under-saturated MCP-mCherry ?

%% Spot Fluorescence comparison (Spot Fluo over time @ specific AP bin)

% Let's start with 0.3 AP bin, which is 13th column of MeanVectorAP
APbin = 12;
NLSAPbin = APbin;
APbin_15 = APbin;

% (mean) Spot fluo of that APbin over time
NLSFluo = DataNLS.MeanVectorAP(:,NLSAPbin);
NoNLSFluo = DataNoNLS.MeanVectorAP(:,APbin);
NoNLSDoubleFluo = DataNoNLSDouble2.MeanVectorAP(:,APbin); 
NoNLS_15Fluo = DataNoNLS_15.MeanVectorAP(:,APbin_15);

% Think about the errorbar 
% I need to change this to SEM
SDNLS = DataNLS.SDVectorAP(:,NLSAPbin);
SDNoNLS = DataNoNLS.SDVectorAP(:,APbin);
SDNoNLSDouble = DataNoNLSDouble2.SDVectorAP(:,APbin); 
SDNoNLS_15 = DataNoNLS_15.SDVectorAP(:,APbin_15);

% Plotting the Spot Fluo over time
figure(3)
hold on
% nc13
errorbar(NLSTime(NLSt13:NLSt14)-NLSTime(NLSt13)-0.6881,NLSFluo(NLSt13:NLSt14),SDNLS(NLSt13:NLSt14),'LineWidth',3)
errorbar(NoNLSTime(NoNLSt13:NoNLSt14)-NoNLSTime(NoNLSt13)-4.109,NoNLSFluo(NoNLSt13:NoNLSt14),SDNoNLS(NoNLSt13:NoNLSt14),'LineWidth',3)
errorbar(NoNLSDoubleTime(NoNLSDoublet13:NoNLSDoublet14)-NoNLSDoubleTime(NoNLSDoublet13)-2.752,NoNLSDoubleFluo(NoNLSDoublet13:NoNLSDoublet14),SDNoNLSDouble(NoNLSDoublet13:NoNLSDoublet14),'LineWidth',3)
errorbar(NoNLSTime_15(NoNLSt13_15:NoNLSt14_15)-NoNLSTime_15(NoNLSt13_15)-2.064,NoNLS_15Fluo(NoNLSt13_15:NoNLSt14_15),SDNoNLS_15(NoNLSt13_15:NoNLSt14_15),'LineWidth',3)
% All the time, just synchronized
% plot((NLSt13-NLSt13:length(NLSFluo)-NLSt13),NLSFluo(NLSt13:end))
% plot((NoNLSt13-NoNLSt13:length(NoNLSFluo)-NoNLSt13),NoNLSFluo(NoNLSt13:end))
% plot((NoNLSDoublet13-NoNLSDoublet13:length(NoNLSDoubleFluo)-NoNLSDoublet13),NoNLSDoubleFluo(NoNLSDoublet13:end))
% plot((NoNLSt13_15-NoNLSt13_15:length(NoNLS_15Fluo)-NoNLSt13_15),NoNLS_15Fluo(NoNLSt13_15:end))

title(['Spot Fluorescence over Time @ AP = ',num2str((APbin-1)*0.025)])
xlabel('Time (min)')
ylabel('Fluorescence (AU)')
legend('MCP-NLS-mCherry','MCP-NoNLS-mCherry','MCP-NoNLS-mCherry(double)','MCP-NoNLS-mCherry(x1.5)')

%% Offset
%NLSOffset = DataNLS.MeanOffsetVector;
% 1x offset
NoNLSOffset = DataNoNLS.MeanOffsetVector;
NoNLSOffset2 = DataNoNLS2.MeanOffsetVector;
% 1x SEM
SENoNLSOffset = DataNoNLS.SDOffsetVector./sqrt(DataNoNLS.NOffsetParticles);
SENoNLSOffset2 = DataNoNLS2.SDOffsetVector./sqrt(DataNoNLS2.NOffsetParticles);

hold on

%plot(1:length(NLSOffset),NLSOffset,'LineWidth',3)
errorbar(15:length(NoNLSOffset),NoNLSOffset(15:end),SENoNLSOffset(15:end),'LineWidth',3)
errorbar(1:length(NoNLSOffset2),NoNLSOffset2,SENoNLSOffset2,'LineWidth',3)
%plot(19:length(NoNLSDoubleOffset),NoNLSDoubleOffset(19:end),'LineWidth',3)
%plot(27:length(NoNLS15Offset),NoNLS15Offset(27:end),'LineWidth',3)

title('Offset')
xlabel('Frames')
ylabel('Fluorescence (AU)')
%legend('1 x NoNLS','2 x NoNLS','1.5 x NoNLS')
%%  1.5x offset
NoNLS15Offset = DataNoNLS_15.MeanOffsetVector;
NoNLS15Offset2 = DataNoNLS_15_2.MeanOffsetVector;
NoNLS15Offset3 = DataNoNLS_15_3.MeanOffsetVector;

% 1.5x SEM
SENoNLS15Offset = DataNoNLS_15.SDOffsetVector./sqrt(DataNoNLS_15.NOffsetParticles);
SENoNLS15Offset2 = DataNoNLS_15_2.SDOffsetVector./sqrt(DataNoNLS_15_2.NOffsetParticles);
SENoNLS15Offset3 = DataNoNLS_15_3.SDOffsetVector./sqrt(DataNoNLS_15_3.NOffsetParticles);

% Plot 1.5x Offset
hold on
errorbar(1:length(NoNLS15Offset),NoNLS15Offset,SENoNLS15Offset,'LineWidth',3)
errorbar(1:length(NoNLS15Offset2),NoNLS15Offset2,SENoNLS15Offset2,'LineWidth',3)
errorbar(1:length(NoNLS15Offset3),NoNLS15Offset3,SENoNLS15Offset3,'LineWidth',3)

%%  2x offset
NoNLSDoubleOffset = DataNoNLSDouble1.MeanOffsetVector;
NoNLSDoubleOffset2 = DataNoNLSDouble2.MeanOffsetVector;

% 2x SEM
SENoNLSDoubleOffset = DataNoNLSDouble1.SDOffsetVector./sqrt(DataNoNLSDouble1.NOffsetParticles);
SENoNLSDoubleOffset2 = DataNoNLSDouble2.SDOffsetVector./sqrt(DataNoNLSDouble2.NOffsetParticles);

hold on
errorbar(1:length(NoNLSDoubleOffset),NoNLSDoubleOffset,SENoNLSDoubleOffset,'LineWidth',3)
errorbar(1:length(NoNLSDoubleOffset2),NoNLSDoubleOffset2,SENoNLSDoubleOffset2,'LineWidth',3)
%% plot
hold on

%plot(1:length(NLSOffset),NLSOffset,'LineWidth',3)
errorbar(15:length(NoNLSOffset),NoNLSOffset(15:end),SENoNLSOffset(15:end),'LineWidth',3)
errorbar(1:length(NoNLSOffset2),NoNLSOffset2,SENoNLSOffset2,'LineWidth',3)
%plot(19:length(NoNLSDoubleOffset),NoNLSDoubleOffset(19:end),'LineWidth',3)
%plot(27:length(NoNLS15Offset),NoNLS15Offset(27:end),'LineWidth',3)

title('Offset')
xlabel('Frames')
ylabel('Fluorescence (AU)')
%legend('1 x NoNLS','2 x NoNLS','1.5 x NoNLS')

%% Standard error of mean (SEM) calculation (10/9/2017. YJK) - x1.5 dosage
% This code calculates the SEM of MeanVectorAP (Spot fluo) for multiple
% datasets at a given (specific) AP bin.

AP = 12;
% First, we have to synchronize the traces. For now, only think about nc 13
Dosage_15(1).Time = (DataNoNLS_15.nc13:DataNoNLS_15.nc14) - DataNoNLS_15.nc13;
Dosage_15(2).Time = (DataNoNLS_15_2.nc13:DataNoNLS_15_2.nc14) - DataNoNLS_15_2.nc13-3;
Dosage_15(3).Time = (DataNoNLS_15_3.nc13+1:DataNoNLS_15_3.nc14) - DataNoNLS_15_3.nc13;

% Second, MeanVectorAP (Spot fluorescence)
Dosage_15(1).SpotFluo = DataNoNLS_15.MeanVectorAP(DataNoNLS_15.nc13:DataNoNLS_15.nc14,:);
Dosage_15(2).SpotFluo = DataNoNLS_15_2.MeanVectorAP(DataNoNLS_15_2.nc13:DataNoNLS_15_2.nc14,:);
Dosage_15(3).SpotFluo = DataNoNLS_15_3.MeanVectorAP(DataNoNLS_15_3.nc13+1:DataNoNLS_15_3.nc14,:);

% Plot to check the synchronization
hold on
plot(Dosage_15(1).Time,Dosage_15(1).SpotFluo(:,AP))
plot(Dosage_15(2).Time,Dosage_15(2).SpotFluo(:,AP))
plot(Dosage_15(3).Time,Dosage_15(3).SpotFluo(:,AP))
legend('1','2','3')
%%
% Third, Number of Particles in AP bins, at each time point
Dosage_15(1).NParticles = DataNoNLS_15.NParticlesAP(DataNoNLS_15.nc13:DataNoNLS_15.nc14,:);
Dosage_15(2).NParticles = DataNoNLS_15_2.NParticlesAP(DataNoNLS_15_2.nc13:DataNoNLS_15_2.nc14,:)
Dosage_15(3).NParticles = DataNoNLS_15_3.NParticlesAP(DataNoNLS_15_3.nc13+1:DataNoNLS_15_3.nc14,:)

% Fourth, Standard deviation of each dataset
Dosage_15(1).SD = DataNoNLS_15.SDVectorAP(DataNoNLS_15.nc13:DataNoNLS_15.nc14,:);
Dosage_15(2).SD = DataNoNLS_15_2.SDVectorAP(DataNoNLS_15_2.nc13:DataNoNLS_15_2.nc14,:)
Dosage_15(3).SD = DataNoNLS_15_3.SDVectorAP(DataNoNLS_15_3.nc13+1:DataNoNLS_15_3.nc14,:)

%% Calculate the SEM for multiple datasets
VarSum = zeros(size(Dosage_15(1).NParticles));
Sum = zeros(size(Dosage_15(1).NParticles));
TotalNumber = zeros(size(Dosage_15(1).NParticles));
Denominator = zeros(size(Dosage_15(1).NParticles));

for i=1:length(Dosage_15)
    % Remove Nans
    Dosage_15(i).SpotFluo(isnan(Dosage_15(i).SpotFluo)) = 0;
    Dosage_15(i).NParticles(isnan(Dosage_15(i).NParticles)) = 0;
    Dosage_15(i).SD(isnan(Dosage_15(i).SD)) = 0;
    
    Sum = Sum + (Dosage_15(i).SpotFluo.*Dosage_15(i).NParticles);
    
    VarSum = VarSum + (Dosage_15(i).SD).^2 .* (Dosage_15(i).NParticles-1);
    TotalNumber = TotalNumber + Dosage_15(i).NParticles;
    Denominator = Denominator + Dosage_15(i).NParticles-1;
end

Mean = Sum./TotalNumber;
SEM = sqrt(VarSum ./ Denominator);

%%
% Now, plot the mean of spot fluo (Mean), with the errorbar of SEM
% First, pick an AP bin,
AP = 14;

hold on
plot(Dosage_15(1).Time,Dosage_15(1).SpotFluo(:,AP))
plot(Dosage_15(2).Time,Dosage_15(2).SpotFluo(:,AP))
plot(Dosage_15(3).Time,Dosage_15(3).SpotFluo(:,AP))
errorbar(Dosage_15(1).Time , Mean(:,AP) , SEM(:,AP),'LineWidth',3)
legend('1','2','3','Mean')

%% Double dosage data - Mean and SEM
% Standard error of mean (SEM) calculation
% This code calculates the SEM of MeanVectorAP (Spot fluo) for multiple
% datasets at a given (specific) AP bin.

AP = 14;
% First, we have to synchronize the traces. For now, only think about nc 13
Dosage_double(1).Time = (DataNoNLSDouble1.nc13:DataNoNLSDouble1.nc14) - DataNoNLSDouble1.nc13;
Dosage_double(2).Time = (DataNoNLSDouble2.nc13+1:DataNoNLSDouble2.nc14) - DataNoNLSDouble2.nc13-2;

% Second, MeanVectorAP (Spot fluorescence)
Dosage_double(1).SpotFluo = DataNoNLSDouble1.MeanVectorAP(DataNoNLSDouble1.nc13:DataNoNLSDouble1.nc14,:);
Dosage_double(2).SpotFluo = DataNoNLSDouble2.MeanVectorAP(DataNoNLSDouble2.nc13+1:DataNoNLSDouble2.nc14,:);

% Plot to check the synchronization
% hold on
% plot(Dosage_double(1).Time,Dosage_double(1).SpotFluo(:,AP))
% plot(Dosage_double(2).Time,Dosage_double(2).SpotFluo(:,AP))
% 
% legend('1','2')

% Third, Number of Particles in AP bins, at each time point
Dosage_double(1).NParticles = DataNoNLSDouble1.NParticlesAP(DataNoNLSDouble1.nc13:DataNoNLSDouble1.nc14,:);
Dosage_double(2).NParticles = DataNoNLSDouble2.NParticlesAP(DataNoNLSDouble2.nc13+1:DataNoNLSDouble2.nc14,:)

% Fourth, Standard deviation of each dataset
Dosage_double(1).SD = DataNoNLSDouble1.SDVectorAP(DataNoNLSDouble1.nc13:DataNoNLSDouble1.nc14,:);
Dosage_double(2).SD = DataNoNLSDouble2.SDVectorAP(DataNoNLSDouble2.nc13+1:DataNoNLSDouble2.nc14,:)

%% Calculate the SEM for x2 dosage dataset
VarSum = zeros(size(Dosage_double(1).NParticles));
Sum = zeros(size(Dosage_double(1).NParticles));
TotalNumber = zeros(size(Dosage_double(1).NParticles));
Denominator = zeros(size(Dosage_double(1).NParticles));

for i=1:length(Dosage_double)
    % Remove Nans
    Dosage_double(i).SpotFluo(isnan(Dosage_double(i).SpotFluo)) = 0;
    Dosage_double(i).NParticles(isnan(Dosage_double(i).NParticles)) = 0;
    Dosage_double(i).SD(isnan(Dosage_double(i).SD)) = 0;
    
    Sum = Sum + (Dosage_double(i).SpotFluo.*Dosage_double(i).NParticles);
    
    VarSum = VarSum + (Dosage_double(i).SD).^2 .* (Dosage_double(i).NParticles-1);
    TotalNumber = TotalNumber + Dosage_double(i).NParticles;
    Denominator = Denominator + Dosage_double(i).NParticles-1;
end

Mean_double = Sum./TotalNumber;
SEM_double = sqrt(VarSum ./ Denominator);

%%
% Now, plot the mean of spot fluo (Mean), with the errorbar of SEM
% First, pick an AP bin,
AP = 14;

hold on
plot(Dosage_double(1).Time,Dosage_double(1).SpotFluo(:,AP))
plot(Dosage_double(2).Time,Dosage_double(2).SpotFluo(:,AP))

errorbar(Dosage_double(1).Time , Mean_double(:,AP) , SEM_double(:,AP),'LineWidth',3)
legend('1','2','Mean')

%% 1x dosage data - Mean and SEM
% Standard error of mean (SEM) calculation
% This code calculates the SEM of MeanVectorAP (Spot fluo) for multiple
% datasets at a given (specific) AP bin.

AP = 14;
% First, we have to synchronize the traces. For now, only think about nc 13
Dosagex1(1).Time = (DataNoNLS.nc13+4:DataNoNLS.nc14) - DataNoNLS.nc13-3;
Dosagex1(2).Time = (DataNoNLS2.nc13:DataNoNLS2.nc14) - DataNoNLS2.nc13-4;

% Second, MeanVectorAP (Spot fluorescence)
Dosagex1(1).SpotFluo = DataNoNLS.MeanVectorAP(DataNoNLS.nc13+4:DataNoNLS.nc14,:);
Dosagex1(2).SpotFluo = DataNoNLS2.MeanVectorAP(DataNoNLS2.nc13:DataNoNLS2.nc14,:);

% Plot to check the synchronization
hold on
plot(Dosagex1(1).Time,Dosagex1(1).SpotFluo(:,AP))
plot(Dosagex1(2).Time,Dosagex1(2).SpotFluo(:,AP))
% 
legend('1','2')
%% 
% Third, Number of Particles in AP bins, at each time point
Dosagex1(1).NParticles = DataNoNLS.NParticlesAP(DataNoNLS.nc13+4:DataNoNLS.nc14,:);
Dosagex1(2).NParticles = DataNoNLS2.NParticlesAP(DataNoNLS2.nc13:DataNoNLS2.nc14,:)

% Fourth, Standard deviation of each dataset
Dosagex1(1).SD = DataNoNLS.SDVectorAP(DataNoNLS.nc13+4:DataNoNLS.nc14,:);
Dosagex1(2).SD = DataNoNLS2.SDVectorAP(DataNoNLS2.nc13:DataNoNLS2.nc14,:)

%% Calculate the SEM for multiple datasets
VarSum = zeros(size(Dosagex1(1).NParticles));
Sum = zeros(size(Dosagex1(1).NParticles));
TotalNumber = zeros(size(Dosagex1(1).NParticles));
Denominator = zeros(size(Dosagex1(1).NParticles));

for i=1:length(Dosagex1)
    % Remove Nans
    Dosagex1(i).SpotFluo(isnan(Dosagex1(i).SpotFluo)) = 0;
    Dosagex1(i).NParticles(isnan(Dosagex1(i).NParticles)) = 0;
    Dosagex1(i).SD(isnan(Dosagex1(i).SD)) = 0;
    
    Sum = Sum + (Dosagex1(i).SpotFluo.*Dosagex1(i).NParticles);
    
    VarSum = VarSum + (Dosagex1(i).SD).^2 .* (Dosagex1(i).NParticles-1);
    TotalNumber = TotalNumber + Dosagex1(i).NParticles;
    Denominator = Denominator + Dosagex1(i).NParticles-1;
end

Mean_x1 = Sum./TotalNumber;
SEM_x1 = sqrt(VarSum ./ Denominator);

%%
% Now, plot the mean of spot fluo (Mean), with the errorbar of SEM
% First, pick an AP bin,
AP = 14;

hold on
plot(Dosagex1(1).Time,Dosagex1(1).SpotFluo(:,AP))
plot(Dosagex1(2).Time,Dosagex1(2).SpotFluo(:,AP))

errorbar(Dosagex1(1).Time , Mean_x1(:,AP) , SEM_x1(:,AP),'LineWidth',3)
legend('1','2','Mean')



%% Plot different dosage data using Mean and SEM

AP = 8;
hold on
errorbar(Dosagex1(1).Time*2/3-2,Mean_x1(:,AP),SEM_x1(:,AP),'LineWidth',3)
errorbar(Dosage_15(1).Time*2/3-1 , Mean(:,AP) , SEM(:,AP),'LineWidth',3)
errorbar(Dosage_double(1).Time*2/3 -2 , Mean_double(:,AP) , SEM_double(:,AP),'LineWidth',3)

legend('x1','x1.5','x2')
title(['Mean Spot Fluorescence','@','APbin=',num2str(AP)])
xlabel('Time (min)')
ylabel('Mean Spot Fluorescence (AU)')

%% Compare the Maximum fluorescence at nc 14

% First, 1x dosage
AP = 14;
NoNLSnc14 = DataNoNLS.nc14;
NoNLS2nc14 = DataNoNLS2.nc14;
%plot(DataNoNLS.ElapsedTime(NoNLSnc14:end),DataNoNLS.MeanVectorAP((NoNLSnc14:end),AP))
Max1x = max(DataNoNLS.MeanVectorAP((NoNLSnc14:end),AP));
Max1x_2 = max(DataNoNLS2.MeanVectorAP((NoNLS2nc14:end),AP));

% 1x SEM
Max1x_SEM = DataNoNLS.SDVectorAP(60,AP)/sqrt(DataNoNLS.NParticlesAP(60,AP));
Max1x_2_SEM = DataNoNLS2.SDVectorAP(70,AP)/sqrt(DataNoNLS2.NParticlesAP(70,AP));

% 1x Mean (of all embryos)
MaxSumTemp = (Max1x*(DataNoNLS.NParticlesAP(60,AP)) + Max1x_2*(DataNoNLS2.NParticlesAP(70,AP)));
TotalNumberTemp = (DataNoNLS.NParticlesAP(60,AP))+(DataNoNLS2.NParticlesAP(70,AP));

MeanMax1x = MaxSumTemp / TotalNumberTemp;
% 1x SEM (of all embryos)
% SEM = sqrt (SD1^2*(n1-1)+SD2^2*(n2-1).../(n1+n2+....-N))
SEMSumTemp = (DataNoNLS.SDVectorAP(60,AP).^2)*(DataNoNLS.NParticlesAP(60,AP)-1)+...
                (DataNoNLS2.SDVectorAP(70,AP).^2)*(DataNoNLS2.NParticlesAP(70,AP)-1);
SEM1x = sqrt(SEMSumTemp / (TotalNumberTemp-2));

% 1.5 x dosage
NoNLS_15_nc14 = DataNoNLS_15.nc14;
NoNLS_15_2_nc14 = DataNoNLS_15_2.nc14;
NoNLS_15_3_nc14 = DataNoNLS_15_3.nc14;


%plot(DataNoNLS_15.ElapsedTime,DataNoNLS_15.MeanVectorAP(:,AP))
Max15x = max(DataNoNLS_15.MeanVectorAP((NoNLS_15_nc14:end),AP));
Max15x_2 = max(DataNoNLS_15_2.MeanVectorAP((NoNLS_15_2_nc14:end),AP));
Max15x_3 = max(DataNoNLS_15_3.MeanVectorAP((NoNLS_15_3_nc14:end),AP));

% SEM
Max15x_SEM = DataNoNLS_15.SDVectorAP(56,AP)/sqrt(DataNoNLS_15.NParticlesAP(56,AP));
Max15x_2_SEM = DataNoNLS_15_2.SDVectorAP(78,AP)/sqrt(DataNoNLS_15_2.NParticlesAP(78,AP));
Max15x_3_SEM = DataNoNLS_15_3.SDVectorAP(67,AP)/sqrt(DataNoNLS_15_3.NParticlesAP(67,AP));

% 2x dosage
NoNLS_Double1_nc14 = DataNoNLSDouble1.nc14;
NoNLS_Double2_nc14 = DataNoNLSDouble2.nc14;

MaxDouble1 = max(DataNoNLSDouble1.MeanVectorAP((NoNLS_Double1_nc14:end),AP));
MaxDouble2 = max(DataNoNLSDouble2.MeanVectorAP((NoNLS_Double2_nc14:end),AP));
% SEM 2x
MaxDouble1_SEM = DataNoNLSDouble1.SDVectorAP(47,AP)/sqrt(DataNoNLSDouble1.NParticlesAP(47,AP));
MaxDouble2_SEM = DataNoNLSDouble2.SDVectorAP(47,AP)/sqrt(DataNoNLSDouble1.NParticlesAP(47,AP));

Dosage = [5.5 5 10.5 10.5 9 14 13]/5;
%Dosage = [1 1 2 2 2 2.5 2.5]; % Dosage matrix , I need to change for different datasets
MaxFluo = [Max1x Max1x_2 Max15x Max15x_2 Max15x_3 MaxDouble1 MaxDouble2];
% Get the SEM at the maximum value
ErrorFluo = [Max1x_SEM Max1x_2_SEM Max15x_SEM Max15x_2_SEM Max15x_3_SEM MaxDouble1_SEM MaxDouble2_SEM];

errorbar(Dosage, MaxFluo,ErrorFluo,'o')
ylim([0 1000])
xlim([0.5 3])
xlabel('MCP-mCherry Dosage')
ylabel('Maximum Fluorescence at nc 14 (AU)')
title(['Maximum Fluorescence for different MCP dosages','@ AP bin =',num2str(AP)])
set(gca,'fontsize',20)
% Maybe I have some errors in FindAPAxisFullEmbryio
% I need to go back to the movies, and double check!
% Also, the variable name sucks, I need to think about a better way
