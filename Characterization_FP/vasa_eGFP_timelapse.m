% vasa-eGFP level (fluorescence intensity) over time
%clear all
% Load the Data
vasaeGFP = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-10-11-vasa-eGFP-His-iRFP\CompiledNuclei.mat')
%vasaeGFP2 = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-09-24-vasa-eGFP-His-iRFP2\CompiledNuclei.mat')
vasaeGFP2 = load('D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-09-24-vasa-eGFP-His-iRFP\CompiledNuclei.mat')


Time = vasaeGFP.ElapsedTime;
AvgFluo = vasaeGFP.MeanVectorAll;
CytoFluo = vasaeGFP.MeanCytoAPProfile{1};
CytoFluo = nanmean(CytoFluo,1);
SDFluo = vasaeGFP.SDVectorAll;
NParticles = vasaeGFP.NParticlesAll;
SEFluo = SDFluo./sqrt(NParticles);
emb1nc12 = vasaeGFP.nc12;
emb1nc13 = vasaeGFP.nc13;
emb1nc14 = vasaeGFP.nc14;

Time2 = vasaeGFP2.ElapsedTime;
AvgFluo2 = vasaeGFP2.MeanVectorAll;
CytoFluo2 = vasaeGFP2.MedianCyto;
SDFluo2 = vasaeGFP2.SDVectorAll;
NParticles2 = vasaeGFP2.NParticlesAll;
SEFluo2 = SDFluo2./sqrt(NParticles2);
emb2nc12 = vasaeGFP2.nc12;
emb2nc13 = vasaeGFP2.nc13;
emb2nc14 = vasaeGFP2.nc14;
%% Plot
hold on
% plot(Time,AvgFluo)

errorbar(Time(emb1nc13:end)-Time(emb1nc13),AvgFluo(emb1nc13:end),SEFluo(emb1nc13:end))
errorbar(Time2(emb2nc13:end)-Time2(emb2nc13),AvgFluo2(emb2nc13:end),SEFluo2(emb2nc13:end))
%plot(Time(emb1nc13:end)-Time(emb1nc13),CytoFluo(emb1nc13:end))

% plot(Time2,AvgFluo2)
%plot(Time2,CytoFluo2)
%ylim([0 100])
title('vasa-eGFP level')
xlabel('Time (min)')
ylabel('Nuclear fluorescence (AU)')
legend('embryo1,15uW,1min','embryo2,20uW','embryo1,cyto fluo')
set(gca,'Fontsize',20)
