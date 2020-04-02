% [MCP-mCherry dosage] Compare the nc 14 fluorescence traces for different dosages

% Select one AP bin
AP = 14;

% Define the spot fluorescence and error bar (1x)
NoNLSTime = DataNoNLS.ElapsedTime;
NoNLSTime_2 = DataNoNLS2.ElapsedTime;

NoNLSnc14 = DataNoNLS.nc14;
NoNLS2nc14 = DataNoNLS2.nc14;

SpotFluo1x = DataNoNLS.MeanVectorAP(:,AP);
SpotFluo1x_2 = DataNoNLS2.MeanVectorAP(:,AP);

ErrorFluo1x = DataNoNLS.MeanVectorAP(:,AP) / sqrt(DataNoNLS.NParticlesAP(:,AP));
ErrorFluo1x_2 = DataNoNLS2.MeanVectorAP(:,AP) / sqrt(DataNoNLS2.NParticlesAP(:,AP));

% Define the spot fluorescence and error bar (2x(previously notated as 1.5x)
% genotype : MCPmcherry(3)/CyO;MCPmCherry(6)
NoNLSTime15 = DataNoNLS_15.ElapsedTime;
NoNLSTime15_2 = DataNoNLS_15_2.ElapsedTime;
NoNLSTime15_3 = DataNoNLS_15_3.ElapsedTime;

NoNLS_15_nc14 = DataNoNLS_15.nc14;
NoNLS_15_2_nc14 = DataNoNLS_15_2.nc14;
NoNLS_15_3_nc14 = DataNoNLS_15_3.nc14;

SpotFluo15x = DataNoNLS_15.MeanVectorAP(:,AP);
SpotFluo15x_2 = DataNoNLS_15_2.MeanVectorAP(:,AP);
SpotFluo15x_3 = DataNoNLS_15_3.MeanVectorAP(:,AP);

ErrorFluo15x = DataNoNLS_15.MeanVectorAP(:,AP) / sqrt(DataNoNLS_15.NParticlesAP(:,AP));
ErrorFluo15x_2 = DataNoNLS_15_2.MeanVectorAP(:,AP) / sqrt(DataNoNLS_15_2.NParticlesAP(:,AP));
ErrorFluo15x_3 = DataNoNLS_15_3.MeanVectorAP(:,AP) / sqrt(DataNoNLS_15_3.NParticlesAP(:,AP));

% Define the spot fluorescence and error bar (2.5x, previously notated as 2x)
% genotype : MCPmCherry(3);MCPmCherry(6)

NoNLSTimeDouble = DataNoNLSDouble1.ElapsedTime;
NoNLSTimeDouble2 = DataNoNLSDouble2.ElapsedTime;

NoNLSDoublenc14 = DataNoNLSDouble1.nc14;
NoNLSDouble2nc14 = DataNoNLSDouble2.nc14;

SpotFluoDouble1 = DataNoNLSDouble1.MeanVectorAP(:,AP);
SpotFluoDouble2 = DataNoNLSDouble2.MeanVectorAP(:,AP);

ErrorFluoDouble1 = DataNoNLSDouble1.MeanVectorAP(:,AP) / sqrt(DataNoNLSDouble1.NParticlesAP(:,AP));
ErrorFluoDouble2 = DataNoNLSDouble2.MeanVectorAP(:,AP) / sqrt(DataNoNLSDouble2.NParticlesAP(:,AP));


% plot all spot fluo traces with error bar (SEM)
hold on
% 1x
errorbar(NoNLSTime(NoNLSnc14:end)-NoNLSTime(NoNLSnc14)-3.5,...
        SpotFluo1x(NoNLSnc14:end),ErrorFluo1x((NoNLSnc14:end),AP),'b')
    
errorbar(NoNLSTime_2(NoNLS2nc14:end)-NoNLSTime_2(NoNLS2nc14)-4,...
        SpotFluo1x_2(NoNLS2nc14:end),ErrorFluo1x_2((NoNLS2nc14:end),AP),'b')
    
% 2x (previously notated as 1.5x)
errorbar(NoNLSTime15(NoNLS_15_nc14:end)-NoNLSTime15(NoNLS_15_nc14),...
        SpotFluo15x(NoNLS_15_nc14:end),ErrorFluo15x((NoNLS_15_nc14:end),AP),'g')

errorbar(NoNLSTime15_2(NoNLS_15_2_nc14:end)-NoNLSTime15_2(NoNLS_15_2_nc14),...
        SpotFluo15x_2(NoNLS_15_2_nc14:end),ErrorFluo15x_2((NoNLS_15_2_nc14:end),AP),'g')
    
errorbar(NoNLSTime15_3(NoNLS_15_3_nc14:end)-NoNLSTime15_3(NoNLS_15_3_nc14),...
        SpotFluo15x_3(NoNLS_15_3_nc14:end),ErrorFluo15x_3((NoNLS_15_3_nc14:end),AP),'g')

% 2.5x (previously notated as 2x)
errorbar(NoNLSTimeDouble(NoNLSDoublenc14:end)-NoNLSTimeDouble(NoNLSDoublenc14)-1.5,...
        SpotFluoDouble1(NoNLSDoublenc14:end),ErrorFluoDouble1((NoNLSDoublenc14:end),AP),'r')
    
errorbar(NoNLSTimeDouble2(NoNLSDouble2nc14:end)-NoNLSTimeDouble2(NoNLSDouble2nc14)-4,...
        SpotFluoDouble2(NoNLSDouble2nc14:end),ErrorFluoDouble2((NoNLSDouble2nc14:end),AP),'r')
    
legend('1x','1x','2x','2x','2x','2.5x','2.5x')
title(['Mean Spot Fluorescence at nc 14',' @ AP bin',num2str(AP)])
xlabel('Time(min)')
ylabel('Fluorescence (AU)')
set(gca,'fontsize',20)
set(gca,'LineWidth',3)