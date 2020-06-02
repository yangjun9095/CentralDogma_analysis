function [thing] = correlatemRNAProtein(Prefix,gene,apBinOfInterest,spotFluorophore,startFrame,endFrame,startAPBin,endAPBin)


%% Load relavent data and define relavent variables
particleFile = strcat('E:\EvanM\LivemRNA\Data\DynamicsResults\',Prefix,'\CompiledParticles.mat');
nucleiFile = strcat('E:\EvanM\LivemRNA\Data\DynamicsResults\',Prefix,'\CompiledNuclei.mat');

compParticles = load(particleFile);
compNuclei = load(nucleiFile);
cutParticleMeanVectorAP = compParticles.MeanVectorAP(startFrame:endFrame,startAPBin:endAPBin);
cutParticleMeanVectorAP(isnan(cutParticleMeanVectorAP)) = 0; %MAYBE THIS SHOULD BE THE OFFSET AND NOT 0
proteinFluo = compNuclei.MeanVectorAP(startFrame:endFrame,startAPBin:endAPBin);
proteinFluo(isnan(proteinFluo)) = 0;
OGNumAPBins = size(compParticles.MeanVectorAP,2);
numAPBins = size(cutParticleMeanVectorAP,2);
numTimeBins = size(cutParticleMeanVectorAP,1);
adjustedElapsedTime = zeros(1,numTimeBins);


% Redefine zero time point.
for i=1:numTimeBins
	adjustedElapsedTime(i) = compParticles.ElapsedTime(1,startFrame+i-1)-...
        compParticles.ElapsedTime(1,startFrame);
end

adjustedElapsedTime = adjustedElapsedTime(:);

%% Assign values to diffusion, etc. constants

 %All constants in units of seconds and microns.

 %Dm: Diffusion constant of mRNA
 %Dp: Diffusion constant of protein
 %Tm: Half-life of mRNA
 %Tp: Half-life of protein
 %rp: Rate of protein synthesis (molecules/mRNA/min)
 
if strcmpi('Hb',gene) || strcmpi('Hunchback',gene)
    Dm = 0.06;
    Dp = 420.0;
    Tm = 60.0;
    Tp = 50.0;
    rp = 2.0;
    
else 
    error(['Gene not found. Either check your spelling or add the gene to this function''s code.' ...
        'Accepted inputs are of the form "Hb" or "Hunchback." Case does not matter.'])
end


%% Calibration of mRNA from fluorescence
 %THIS SHOULD BE CHANGED TO ACCURATELY SCALE BASED ON HERNAN'S 2013 PAPER
 %Approximate max fluorescence is 1200 and the corresponding number of mRNA 
 %loaded is 100. So an approximate scaling is 12.
 
 %FOR NOW THIS IS BEING USED TO SYNC PROTEIN FLUORESCENCE WITH RNA FLUO

RNAProductionRate = cutParticleMeanVectorAP/12;

for i = 23:numAPBins
    for j = 1:numTimeBins
        RNAProductionRate(j,i) = 0;
    end
end

%% Plot RNA production rate

color = jet(numTimeBins);

figure(1)
hold on


for i=1:numTimeBins
	plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,RNAProductionRate(i,:),'color',color(i,:))
end

title(['raw RNA Production Rate: Late nc 12 to nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('Protein Fluorescence (AU)')


%% Plot average spot fluorescence trace for apBinOfInterest

% hold on
% plot(adjustedElapsedTime,cutParticleMeanVectorAP(:,apBinOfInterest))
% hold off
% 
% title('MS2 spot fluorescence for apBinOfInterest')
% xlabel('Time(min)')
% ylabel('Spot Fluorescence (AU)')
% legend(spotFluorophore)


%% Plot offset




%% Temporal Interpolation of spot fluorescence

%Similar to the spatial interpolation, except time points are given by 
%variable 'adjustedElapsedTime'. scaledT divides the intervals of adjustedElapsedTime by 10.

% totalTime = max(adjustedElapsedTime);
% dt = 0.1;
% scaledT = totalTime/numTimeBins/10/2:totalTime/numTimeBins/10:totalTime-totalTime/numTimeBins/10/2;

% totalTime = max(adjustedElapsedTime);
% dt = 0.5;
% scaledT = adjustedElapsedTime(1):dt*adjustedElapsedTime(2):totalTime;
% scaledT = scaledT(:);

dt = 0.25;
x = 1:1:length(adjustedElapsedTime)
xx = 1:dt:length(adjustedElapsedTime)
scaledT = pchip(x,adjustedElapsedTime,xx);


tempInterRNAProductionRate = zeros(length(scaledT),numAPBins);

for j=1:numAPBins
    tempInterRNAProductionRate(:,j) = pchip(adjustedElapsedTime,RNAProductionRate(:,j),scaledT);
end

tempExamp = pchip(adjustedElapsedTime,RNAProductionRate(:,4),scaledT);


%% Plot Temporal Interpolation of spot

color = jet(length(scaledT));

figure(77)
hold on


for i=1:length(scaledT)
	plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,tempInterRNAProductionRate(i,:),'color',color(i,:))
end

title(['temp RNA Production Rate: Late nc 12 to nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('Protein Fluorescence (AU)')




%% Accumulated mRNA from raw data

%Find out how to get embryo length.

RNAProductionRate(isnan(RNAProductionRate))=0;
rawAccumulatedRNA = zeros(numTimeBins,numAPBins);
dt = adjustedElapsedTime(2);
dx = 500/numAPBins;

rawAccumulatedRNA(1,:) = dt*RNAProductionRate(1,:);

for i=2:numTimeBins
	for j=1:numAPBins
        
        %First AP bin. No diffusion to the left.
		if j==1
			rawAccumulatedRNA(i,j) = rawAccumulatedRNA(i-1,j)+... %previously accumulated mRNA
            dt*RNAProductionRate(i,j)-... %synthesized mRNA
            dt*log(2)/Tm*rawAccumulatedRNA(i-1,j)+... %degraded mRNA
            dt*Dm/(dx^2)*(rawAccumulatedRNA(i-1,j+1)-rawAccumulatedRNA(i-1,j)); %diffused mRNA

        %Last AP bin. No diffusion to the right.
        elseif j==numAPBins
            rawAccumulatedRNA(i,j) = rawAccumulatedRNA(i-1,j)+... %previously accumulated mRNA
            dt*RNAProductionRate(i,j)-... %synthesized mRNA
            dt*log(2)/Tm*rawAccumulatedRNA(i-1,j)+... %degraded mRNA
            dt*Dm/(dx^2)*(rawAccumulatedRNA(i-1,j-1)-rawAccumulatedRNA(i-1,j)); %diffused mRNA

        %All other AP bins.
        else
            rawAccumulatedRNA(i,j) = rawAccumulatedRNA(i-1,j)+... %previously accumulated mRNA
            dt*RNAProductionRate(i,j)-... %synthesized mRNA
            dt*log(2)/Tm*rawAccumulatedRNA(i-1,j)+... %degraded mRNA
            dt*Dm/(dx^2)*(rawAccumulatedRNA(i-1,j+1)+rawAccumulatedRNA(i-1,j-1)-... %diffused mRNA
            2*rawAccumulatedRNA(i-1,j));

        end
    end
end


%% Plot accumulated mRNA from raw data

color = jet(numTimeBins);

figure(70)
hold on


for i=1:numTimeBins
	plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,rawAccumulatedRNA(i,:),'color',color(i,:))
end

title(['Accumulated mRNA: nc 13 to ~40 min into nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('mRNA (AU)')



%% Accumulated mRNA from temporal interpolation

tempInterRNAProductionRate(isnan(tempInterRNAProductionRate))=0;
tempInterAccumulatedRNA = zeros(length(scaledT),numAPBins);
dt = scaledT(2);
dx = 500/numAPBins;

tempInterAccumulatedRNA(1,:) = dt*tempInterRNAProductionRate(1,:);

for i=2:length(scaledT)
	for j=1:numAPBins
        
        %First AP bin. No diffusion to the left.
		if j==1
			tempInterAccumulatedRNA(i,j) = tempInterAccumulatedRNA(i-1,j)+... %previously accumulated mRNA
            dt*tempInterRNAProductionRate(i,j)-... %synthesized mRNA
            dt*log(2)/Tm*tempInterAccumulatedRNA(i-1,j)+... %degraded mRNA
            dt*Dm/(dx^2)*(tempInterAccumulatedRNA(i-1,j+1)-tempInterAccumulatedRNA(i-1,j)); %diffused mRNA
        
        %Last AP bin. No diffusion to the right.
        elseif j==numAPBins
            tempInterAccumulatedRNA(i,j) = tempInterAccumulatedRNA(i-1,j)+... %previously accumulated mRNA
            dt*tempInterRNAProductionRate(i,j)-... %synthesized mRNA
            dt*log(2)/Tm*tempInterAccumulatedRNA(i-1,j)+... %degraded mRNA
            dt*Dm/(dx^2)*(tempInterAccumulatedRNA(i-1,j-1)-tempInterAccumulatedRNA(i-1,j)); %diffused mRNA

        %All other AP bins.
        else
            tempInterAccumulatedRNA(i,j) = tempInterAccumulatedRNA(i-1,j)+... %previously accumulated mRNA
            dt*tempInterRNAProductionRate(i,j)-... %synthesized mRNA
            dt*log(2)/Tm*tempInterAccumulatedRNA(i-1,j)+... %degraded mRNA
            dt*Dm/(dx^2)*(tempInterAccumulatedRNA(i-1,j+1)+tempInterAccumulatedRNA(i-1,j-1)-... %diffused mRNA
            2*tempInterAccumulatedRNA(i-1,j));

        end
    end
end

%% Plot accumulated mRNA from temp interpolation

color = jet(length(scaledT));

figure(50)
hold on


for i=1:length(scaledT)
	plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,tempInterAccumulatedRNA(i,:),'color',color(i,:))
end

title(['Accumulated mRNA: Late nc 12 to nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('Protein Fluorescence (AU)')

%% Plot accumulated mRNA from temp interpolation at certain time point


figure(60)
hold on

plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,tempInterAccumulatedRNA(60,:))


title(['Accumulated mRNA: Late nc 12 to nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('Protein Fluorescence (AU)')


%% Calculation of temp inter Accumulated mRNA for different Dm values

%RNAProductionRate(isnan(RNAProductionRate))=0;
DmVar = logspace(-1,2.5,10);
%DmVar = logspace(1,2.4,10);
tempInterAccumulatedRNAVar = zeros(length(scaledT),numAPBins,length(DmVar));
dt = scaledT(2);
dx = 500/numAPBins;

for k=1:length(DmVar)
    tempInterAccumulatedRNAVar(1,:,k) = dt*tempInterRNAProductionRate(1,:);

    for i=2:length(scaledT)
        for j=1:numAPBins
        
            %First AP bin. No diffusion to the left.
            if j==1
                tempInterAccumulatedRNAVar(i,j,k) = tempInterAccumulatedRNAVar(i-1,j,k)+... %previously accumulated RNA
                dt*tempInterRNAProductionRate(i,j)-... %synthesized RNA
                dt*log(2)/Tm*tempInterAccumulatedRNAVar(i-1,j,k)+... %degraded RNA
                dt*DmVar(k)/(dx^2)*(tempInterAccumulatedRNAVar(i-1,j+1,k)-tempInterAccumulatedRNAVar(i-1,j,k)); %diffused RNA

            %Last AP bin. No diffusion to the right.
            elseif j==numAPBins
                tempInterAccumulatedRNAVar(i,j,k) = tempInterAccumulatedRNAVar(i-1,j,k)+... %previously accumulated RNA
                dt*tempInterRNAProductionRate(i,j)-... %synthesized RNA
                dt*log(2)/Tm*tempInterAccumulatedRNAVar(i-1,j,k)+... %degraded RNA
                dt*DmVar(k)/(dx^2)*(tempInterAccumulatedRNAVar(i-1,j-1,k)-tempInterAccumulatedRNAVar(i-1,j,k)); %diffused RNA

            %All other AP bins.
            else
                tempInterAccumulatedRNAVar(i,j,k) = tempInterAccumulatedRNAVar(i-1,j,k)+... %previously accumulated RNA
                dt*tempInterRNAProductionRate(i,j)-... %synthesized RNA
                dt*log(2)/Tm*tempInterAccumulatedRNAVar(i-1,j,k)+... %degraded RNA
                dt*DmVar(k)/(dx^2)*(tempInterAccumulatedRNAVar(i-1,j+1,k)+tempInterAccumulatedRNAVar(i-1,j-1,k)-... %diffused RNA
                2*tempInterAccumulatedRNAVar(i-1,j,k));

            end
        end
    end
end



%% Calculation of Protein from raw Accumulated mRNA

%RNAProductionRate(isnan(RNAProductionRate))=0;
rawProtein = zeros(numTimeBins,numAPBins);
dt = adjustedElapsedTime(2);
dx = 500/numAPBins;

rawProtein(1,:) = dt*rp*rawAccumulatedRNA(1,:);

for i=2:numTimeBins
	for j=1:numAPBins
        
        %First AP bin. No diffusion to the left.
		if j==1
			rawProtein(i,j) = rawProtein(i-1,j)+... %previously accumulated protein
            dt*rp*rawAccumulatedRNA(i,j)-... %synthesized protein
            dt*log(2)/Tp*rawProtein(i-1,j)+... %degraded protein
            dt*Dp/(dx^2)*(rawProtein(i-1,j+1)-rawProtein(i-1,j)); %diffused protein

        %Last AP bin. No diffusion to the right.
        elseif j==numAPBins
            rawProtein(i,j) = rawProtein(i-1,j)+... %previously accumulated protein
            dt*rp*rawAccumulatedRNA(i,j)-... %synthesized protein
            dt*log(2)/Tp*rawProtein(i-1,j)+... %degraded protein
            dt*Dp/(dx^2)*(rawProtein(i-1,j-1)-rawProtein(i-1,j)); %diffused protein

        %All other AP bins.
        else
            rawProtein(i,j) = rawProtein(i-1,j)+... %previously accumulated protein
            dt*rp*rawAccumulatedRNA(i,j)-... %synthesized protein
            dt*log(2)/Tp*rawProtein(i-1,j)+... %degraded protein
            dt*Dp/(dx^2)*(rawProtein(i-1,j+1)+rawProtein(i-1,j-1)-... %diffused protein
            2*rawProtein(i-1,j));

        end
    end
end


%% Plot Protein accumulation from raw data

color = jet(numTimeBins);

figure(8)
hold on

%(500/(numAPBins-9)):500/(numAPBins-9:500

for i=1:numTimeBins
	plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,rawProtein(i,:),'color',color(i,:))
end

title(['Accumulated Protein from Accumulated mRNA: Late nc 12 to nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('Accumulated protein (AU)')

%% Calculation of Protein from temporally interpolated Accumulated mRNA

%RNAProductionRate(isnan(RNAProductionRate))=0;
tempInterProtein = zeros(length(scaledT),numAPBins);
dt = scaledT(2);
dx = 500/numAPBins;

tempInterProtein(1,:) = dt*rp*tempInterAccumulatedRNA(1,:);
thing = tempInterAccumulatedRNA(1,:);

for i=2:length(scaledT)
	for j=1:numAPBins
        
        %First AP bin. No diffusion to the left.
		if j==1
			tempInterProtein(i,j) = tempInterProtein(i-1,j)+... %previously accumulated protein
            dt*rp*tempInterAccumulatedRNA(i,j)-... %synthesized protein
            dt*log(2)/Tp*tempInterProtein(i-1,j)+... %degraded protein
            dt*Dp/(dx^2)*(tempInterProtein(i-1,j+1)-tempInterProtein(i-1,j)); %diffused protein

        %Last AP bin. No diffusion to the right.
        elseif j==numAPBins
            tempInterProtein(i,j) = tempInterProtein(i-1,j)+... %previously accumulated protein
            dt*rp*tempInterAccumulatedRNA(i,j)-... %synthesized protein
            dt*log(2)/Tp*tempInterProtein(i-1,j)+... %degraded protein
            dt*Dp/(dx^2)*(tempInterProtein(i-1,j-1)-tempInterProtein(i-1,j)); %diffused protein

        %All other AP Bins.
        else
            tempInterProtein(i,j) = tempInterProtein(i-1,j)+... %previously accumulated protein
            dt*rp*tempInterAccumulatedRNA(i,j)-... %synthesized protein
            dt*log(2)/Tp*tempInterProtein(i-1,j)+... %degraded protein
            dt*Dp/(dx^2)*(tempInterProtein(i-1,j+1)+tempInterProtein(i-1,j-1)-... %diffused protein
            2*tempInterProtein(i-1,j));

        end
    end
end

%% Plot Protein accumulation from temporally interpolated mRNA

color = jet(length(scaledT));

figure(15)
hold on


for i=1:length(scaledT)
	plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,tempInterProtein(i,:),'color',color(i,:))
end

title(['Predicted Protein Gradient: nc 13 to ~40 min into nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('Protein (AU)')



%% Plot Protein from Spatially interpolated accumulated mRNA


% color = jet(numTimeBins);
% 
% figure(13)
% hold on
% 
% for i=1:numTimeBins
% 	plot(1:1:length(scaledX)-99,spatInterProtein(i,100:length(scaledX)),'color',color(i,:))
% end


%% Calculation of Protein from temp inter Accumulated mRNA for different Dp values and half life values


%RNAProductionRate(isnan(RNAProductionRate))=0;
DpVar = logspace(-1,2.5,10);
TpVar = logspace(0,2,10);
tempInterProteinVar = zeros(length(scaledT),numAPBins,length(DmVar),length(DpVar),length(TpVar));
dt = scaledT(2);
dx = 500/numAPBins;

for k=1:length(DmVar)
    for l=1:length(DpVar)
        for m=1:length(TpVar)
           tempInterProteinVar(1,:,k,l,m) = dt*rp*tempInterAccumulatedRNAVar(1,:,k);
           kp = DpVar(l)/(dx^2);

            for i=2:length(scaledT)
                for j=1:numAPBins
        
                    %First AP bin. No diffusion to the left.
                    if j==1
                        tempInterProteinVar(i,j,k,l,m) = tempInterProteinVar(i-1,j,k,l,m)+... %previously accumulated protein
                        dt*rp*tempInterAccumulatedRNAVar(i,j,k)-... %synthesized protein
                        dt*log(2)/TpVar(m)*tempInterProteinVar(i-1,j,k,l,m)+... %degraded protein
                        dt*kp*(tempInterProteinVar(i-1,j+1,k,l,m)-tempInterProteinVar(i-1,j,k,l,m)); %diffused protein

                    %Last AP bin. No diffusion to the right.
                    elseif j==numAPBins
                        tempInterProteinVar(i,j,k,l,m) = tempInterProteinVar(i-1,j,k,l,m)+... %previously accumulated protein
                        dt*rp*tempInterAccumulatedRNAVar(i,j,k)-... %synthesized protein
                        dt*log(2)/TpVar(m)*tempInterProteinVar(i-1,j,k,l,m)+... %degraded protein
                        dt*kp*(tempInterProteinVar(i-1,j-1,k,l,m)-tempInterProteinVar(i-1,j,k,l,m)); %diffused protein

                    %All other AP bins.
                    else
                        tempInterProteinVar(i,j,k,l,m) = tempInterProteinVar(i-1,j,k,l,m)+... %previously accumulated protein
                        dt*rp*tempInterAccumulatedRNAVar(i,j,k)-... %synthesized protein
                        dt*log(2)/TpVar(m)*tempInterProteinVar(i-1,j,k,l,m)+... %degraded protein
                        dt*kp*(tempInterProteinVar(i-1,j+1,k,l,m)+tempInterProteinVar(i-1,j-1,k,l,m)-... %diffused protein
                        2*tempInterProteinVar(i-1,j,k,l,m));

                    end
                end
            end
        end
    end
end



%% Measured Protein fluorescence removing background

backGroundProteinFluo = zeros(numTimeBins);
cutProteinFluo = zeros(numTimeBins,numAPBins);

for i = 1:numTimeBins
    backGroundProteinFluo(i) = proteinFluo(i,numAPBins);
    cutProteinFluo(i,:) = proteinFluo(i,:) - backGroundProteinFluo(i);
end

%% Temporal interpolation of measured Protein data

% totalTime = max(adjustedElapsedTime);
% scaledT = totalTime/numTimeBins/10/2:totalTime/numTimeBins/10:totalTime-totalTime/numTimeBins/10/2;

% totalTime = max(adjustedElapsedTime);
% dt = 0.1;
% scaledT = adjustedElapsedTime(1):dt*adjustedElapsedTime(2):totalTime;

tempInterMeasProtein = zeros(length(scaledT),numAPBins);

for j=1:numAPBins
    tempInterMeasProtein(:,j) = pchip(adjustedElapsedTime,cutProteinFluo(:,j),scaledT);
end

%% Syncing measured and accumulated at t=1

offsetMeasProtein = tempInterMeasProtein(1,:);

for i = 1:size(tempInterMeasProtein,1)
    
    tempInterMeasProtein(i,:) = tempInterMeasProtein(i,:) - offsetMeasProtein;
end




%% Plot measured protein fluorescence

color = jet(numTimeBins);

figure(9)
hold on


for i=1:numTimeBins
	plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,cutProteinFluo(i,:),'color',color(i,:))
end

title(['Measured Protein: nc 13 to ~40 min into nc 14'])
xlabel('Embryo Length (microns)')
ylabel('Protein Fluorescence (AU)')

%% Plot temporally interpolated measured Protein fluorescence

color = jet(length(scaledT));

figure(119)
hold on


for i=1:length(scaledT)
	plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,tempInterMeasProtein(i,:),'color',color(i,:))
end

title(['Measured Protein: nc 13 to ~40 min into nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('Accumulated protein (AU)')


length(scaledT)
% 
% %% Calculate matrix of boundary position/width for temporal interpolation
% 
% for t = linspace(152,336,4)
%     k = 1;
% %t = 241;
% length(scaledT)
% 
% [measBound, measWid] = singleTimeCalculateBoundaryFeatures(tempInterMeasProtein(:,:),scaledT,t,numAPBins,startAPBin,endAPBin,OGNumAPBins);
% boundaryNorm = zeros(length(DmVar),length(DpVar),length(TpVar));
% widthNorm = zeros(length(DmVar),length(DpVar),length(TpVar));
% widBoundNorm = zeros(length(DmVar),length(DpVar),length(TpVar));
% counter = 0;
% 
% for i=1:length(DmVar)
%     for j=1:length(DpVar)
%         for k=1:length(TpVar)
%         [predBound, predWid] = singleTimeCalculateBoundaryFeatures(tempInterProteinVar(:,:,i,j,k),scaledT,t,numAPBins,startAPBin,endAPBin,OGNumAPBins)
%             boundaryNorm(i,j,k) = abs(measBound-predBound);
%             widthNorm(i,j,k) = abs(measWid-predWid);
%             widBoundNorm = boundaryNorm .* widthNorm;
%             counter = counter +1
%         end
%     end
% end
% 
% boundaryMin = max(boundaryNorm);
% 
% % for k =1:length(DmVar)
% %     for l = 1:length(DpVar)
% %         for m = 1:length(TpVar)
% %             if boundaryNorm(k,l,m) < boundaryMin
% %                 boundaryMin = boundaryNorm(k,l,m);
% %                 DmBoundMin = k;
% %                 DpBoundMin = l;
% %                 TpBoundMin = m;
% %             end
% %         end
% %     end
% % end
% % 
% % MinBoundParam = [DmVar(DmBoundMin), DpVar(DpBoundMin), TpVar(TpBoundMin)]
% 
% bcutoff = 0.032
% boundaryNormSub = NaN(length(DmVar),length(DpVar),length(TpVar));
% 
% for i = 1:length(DmVar)
%     for j = 1:length(DpVar)
%         for k = 1:length(TpVar)
%             if boundaryNorm(i,j,k) <= bcutoff
%                 boundaryNormSub(i,j,k) = boundaryNorm(i,j,k);
%             end
%         end
%     end
% end
% 
% wcutoff = 0.05
% widthNormSub = NaN(length(DmVar),length(DpVar),length(TpVar));
% 
% for i = 1:length(DmVar)
%     for j = 1:length(DpVar)
%         for k = 1:length(TpVar)
%             if widthNorm(i,j,k) <= wcutoff
%                 widthNormSub(i,j,k) = widthNorm(i,j,k);
%             end
%         end
%     end
% end
% 
% %% Calculate matrix of least squares over each AP bin for one time point, for all values of Dm,Dp,Tp
% 
% lsNorm = zeros(length(DmVar),length(DpVar),length(TpVar));
% 
% 
% for k = 1:length(DmVar)
%     for l = 1:length(DpVar)
%         for m = 1:length(TpVar)
%             lsNorm(k,l,m) = lsfit(tempInterMeasProtein,tempInterProteinVar(:,:,k,l,m),t);
%         end
%     end
% end
% 
% lsMin = max(lsNorm);
% 
% for k =1:length(DmVar)
%     for l = 1:length(DpVar)
%         for m = 1:length(TpVar)
%             if lsNorm(k,l,m) < lsMin
%                 lsMin = lsNorm(k,l,m);
%                 DmSqMin = k;
%                 DpSqMin = l;
%                 TpSqMin = m;
%             end
%         end
%     end
% end
% 
% MinSqParam = [DmVar(DmSqMin), DpVar(DpSqMin), TpVar(TpSqMin)]
% 
% 
% %% Plot boundary position matrix as heat map
% 
% figure(11)
% [X,Y,Z] = ndgrid(DmVar,DpVar,TpVar);
% scatter3(X(:),Y(:),Z(:),40,boundaryNorm(:),'filled')
% set(gca,'Xscale','log','Yscale','log','Zscale','log')
% 
% title('Difference in boundary position between predicted and measured data')
% xlabel('mRNA Diffusion Coefficient [\mu m^2/min]')
% ylabel('Protein Diffusion Coefficient [\mu m^2/min]')
% zlabel('Protein Half-Life [min]')
% 
% colormap(jet(256));                                     % create and label the colorbar
% colorbar;
% M(k) = getframe;
% %cb.Label.String = '\Delta Boundary Position';
% 
% %% Plot boundary position subspace matrix as heat map
% 
% figure(122)
% [X,Y,Z] = ndgrid(DmVar,DpVar,TpVar);
% scatter3(X(:),Y(:),Z(:),40,boundaryNormSub(:),'filled')
% set(gca,'Xscale','log','Yscale','log','Zscale','log')
% 
% title('Subspace Difference in boundary position between predicted and measured data')
% xlabel('mRNA Diffusion Coefficient [\mu m^2/min]')
% ylabel('Protein Diffusion Coefficient [\mu m^2/min]')
% zlabel('Protein Half-Life [min]')
% 
% colormap(jet(256));                                     % create and label the colorbar
% colorbar;
% %cb.Label.String = '\Delta Boundary Position';
% 
% %% Plot boundary width subspace matrix
% 
% figure(123)
% [X,Y,Z] = ndgrid(DmVar,DpVar,TpVar);
% scatter3(X(:),Y(:),Z(:),40,widthNormSub(:),'filled')
% set(gca,'Xscale','log','Yscale','log','Zscale','log')
% 
% title('Difference in boundary width < 10% of embryo length')
% xlabel('mRNA Diffusion Coefficient [\mu m^2/min]')
% ylabel('Protein Diffusion Coefficient [\mu m^2/min]')
% zlabel('Protein Half-Life [min]')
% 
% % colormap(jet(256));                                     % create and label the colorbar
% % colorbar;
% %cb.Label.String = '\Delta Boundary Position';
% 
% %% Plot boundary width*position matrix as heat map
% 
% figure(30)
% [X,Y,Z] = ndgrid(DmVar,DpVar,TpVar);
% scatter3(X(:),Y(:),Z(:),40,widBoundNorm(:),'filled')
% set(gca,'Xscale','log','Yscale','log','Zscale','log')
% 
% title('bound width times bound')
% xlabel('mRNA Diffusion Coefficient [\mu m^2/min]')
% ylabel('Protein Diffusion Coefficient [\mu m^2/min')
% zlabel('Protein Half-Life [min]')
% 
% colormap(jet(256));                                     % create and label the colorbar
% colorbar;
% %cb.Label.String = '\Delta Boundary Width * \Delta Boundary Position';
% 
% %% Plot boundary width matrix as heat map
% 
% figure(20)
% [X,Y,Z] = ndgrid(DmVar,DpVar,TpVar);
% scatter3(X(:),Y(:),Z(:),40,widthNorm(:),'filled')
% set(gca,'Xscale','log','Yscale','log','Zscale','log')
% 
% title('Difference in boundary width between predicted and measured data')
% xlabel('mRNA Diffusion Coefficient [\mu m^2/min]')
% ylabel('Protein Diffusion Coefficient [\mu m^2/min]')
% zlabel('Protein Half-Life [min]')
% 
% colormap(jet(256));                                     % create and label the colorbar
% colorbar;
% 
% N(k) = getframe;
% %cb.Label.String = '\Delta Boundary Width';
% 
% k = k + 1;
% end
% 
% figure(333)
% movie(M,1)
% 
% figure(232)
% movie(N,1)
% 
% %% Plot least square difference as heat map between measured and predicted protein with temporal interpolation
% 
% figure(17)
% [X,Y,Z] = ndgrid(DmVar,DpVar,TpVar);
% scatter3(X(:),Y(:),Z(:),40,lsNorm(:),'filled')
% set(gca,'Xscale','log','Yscale','log','Zscale','log')
% 
% xlabel('mRNA Diffusion Coefficient')
% ylabel('Protein Diffusion Coefficient')
% zlabel('Protein Half-Life')
% 
% colormap(jet(256));                                     % create and label the colorbar
% colorbar;
% %cb.Label.String = 'Least squares predicted/measured Protein';
% 
%% Plot meas and pred at given time point

figure(3000)
hold on



plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,tempInterMeasProtein(160,:))
plot(startAPBin/OGNumAPBins:1/OGNumAPBins:endAPBin/OGNumAPBins,tempInterProteinVar(160,:,5,5,5))


title(['Predicted and measured protein ~30 min into nc 14'])
xlabel('Fraction of Embryo Length')
ylabel('Accumulated protein (AU)')

% %% Plot the Average mRNA and Protein over time
% % a = 50; %constant for scaling of protein
% % 
% % figure(11)
% % hold on
% % plot(adjustedElapsedTime,AccumulatedmRNA,'LineWidth',3)
% % plot(adjustedElapsedTime,a*cutProteinFluoAll,'LineWidth',3)
% % 
% % title('Averaged protein level over time')
% % xlabel('Time (min)')
% % ylabel('Fluorescence (AU)')
% % legend('Protein (NB signal)')
% 
% %% Comparing inferred Protein with measured Protein
% % 
% % b = 10; %protein scaling 
% % 
% % figure(12)
% % hold on
% % for k=1:length(DpVar)
% %     DpVarVect= DpVar(k)*ones(size(rawProteinVar(2)));
% %     
% %     for i=20:floor(numTimeBins/20):20*floor(numTimeBins/20)
% %         scatter3(rawProteinVar(i,:,k),cutProteinFluo(i,:),DpVarVect) %DpVar must be length of rawProteinVar
% %     end
% %     
% % end
% %         


%% Parameter fitting estimation with boundary position and width

tempInterMeasProtein(1,:)

length(scaledT)
save('E:\EvanM\RandomVariables\tempInterMeasProtein.mat','tempInterMeasProtein')
save('E:\EvanM\RandomVariables\scaledT.mat','scaledT')
save('E:\EvanM\RandomVariables\numAPBins.mat','numAPBins')
save('E:\EvanM\RandomVariables\startAPBin.mat','startAPBin')
save('E:\EvanM\RandomVariables\endAPBin.mat','endAPBin')
save('E:\EvanM\RandomVariables\OGNumAPBins.mat','OGNumAPBins')
save('E:\EvanM\RandomVariables\tempInterAccumulatedRNA.mat','tempInterAccumulatedRNA')
save('E:\EvanM\RandomVariables\RNAProductionRate.mat','RNAProductionRate')
save('E:\EvanM\RandomVariables\rawAccumulatedRNA.mat','rawAccumulatedRNA')
save('E:\EvanM\RandomVariables\tempInterRNAProductionRate.mat','tempInterRNAProductionRate')
save('E:\EvanM\RandomVariables\adjustedElapsedTime.mat','adjustedElapsedTime')
save('E:\EvanM\RandomVariables\tempExamp.mat','tempExamp')
save('E:\EvanM\RandomVariables\tempInterProteinVar.mat','tempInterProteinVar')



end