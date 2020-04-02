% Use the CorrelatemRNAProtein function to plot different combinations of
% Dm, Dp, and gammaP

[AccumulatedmRNA,Protein,RNAPLoadingRate3,InferredProtein, NBProteinFluo] = ...
                CorrelatemRNAProtein_function(Prefix1,Prefix2,Dm,Dp,Tm,Tp,rp)
            
%% Correlation for different Dm (mRNA diffusion coefficient)
Dm = [0 0.1 1 10 100];
Dp = 5;
Tm = 60;
Tp = 50;
rp = 2;
mRNADiffusion={};

for i=1:length(Dm)
    clear InferredProtein
    clear NBProteinFluo
    [AccumulatedmRNA,Protein,RNAPLoadingRate3,InferredProtein, NBProteinFluo] = ...
                CorrelatemRNAProtein_function(Prefix1,Prefix2,Dm(i),Dp,Tm,Tp,rp)
    mRNADiffusion(i).InferredProtein = InferredProtein;
    mRNADiffusion(i).NBProteinFluo = NBProteinFluo;
    mRNADiffusion(i).Dm = Dm(i);
    
end 

%% Plot the correlation for Dm
%ap = 11;
% Color code for time points
iStart=33; %Start time
iEnd=93; %End time
colormap(jet(256));
cmap=colormap ;
Color=cmap(round(((iStart:iEnd)-iStart)/(iEnd-iStart)*255)+1,:);

hold on
for i=1:length(Dm)
    for tpoint = 33:15:93
        clear InferredProtein
        clear NBProteinFluo
        clear Z

        InferredProtein = mRNADiffusion(i).InferredProtein;
        NBProteinFluo = mRNADiffusion(i).NBProteinFluo;
        Z = Dm(i)*ones(size(InferredProtein));
        Z = log10(Z);
        scatter3(InferredProtein(tpoint,:),Z(tpoint,:),NBProteinFluo(tpoint,:))
        
        title('Prediction vs Measured protein')
        xlabel('Predicted Protein (AU)')
        ylabel('log(Dm)')
        zlabel('Measured Protein (AU)')
        %zlim([])
        pause
    end
end
%% Correlation for different Dp (protein diffusion coefficient)
Dm = 0;
Dp = [0 0.1 1 10 100];
Tm = 60;
Tp = 50;
rp = 2;
ProteinDiffusion={};

for i=1:length(Dp)
    clear InferredProtein
    clear NBProteinFluo
    [AccumulatedmRNA,Protein,RNAPLoadingRate3,InferredProtein, NBProteinFluo] = ...
                CorrelatemRNAProtein_function(Prefix1,Prefix2,Dm,Dp(i),Tm,Tp,rp)
    ProteinDiffusion(i).InferredProtein = InferredProtein;
    ProteinDiffusion(i).NBProteinFluo = NBProteinFluo;
    ProteinDiffusion(i).Dp = Dp(i);
    
end 

%% Plot the correlation for Dp
%ap = 12;
hold on
for i=1:length(Dp)
    for tpoint = 33:15:93
        clear InferredProtein
        clear NBProteinFluo
        clear Z

        InferredProtein = ProteinDiffusion(i).InferredProtein;
        NBProteinFluo = ProteinDiffusion(i).NBProteinFluo;
        Z = Dp(i)*ones(size(InferredProtein));
        Z = log10(Z);
        scatter3(InferredProtein(tpoint,:),Z(tpoint,:),NBProteinFluo(tpoint,:))
        
        title('Prediction vs Measured protein')
        xlabel('Predicted Protein (AU)')
        ylabel('log(Dp)')
        zlabel('Measured Protein (AU)')
        %zlim([])
        pause
    end
end
%% Correlation for different Tp (protein half-life)
Dm = 0;
Dp = 5;
Tm = 60;
Tp = [1,5,10,50,100,1000];
rp = 2;
ProteinLifeTime={};

for i=1:length(Tp)
    clear InferredProtein
    clear NBProteinFluo
    [AccumulatedmRNA,Protein,RNAPLoadingRate3,InferredProtein, NBProteinFluo] = ...
                CorrelatemRNAProtein_function(Prefix1,Prefix2,Dm,Dp,Tm,Tp(i),rp)
    ProteinLifeTime(i).InferredProtein = InferredProtein;
    ProteinLifeTime(i).NBProteinFluo = NBProteinFluo;
    ProteinLifeTime(i).Tp = Tp(i);
    
end 

%% Plot the correlation for Tp
%ap = 12;
hold on
for i=1:length(Tp)
    for tpoint = 33:15:93
        clear InferredProtein
        clear NBProteinFluo
        clear Z

        InferredProtein = ProteinLifeTime(i).InferredProtein;
        NBProteinFluo = ProteinLifeTime(i).NBProteinFluo;
        Z = Tp(i)*ones(size(InferredProtein));
        Z = log10(Z);
        scatter3(InferredProtein(tpoint,:),Z(tpoint,:),NBProteinFluo(tpoint,:))
        
        title('Prediction vs Measured protein')
        xlabel('Predicted Protein (AU)')
        ylabel('log(Tp)')
        zlabel('Measured Protein (AU)')
        %zlim([])
        pause
    end
end

%% save the result variables
save(['D:\Data\YangJoon\LivemRNA\Data\DynamicsResults\2017-03-02-Hb-P2P-MS2V5-MCP-GFP\ComparePredictiontoMeasurement.mat'],...
    'mRNADiffusion','ProteinDiffusion','ProteinLifeTime')