%% Plot vasa-eGFP intensities compared to eGFP23

% The data is acquired from each line by measruing the fluorescence of
% 4~5 embryos, then averaging the fluorescence from the rectangular region.
% I need to double check for one or two datasets just as a sanity check.
% The lines could be either hetero or homozygous, we could only know when
% there are heterogenous populations.
% I need to measure again for 
Lines = {'eGFP23','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15'};
% Hets = [2,6,7,8,9];
MeanFluoLines = [1
1.074909206
0.566998263
0.653124288
0.577261961
0.733901784
0.441654824
0.560760189
0.879414445
0.552254998
0.601812017
0.592834633
0.631916103
0.471663998
0.799270905
0.938730175
]; % From Meghan's dataset in Excel file

hold on
plot(1:16,MeanFluoLines,'o','MarkerFaceColor','b')

title('Dosage of vasa-eGFP')
xlabel('Lines')
ylabel('dosage normalized to eGFP23')

%set(gca,'XtickLabel',Lines)
xtick = (1:16)
xticklabels(Lines)
%xlim([-1 10]);
%hold off
%% Plot using Simon's code
%%
Lines = {'1','2','3A','3B','4','5','6','7','8','9','10','11','12','13A','13B','14A',''};
Hets = [2,6,7,8,9];
MeanFluoLines = zeros(1,length(Lines));
hold on
for line = 1:length(Lines)
    MeanFluoColumn = Data(:,line*3);
    EmbryoMeans =[];
    for embryo = 1:8
        SingleEmbryoMean = nanmean(MeanFluoColumn([embryo:embryo+3]));
        EmbryoMeans = [EmbryoMeans SingleEmbryoMean]
    end
        
    if any(Hets==line)
        EmbryoMeans = EmbryoMeans*2;
        MeanFluoLines(line) = MeanFluoLines(line)*2;
    end
    
    MeanEmbryoMeans = nanmean(EmbryoMeans);
    MeanFluoLines(line) = MeanEmbryoMeans;

    plot(line,MeanFluoLines(line),'o','MarkerFaceColor','b')
    plot(line*ones(1,8),EmbryoMeans,'ro')
end
%xticks(Lines)
xlim([-1 10]);
hold off

