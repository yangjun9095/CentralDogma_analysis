function compareOffsets(Prefix1,Prefix2)

%%Load Data

compPart1 = dummy1
compPart2 = dummy2

offset1 = compPart1.Off
offset2 = compPart2.Off


%%Compare Offsets

%13 start of offset1
%38 start of offset2


something1 = nan(13,300);
something2 = nan(13,500);

adder1 = ones(13);
adder2 = ones(13);

for i = 1:length(compPart1.OriginalParticle)
    for j = 1:length(compPart1.Frame(i))
        something1(compPart1.Frame(i,j),adder1(compPart1.Frame(i,j))) = compPart1.Off(i,j);
        adder1(compPart1.Frame(i,j)) = adder1(compPart1.Frame(i,j)) + 1;
    end
end

for i = 1:length(compPart2.OriginalParticle)
    for j = 1:length(compPart2.Frame(i))
        something2(compPart2.Frame(i,j),adder2(compPart2.Frame(i,j))) = compPart2.Off(i,j);
        adder2(compPart2.Frame(i,j)) = adder2(compPart2.Frame(i,j)) + 1;
    end
end

means1 = nanmean(something1);
means2 = nanmean(something2);

stdDev1 = zeros(13);
stdDev2 = zeros(13);

for i = 1:13
    stdDev1(i) = std(something1(i,:));
    stdDev2(i) = std(something2(i,:));
end

difOfMeans = means1 - means2;


%%Plot comparisons

