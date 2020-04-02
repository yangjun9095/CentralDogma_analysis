function boundFeaturesTestScript

tempInterMeasProteinFile = load('E:\EvanM\RandomVariables\tempInterMeasProtein','tempInterMeasProtein');
tempInterMeasProtein = tempInterMeasProteinFile.tempInterMeasProtein;
scaledTFile = load('E:\EvanM\RandomVariables\scaledT','scaledT');
scaledT = scaledTFile.scaledT;
scaledT = scaledT(:);
numAPBinsFile = load('E:\EvanM\RandomVariables\numAPBins','numAPBins');
numAPBins = numAPBinsFile.numAPBins;
startAPBinFile = load('E:\EvanM\RandomVariables\startAPBin','startAPBin');
startAPBin = startAPBinFile.startAPBin;
endAPBinFile = load('E:\EvanM\RandomVariables\endAPBin','endAPBin');
endAPBin = endAPBinFile.endAPBin;
OGNumAPBinsFile = load('E:\EvanM\RandomVariables\OGNumAPBins','OGNumAPBins');
OGNumAPBins = OGNumAPBinsFile.OGNumAPBins;
tempInterProteinVarFile = load('E:\EvanM\RandomVariables\tempInterProteinVar','tempInterProteinVar');
tempInterProteinVar = tempInterProteinVarFile.tempInterProteinVar;

length(tempInterMeasProtein)
length(scaledT)

singleTimeCalculateBoundaryFeatures(tempInterMeasProtein,scaledT,800,numAPBins,startAPBin,endAPBin,OGNumAPBins);
singleTimeCalculateBoundaryFeatures(tempInterProteinVar(:,:,1,1,1),scaledT,800,numAPBins,startAPBin,endAPBin,OGNumAPBins);
%CalculateBoundaryFeatures(tempInterMeasProtein,scaledT,4,18,38);