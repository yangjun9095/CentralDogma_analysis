% New analysis code, example
%% Load the dataset
mRNA = load('CompiledParticles.mat')
Protein = load('CompiledNuclei.mat')

%% Protein Boundary Calculation along the time 
%Extract the MeanVectorAP data as Protein at diff.time points, at different AP bins
protein=Protein.MeanVectorAP;
ncTime=Protein.ElapsedTime;
nc14=Protein.nc14;

for i=nc14:length(ncTime)
    Pos(i)=sum(protein(i,:)>max(protein(i,:))*0.5);
    Bpos(i)=Pos(i)*0.025;
end

plot(ncTime(nc14:end)-ncTime(nc14),Bpos(nc14:end))