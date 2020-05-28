% projectName = 'CentralDogma'
% Last updated : 5/22/2020
% Step1. Process the datasets in different ways with options so that it's
% ready to be plugged into the next steps for calculating the cytoplasmic
% mRNA, and Protein.

% Inputs : 1) single embryo (Prefix), 2) DataType (LoadMS2Sets) or
% AveragedDatasets, 3) individual embryos (optional...)

% Outputs : 1) MS2 Particles : Time, MeanVectorAP, SDVectorAP, SEVectorAP, NParticlesAP, etc.
%          2) Nuclear fluo :  MeanVectorAP, SDVectorAP, SEVectorAP, NParticlesAP, etc.
% Save those into two different structures, cp.mat, cn.mat
% (compiledparticles, and compilednuclei)

function process_datasets_multipleEmbryos(DataType,varargin)
end
