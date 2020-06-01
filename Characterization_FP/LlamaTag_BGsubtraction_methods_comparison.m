% projectName = 'CentralDogma'
% Last updated : 5/29/2020
% Description : For a given Prefix (single embryo), subtract the background
% fluo from the nuclear fluo, for the LlamaTag datasets. 
% Note that the cytoplasmic/nuclear ratio of TF can be TF specific.

% Inputs : 1) single embryo (Prefix), 2) DataType (LoadMS2Sets) or
% AveragedDatasets, 3) individual embryos (optional...)

%% Step1. Subtract the nuc fluo from the most posterior AP bins

%% Step2. 