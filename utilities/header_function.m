function [RawResultsPath, DataPath, FigureRoot] =   header_function(DropboxFolder) 
    RawResultsPath = [DropboxFolder, filesep, 'CentralDogmaResults\'];
    DataPath = [DropboxFolder, filesep, 'CentralDogmaResults\'];
    FigureRoot = [DropboxFolder, filesep, 'Garcia Lab\Figures\CentralDogmaFigures\'];
    %mkdir(FigureRoot);
end