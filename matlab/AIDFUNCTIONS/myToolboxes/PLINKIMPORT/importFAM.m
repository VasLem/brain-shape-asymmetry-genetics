function [FID,IID,PID,MID,Sex,nSamples] = importFAM(bfile)
% IMPORTPLINKDATA Import Plink binary data
%   IMPORTPLINKDATA(obj, bfile) import plink binary data
%   (bfile.bed, bfile.bim, bfile.fam) into obj
    filename = [bfile,'.fam'];
    formatSpec = '%s%s%s%s%f%s';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    FID = dataArray{1};
    IID = dataArray{2};
    PID = dataArray{3};
    MID = dataArray{4};
    Sex = dataArray{5};
    nSamples = length(FID);
end