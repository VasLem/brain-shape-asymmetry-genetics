function [CHR,POS,A1,A2,RSID,nSNP] = importBIM(bfile)
% IMPORTPLINKDATA Import Plink binary data
%   IMPORTPLINKDATA(obj, bfile) import plink binary data
%   (bfile.bed, bfile.bim, bfile.fam) into obj
    filename = [bfile,'.bim'];
    formatSpec = '%s%s%f%f%s%s';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    chr = dataArray{1};
    CHR = str2double(chr);
    CHR(strcmp(chr,'X')) = 23;
    CHR(strcmp(chr,'Y')) = 24;
    CHR(strcmp(chr,'XY')) = 25;
    CHR(strcmp(chr,'MT')) = 26;
    RSID = dataArray{2};
    POS = dataArray{4};
    A1 = dataArray{5};
    A2 = dataArray{6};
    nSNP = length(RSID);
end