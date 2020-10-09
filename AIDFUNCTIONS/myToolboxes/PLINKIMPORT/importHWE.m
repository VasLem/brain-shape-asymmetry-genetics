function [GTCount, Ohet, Ehet, P]  = importHWE(bfile)
%function [CHR, RSID, Geno, Ohet, Ehet, P]  = importHWE(bfile)
filename = [bfile,'.hwe']; % import hwe-file 
formatSpec = '%s%s%s%s%s%s%f%f%f';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'HeaderLines',1,'ReturnOnError', false);
fclose(fileID);
% chr = dataArray{1};
% CHR = str2double(chr);
% CHR(strcmp(chr,'X')) = 23;
% CHR(strcmp(chr,'Y')) = 24;
% CHR(strcmp(chr,'XY')) = 25;
% CHR(strcmp(chr,'MT')) = 26;
% RSID = dataArray{2};
GTCount = dataArray{6};    % Genotype counts
Ohet = dataArray{7};    % Observed heterozygosity
Ehet = dataArray{8};    % Expected heterozygosity
P = dataArray{9};       % HWE test P-value 
end
