function [NMISSsnp, NGENOsnp, FMISsnp] = importLMISS(bfile)
% function [CHR, RSID,NMISSsnp, NGENOsnp, FMISsnp] = importLMISS
filename = [bfile,'.lmiss']; % import missing-file 
formatSpec = '%s%s%f%f%f';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'HeaderLines',1,'ReturnOnError', false);
fclose(fileID);
% chr = dataArray{1};
% CHR = str2double(chr);
% CHR(strcmp(chr,'X')) = 23;l
% CHR(strcmp(chr,'Y')) = 24;
% CHR(strcmp(chr,'XY')) = 25;
% CHR(strcmp(chr,'MT')) = 26;
% RSID = dataArray{2};
NMISSsnp = dataArray{3};    % Number of individuals missing this SNP
NGENOsnp = dataArray{4};    % Number of non-obligatory (certain genotypes were never attempted- they were obligatory missing) missing genotypes
FMISsnp = dataArray{5};     % Proportion of sample missing for this SNP
end
