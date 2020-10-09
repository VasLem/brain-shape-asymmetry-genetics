function [MAFplink, NCHROBS] = importMAF(bfile)
%function [CHR, SNP, A1, A2, MAFplink, NCHROBS] = importMAF(bfile)
filename = [bfile,'.frq']; 
formatSpec = '%s%s%s%s%s%f'; %import freq-file
fileID = fopen(filename,'r'); 
dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true,'EmptyValue' ,NaN, 'HeaderLines',1,'ReturnOnError', false); 
fclose(fileID);
% chr = dataArray{1};
% CHR = str2double(chr);
% CHR(strcmp(chr,'X')) = 23;
% CHR(strcmp(chr,'Y')) = 24;
% CHR(strcmp(chr,'XY')) = 25;
% CHR(strcmp(chr,'MT')) = 26;
% RSID = dataArray{2};
% A1 = dataArray{3};
% A2 = dataArray{4};
MAFplink = dataArray{5};    % Minor Allele Frequency
NCHROBS = dataArray{6};     % Non-missing allele count 
bad = strcmp(MAFplink,'NA');
MAFplink(bad) = {'0'};
MAFplink = str2double(MAFplink);
end


