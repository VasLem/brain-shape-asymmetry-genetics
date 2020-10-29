function [MISSPHENO, NMISS, NGENO, FMISS] = importIMISS(bfile)
%function [FID,IID,MISSPHENO, NMISS, NGENO, FMISS] = importIMISS(bfile)
filename = [bfile,'.imiss']; % import missing-file 
formatSpec = '%s%s%s%f%f%f';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'HeaderLines',1,'ReturnOnError', false);
fclose(fileID);
% FID = dataArray{1};
% IID = dataArray{2};
MISSPHENO = dataArray{3};   % Missing phenotype? (Y/N)
NMISS = dataArray{4};       % Number of missing SNPs
NGENO = dataArray{5};       % Number of non-obligatory missing genotypes - denominator for the genotyping rate calculations
FMISS = dataArray{6};       % Proportion of missing SNPs 
end
