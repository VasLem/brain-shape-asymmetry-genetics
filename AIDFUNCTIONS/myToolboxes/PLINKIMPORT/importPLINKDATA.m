function obj = importPLINKDATA(bfile)
% IMPORTPLINKDATA Import Plink binary data
%   IMPORTPLINKDATA(obj, bfile) import plink binary data
%   (bfile.bed, bfile.bim, bfile.fam) into obj
filename = [bfile,'.bim'];
formatSpec = '%s%s%f%f%s%s';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
chr = dataArray{1};
obj.CHR = str2double(chr);
obj.CHR(strcmp(chr,'X')) = 23;
obj.CHR(strcmp(chr,'Y')) = 24;
obj.CHR(strcmp(chr,'XY')) = 25;
obj.CHR(strcmp(chr,'MT')) = 26;
obj.RSID = dataArray{2};
obj.POS = dataArray{4};
obj.A1 = dataArray{5};
obj.A2 = dataArray{6};
filename = [bfile,'.fam'];
formatSpec = '%s%s%s%s%f%s';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
fclose(fileID);
obj.FID = dataArray{1};
obj.IID = dataArray{2};
obj.PID = dataArray{3};
obj.MID = dataArray{4};
obj.Sex = dataArray{5};
fid = fopen([bfile '.bed'],'rb');
bin = fread(fid,inf,'uint8=>uint8');
fclose(fid);
obj.nSamples = length(obj.FID);
obj.nSNPs = length(obj.CHR);
L = ceil(obj.nSamples/4);
GENO = reshape(bin(4:end),[L obj.nSNPs]);
obj.GENO = UnpackGeno_(GENO,obj.nSamples);
end