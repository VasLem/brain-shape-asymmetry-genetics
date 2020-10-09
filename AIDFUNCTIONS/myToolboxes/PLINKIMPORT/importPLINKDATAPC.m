function out  = importPLINKDATA(bfile,selmaf,selcallrate,seliid)
% INIT    
    if nargin<3,selcallrate = 0;end
    if nargin<2,selmaf = 0;end
    if nargin<4, seliid = []; end
%% Read SNP info
    filename = [bfile,'.bim'];
    %delimiter = '\t';
    %formatSpec = '%*s%s%*s%f%*s%*s%[^\n\r]';
    formatSpec = '%f%s%f%f%s%s';
    fileID = fopen(filename,'r');
    %dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    CHR = dataArray{1};
    RSID = dataArray{2};
    %out.RS = cleanUpRSArray(RSID);
    POS = dataArray{4};
    A1 = dataArray{5};
    A2 = dataArray{6};
%% Read Individual info
    %delimiter = ' ';
    filename = [bfile,'.fam'];
    %formatSpec = '%s%s%*s%*s%*s%*s%[^\n\r]';
    formatSpec = '%s%s%s%s%s%s';
    fileID = fopen(filename,'r');
    %dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    dataArray = textscan(fileID, formatSpec, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    FID = dataArray{1};
    IID = dataArray{2};
    N = length(FID);
    if isempty(seliid)
       ind12 = 1:N;ind21 = 1:N;
    else
       [ind12,ind21] = vlookupFast(seliid,IID);
    end
    IID = IID(ind12);
    FID = FID(ind12);
%% Read Genotypes
    nSNP = length(RSID);
    if ~ismac
        fid = fopen([bfile '.bed']);
        bin = fread(fid,inf,'uint8=>uint8');
        fclose(fid);
        bin = bin(4:end);
        bin = de2bi(bin,8);
        geno = zeros(size(bin,1),4,'int8');
        geno(:,1) = bi2de(bin(:,[1 2]));
        geno(:,2) = bi2de(bin(:,[3 4]));
        geno(:,3) = bi2de(bin(:,[5 6]));
        geno(:,4) = bi2de(bin(:,[7 8]));
        geno = geno';
        nBytes = ceil((N*2)/8);
        geno = reshape(geno,nBytes*4,nSNP);
        Genotypes = geno;
        Genotypes(geno==3) = 2;
        Genotypes(geno==1) = -1;
        Genotypes(geno==2) = 1;
        Genotypes = Genotypes(1:N,:);
    else
        geno = read_genotypes(bfile);
        Genotypes = unpack_geno(geno,N,feature('numcores'));
    end
    Genotypes = Genotypes(ind12,:);
    N = length(FID);
%% CALCULATING CALLRATE AND MAF and select accordingly
    Ne = N-sum(Genotypes==-1,1);
    CALL = Ne./N;
    tmp = Genotypes;
    tmp(Genotypes==-1)=0;
    Na = sum(tmp);
    p = Na./Ne./2;
    MAF = min([p;1-p],[],1);
    ind = find(CALL>=selcallrate&MAF>=selmaf);
    out.CHR = CHR(ind);
    out.RSID = RSID(ind);
    out.POS = POS(ind);
    out.SNP = Genotypes(:,ind);
    out.CALL = CALL(ind);
    out.MAF = MAF(ind);
    out.A1 = A1(ind);
    out.A2 = A2(ind);
    out.FID = FID;
    out.IID = IID;
    out.nSNP = length(out.POS);
    out.ind12 = ind12;
    out.ind21 = ind21;
end



