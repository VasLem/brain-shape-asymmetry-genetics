function [Genotypes,CALL,MAF] = importBED(bfile,nSNP,nSamples)
% IMPORTPLINKDATA Import Plink binary data
%   IMPORTPLINKDATA(obj, bfile) import plink binary data
%   (bfile.bed, bfile.bim, bfile.fam) into obj
    fid = fopen([bfile '.bed'],'rb'); 
    bin = fread(fid,inf,'uint8=>uint8');
    fclose(fid);
    if ~ismac
        bin = bin(4:end);
        bin = de2bi(bin,8);
        geno = zeros(size(bin,1),4,'int8');
        geno(:,1) = bi2de(bin(:,[1 2]));
        geno(:,2) = bi2de(bin(:,[3 4]));
        geno(:,3) = bi2de(bin(:,[5 6]));
        geno(:,4) = bi2de(bin(:,[7 8]));
        geno = geno';
        nBytes = ceil((nSamples*2)/8);
        geno = reshape(geno,nBytes*4,nSNP);
        Genotypes = geno;
        Genotypes(geno==3) = 2;
        Genotypes(geno==1) = -1;
        Genotypes(geno==2) = 1;
        Genotypes = Genotypes(1:nSamples,:);
    else
        L = ceil(nSamples/4);
        geno = reshape(bin(4:end),[L nSNP]);
        Genotypes = UnpackGeno_(geno,nSamples);
    end
    if nargout==1, return; end
%% CALCULATING CALLRATE AND MAF and select accordingly
    Ne = nSamples-sum(Genotypes==-1,1);
    CALL = Ne./nSamples;
    tmp = Genotypes;
    tmp(Genotypes==-1)=0;
    Na = sum(tmp);
    p = Na./Ne./2;
    MAF = min([p;1-p],[],1);
end
