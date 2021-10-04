function GWAS = LoadBRAINGWASScriptCHIv2SWH
OUTPUTPATH = '/IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/GWASOUTPUT2020/SWH_PA80_20200513/';
INPUTPATH = '/IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/SNPJOBSv500/';
   
in = load('/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/SWH/AUX_SWH_PA80_20200513.mat');
IID = in.DB.IID;
nSegments = in.PA.nSegments;
%nSubjects = length(IID);
nSNPperJOB = 500;
clear in;
%%
try
    parpool('LocalSingle',20);
catch
end
%% PREPARE DATA CONTAINERS
cd(INPUTPATH);files = dir('**/*.mat');% retrieve all job folders recursively
nJOBS = length(files); 

RS = cell(nSNPperJOB,nJOBS);
POS = zeros(nSNPperJOB,nJOBS,'uint32');
CHR = zeros(nSNPperJOB,nJOBS,'uint8');
P = zeros(nSNPperJOB,nSegments,nJOBS,'uint32');% p value 
JOBIID = zeros(nSNPperJOB,nJOBS,'uint16');
%
A1 = cell(nSNPperJOB,nJOBS);
A2 = cell(nSNPperJOB,nJOBS);
MAF = zeros(nSNPperJOB,nJOBS,'uint8');
N = zeros(nSNPperJOB,2,nJOBS,'uint16');
CHI= zeros(nSNPperJOB,nSegments,nJOBS,'single');
SIGN = zeros(nSNPperJOB,nSegments,nJOBS,'int8');

factor = 10000;    
%% HOUSEHOLDING, making sure that the JOBS will be loaded in order
CHRID = zeros(1,nJOBS);
JOBID = 1:nJOBS;
for i=1:nJOBS
    fpart = split(files(i).name,["_","."]); 
    CHRID(i) = str2double(fpart{3,1});
end
CHRID = sort(CHRID);
%% LOAD RESULTS
[path,ID] = setupParForProgress(nJOBS);tic;
parfor i= 1:nJOBS
    % i=1;
    inputfile = [INPUTPATH 'JOB_' num2str(JOBID(i)) '_' num2str(CHRID(i)) '.mat'];
    %outputfile = [OUTPUTPATH 'resJOB_' num2str(JOBID(i)) '_' num2str(CHRID(i)) '.mat'];
    outputfile = [OUTPUTPATH 'JOB_' num2str(JOBID(i)) '_' num2str(CHRID(i)) '.mat'];
    if isempty(outputfile), continue; end
    in = load(inputfile);
    JOB = in.JOB;
    nSNP = length(JOB.POS);
    [ind12,~] = vlookupFast(IID,JOB.IID);%
    JOB.SNP = JOB.SNP(ind12,:); % extract SNP data from overlapping subjects   
    n = sum(JOB.SNP>=0,1);
    n = n(:);
    maf = zeros(nSNP,1);
    a1 = cell(nSNP,1);
    a2 = cell(nSNP,1);
    for s=1:1:nSNP
        % s=1;
        geno = JOB.SNP(:,s);  
        [geno,~] = codetoM(geno);
        a1{s} = JOB.A1{s};
        a2{s} = JOB.A2{s};
        maf(s) = getMAF(geno);
        JOB.SNP(:,s) = geno;
    end
    JOBIID(:,i) = uint16(i);
    if nSNP==nSNPperJOB      
      % SNP(:,:,i) = JOB.SNP';
       RS(:,i) = in.JOB.RSID(:);
       A1(:,i) = a1(:);
       A2(:,i) = a2(:);
       POS(:,i) = uint32(in.JOB.POS(:));
       CHR(:,i) = in.JOB.CHRID(:);
       N(:,i) = uint16(n);
       MAF(:,i) = uint8(100*maf(:));
     else
       tmpRS = cell(nSNPperJOB,1);tmpRS(1:nSNP) = in.JOB.RSID;RS(:,i) = tmpRS(:);
       tmpA1 = cell(nSNPperJOB,1);tmpA1(1:nSNP) = a1;A1(:,i) = tmpA1(:);
       tmpA2 = cell(nSNPperJOB,1);tmpA2(1:nSNP) = a2;A2(:,i) = tmpA2(:);
       tmpPOS = zeros(nSNPperJOB,1);tmpPOS(1:nSNP) = in.JOB.POS;POS(:,i) = uint32(tmpPOS(:));
       tmpCHR = zeros(nSNPperJOB,1);tmpCHR(1:nSNP) = in.JOB.CHRID;CHR(:,i) = uint8(tmpCHR(:));
       tmpMAF = zeros(nSNPperJOB,1);tmpMAF(1:nSNP) = uint8(100*maf(:));MAF(:,i) = tmpMAF;
       tmpN = zeros(nSNPperJOB,1);tmpN(1:nSNP) = uint16(n);N(:,i) = tmpN;
       %tmpSNP = zeros(nSNPperJOB,nSubjects,'int8');tmpSNP(1:nSNP,:) = JOB.SNP';SNP(:,:,i) = tmpSNP;
    end
    in = load(outputfile);
    if nSNP==nSNPperJOB
       pvalues = -log10(in.res.pvalues);
       chi = in.res.CHI2;
       effectsign = int8(sign(in.res.A));
    else
       pvalues = nan*zeros(nSNPperJOB,size(in.res.pvalues,2));
       pvalues(1:nSNP,:) = -log10(in.res.pvalues);
       chi = nan*zeros(nSNPperJOB,size(in.res.CHI2,2));
       chi(1:nSNP,:) = in.res.CHI2;
       effectsign = nan*zeros(nSNPperJOB,size(in.res.CHI2,2));
       effectsign(1:nSNP,:) = int8(sign(in.res.A));
    end            
    P(:,:,i) = uint32(pvalues*factor);
    CHI(:,:,i) = chi;
    SIGN(:,:,i) = effectsign;
    parfor_progress; 
end
closeParForProgress(path,ID);toc;
    %% REDUCING DATA
    TEST = sum(POS,1);
    FinishedJOBs = find(TEST);
    nFinishedJOBs = length(FinishedJOBs);
    POS = POS(:,FinishedJOBs);POS = POS(:);index = find(POS);POS = POS(index);
    nSNP = length(POS);
    disp([num2str(nFinishedJOBs) ' Jobs Finished / ' num2str(nSNP) ' SNPs analyzed']);

    RS = RS(:,FinishedJOBs);RS = RS(:);RS = RS(index);
    A1 = A1(:,FinishedJOBs);A1 = A1(:);A1 = A1(index);
    A2 = A2(:,FinishedJOBs);A2 = A2(:);A2 = A2(index);
    CHR = CHR(:,FinishedJOBs);CHR = CHR(:);CHR = CHR(index);
    MAF = MAF(:,FinishedJOBs);MAF = MAF(:);MAF = MAF(index);
    N = N(:,FinishedJOBs);N = N(:);N = N(index);
    JOBIID = JOBIID(:,FinishedJOBs);JOBIID = JOBIID(:);JOBIID = JOBIID(index);
   
    P = P(:,:,FinishedJOBs);
    P = permute(P,[1 3 2]);
    P = reshape(P,size(P,1)*size(P,2),size(P,3));
    P = P(index,:);
    P = single(P)/factor;
%% 
%SNP = permute(SNP,[1 3 2]);
%SNP = reshape(SNP,size(SNP,1)*size(SNP,2),size(SNP,3)); 

    CHI = CHI(:,:,FinishedJOBs);
    CHI = permute(CHI,[1 3 2]);
    CHI = reshape(CHI,size(CHI,1)*size(CHI,2),size(CHI,3));
    CHI = CHI(index,:);
    
    SIGN = SIGN(:,:,FinishedJOBs);
    SIGN = permute(SIGN,[1 3 2]);
    SIGN = reshape(SIGN,size(SIGN,1)*size(SIGN,2),size(SIGN,3));
    SIGN = SIGN(index,:);

    GWAS.CHR = CHR;
    GWAS.POS = POS;
    GWAS.RS = RS;
    GWAS.A1 = A1;
    GWAS.A2 = A2;
    GWAS.MAF = MAF;
    GWAS.N = N;
    GWAS.P = P;
    GWAS.CHI = CHI;
    GWAS.SIGN = SIGN;
    %GWAS.SNP = SNP;
    GWAS.JOBID = JOBIID;
    clear CHR POS RS P A1 A2 MAF N CHI SIGN;
end