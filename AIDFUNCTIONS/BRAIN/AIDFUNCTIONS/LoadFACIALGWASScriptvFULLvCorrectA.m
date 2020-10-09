%% LOADING FACE GWAS
% LOAD the results over izquirda y derecha hemisphere, and plot Miami style
%close all;clear all;
function GWAS = LoadFACIALGWASScriptvFULLvCorrectA
    studypath = '/home/pclaes4/Documents/MATLAB/FACEGWAS/';
    if isfolder(studypath)
        disp('LOADING LOCAL VERSION on HOME');
        INPUTPATH = [studypath 'SNPJOBSINTERSECTpoint8/'];
        OUTPUTPATH = [studypath 'OUTPUTINTERSECT_FULLSTATS/'];
    else
        disp('LOADING AVALOK VERSION');
        studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2018/METAEUROGWAS/DATA/';
        INPUTPATH = [studypath 'GENOTYPES/SNPJOBSINTERSECTpoint8/'];
        OUTPUTPATH = [studypath 'GENOTYPES/OUTPUTINTERSECT_FULLSTATS/'];
    end
    % change back to avalok
    studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2018/METAEUROGWAS/DATA/';
    PHENOPATH = [studypath 'PHENOTYPES/'];
    in = load([PHENOPATH 'DATAB_0727']);
    DB = cell(1,2);
    nDB = 2;
    DB{1}.nSegments = 63;
    DB{2}.nSegments = 63;
    DB{1}.IID = in.DATABASE{1}.IID;DB{1}.nSubj = length(DB{1}.IID);
    DB{2}.IID = in.DATABASE{2}.IID;DB{2}.nSubj = length(DB{2}.IID);
    clear in;
    cd(studypath);
    %% OPEN WORKERS IF YOU DON'T HAVE ANY YET
    try
        parpool('LocalSingle',30);
    catch
    end
    %% LOAD FACIAL GWAS US ONLY
    db = 1;
    %% PREPARE DATA CONTAINERS
    cd(INPUTPATH);files = dir('**/*.mat');% retrieve all job folders recursively
    nJOBS = length(files);
    nSNPperJOB = 2000;
    RS = cell(nSNPperJOB,nJOBS);
    POS = zeros(nSNPperJOB,nJOBS,'uint32');
    A1 = cell(nSNPperJOB,nDB,nJOBS);
    A2 = cell(nSNPperJOB,nDB,nJOBS);
    MAF = zeros(nSNPperJOB,nJOBS,'uint8');
    HWE = zeros(nSNPperJOB,nJOBS,'uint32');
    CHR = zeros(nSNPperJOB,nJOBS,'uint8');
    CHI = zeros(nSNPperJOB,DB{db}.nSegments,nDB,nJOBS,'single');
    SIGN = zeros(nSNPperJOB,DB{db}.nSegments,nDB,nJOBS,'int8');
    P = zeros(nSNPperJOB,DB{db}.nSegments,nDB,nDB+1,nJOBS,'uint32');
    N = zeros(nSNPperJOB,2,nJOBS,'uint16');
    factor = 10000;
    %% HOUSEHOLDING, making sure that the JOBS will be loaded in order
    CHRID = zeros(1,nJOBS);
    JOBID = zeros(1,nJOBS);
    for i=1:nJOBS
        % i=1;
        ind = strfind(files(i).folder,'CHR');
        CHRID(i) = str2double(files(i).folder(ind+3:end));
        ind = strfind(files(i).name,'.');
        JOBID(i) = str2double(files(i).name(4:ind-1));
    end
    [~,index] = sort(CHRID,'ascend');% order in chromosome
    files = files(index);
    CHRID = zeros(1,nJOBS);
    JOBID = zeros(1,nJOBS);
    for i=1:nJOBS
        % i=1;
        ind = strfind(files(i).folder,'CHR');
        CHRID(i) = str2double(files(i).folder(ind+3:end));
        ind = strfind(files(i).name,'.');
        JOBID(i) = str2double(files(i).name(4:ind-1));
    end
    for i=1:1:max(CHRID)% order within each chromosome in ascending JOB order
       ind = find(CHRID==i);
       [~,index] = sort(JOBID(ind),'ascend');
       files(ind) = files(ind(index));
    end
    CHRID = zeros(1,nJOBS);
    JOBID = zeros(1,nJOBS);
    for i=1:nJOBS
        % i=1;
        ind = strfind(files(i).folder,'CHR');
        CHRID(i) = str2double(files(i).folder(ind+3:end));
        ind = strfind(files(i).name,'.');
        JOBID(i) = str2double(files(i).name(4:ind-1));
    end
    %% LOAD RESULTS
    [path,ID] = setupParForProgress(nJOBS);tic;
    parfor i=1:nJOBS
        % i=1;
        inputpath = [INPUTPATH 'CHR' num2str(CHRID(i)) '/'];
        outputpath = [OUTPUTPATH 'CHR' num2str(CHRID(i)) '/'];
        inputfile = [inputpath 'JOB' num2str(JOBID(i)) '.mat'];
        outputfile = dir([outputpath num2str(JOBID(i)) '_*']);
        if isempty(outputfile), continue; end
        outputfile = [outputpath outputfile(1).name];
        in = load(inputfile);
        % COUNTING INDIVIDUALS DATASET 1
        JOB = in.JOB{1};
        [ind12,~] = vlookupFast(DB{1}.IID,JOB.IID);% match subjects having both genotypes and phenotypes
        JOB.SNP = JOB.SNP(ind12,:);% extract SNP data from overlapping subjects
        N1 = sum(JOB.SNP>=0,1);
        A11 = JOB.A1;
        A12 = JOB.A2;
        % COUNTING INDIVIDUALS DATASET 2
        JOB = in.JOB{2};
        [ind12,~] = vlookupFast(DB{2}.IID,JOB.IID);% match subjects having both genotypes and phenotypes
        JOB.SNP = JOB.SNP(ind12,:);% extract SNP data from overlapping subjects
        N2 = sum(JOB.SNP>=0,1);
        A21 = JOB.A1;
        A22 = JOB.A2;
        
        JOB = in.JOB{db};
        JOB.RSID = in.JOB{2}.RSID;% TAKE UK AS REFERENCE, since it contains proper RS numbers
        [ind12,~] = vlookupFast(DB{db}.IID,JOB.IID);% match subjects having both genotypes and phenotypes
        JOB.SNP = JOB.SNP(ind12,:);% extract SNP data from overlapping subjects
        nSNP = length(JOB.POS);
        maf = zeros(nSNP,1);
        a1 = cell(nSNP,nDB);
        a2 = cell(nSNP,nDB);
        hwe = zeros(nSNP,1);
        for s=1:1:nSNP
           % s=1;
           geno = JOB.SNP(:,s);  
           [geno,flip] = codetoM(geno);
%            if flip
%                a1{s} = JOB.A2{s};
%                a2{s} = JOB.A1{s};
%            else
%                a1{s} = JOB.A1{s};
%                a2{s} = JOB.A2{s};
%            end
           a1{s,1} = A11{s};
           a2{s,1} = A12{s};
           a1{s,2} = A21{s};
           a2{s,2} = A22{s};
           maf(s) = getMAF(geno);
           hwe(s) = -log10(myHWE(geno));
           JOB.SNP(:,s) = geno;
        end   
        if nSNP==nSNPperJOB
           RS(:,i) = JOB.RSID(:);
           A1(:,:,i) = a1;
           A2(:,:,i) = a2;
           POS(:,i) = uint32(JOB.POS(:));
           CHR(:,i) = JOB.CHRID(:);
           N(:,:,i) = uint16([N1(:), N2(:)]);
           MAF(:,i) = uint8(100*maf(:));
           HWE(:,i) = uint32(factor*hwe(:));
        else
           tmpRS = cell(nSNPperJOB,1);tmpRS(1:nSNP) = JOB.RSID;RS(:,i) = tmpRS(:);
           tmpA1 = cell(nSNPperJOB,nDB,1);tmpA1(1:nSNP,:) = a1;A1(:,:,i) = tmpA1;
           tmpA2 = cell(nSNPperJOB,nDB,1);tmpA2(1:nSNP,:) = a2;A2(:,:,i) = tmpA2;
           tmpPOS = zeros(nSNPperJOB,1);tmpPOS(1:nSNP) = JOB.POS;POS(:,i) = uint32(tmpPOS(:));
           tmpCHR = zeros(nSNPperJOB,1);tmpCHR(1:nSNP) = JOB.CHRID;CHR(:,i) = uint8(tmpCHR(:));
           tmpMAF = zeros(nSNPperJOB,1);tmpMAF(1:nSNP) = uint8(100*maf(:));MAF(:,i) = tmpMAF;
           tmpHWE = zeros(nSNPperJOB,1);tmpHWE(1:nSNP) = uint32(factor*hwe(:));HWE(:,i) = tmpHWE;
           tmpN = zeros(nSNPperJOB,2);tmpN(1:nSNP,:) = uint16([N1(:), N2(:)]);N(:,:,i) = tmpN;
        end
        in = load(outputfile);
        in = in.in;
        p = in.pvalues;p(:,:,:,3) = in.pfastvalues;
        chi = squeeze(in.FULLSTATS(:,1,:,:));
        effectsign = int8(sign(squeeze(in.FULLSTATS(:,8,:,:))));
        if nSNP==nSNPperJOB
           CHI(:,:,:,i) = single(chi);
           SIGN(:,:,:,i) = effectsign;
           P(:,:,:,:,i) = uint32(factor*-log10(p));
        else
           tmpCHI = zeros(nSNPperJOB,DB{db}.nSegments,nDB,'single');tmpCHI(1:nSNP,:,:) = single(chi);CHI(:,:,:,i)=tmpCHI;
           tmpSIGN = zeros(nSNPperJOB,DB{db}.nSegments,nDB,'int8');tmpSIGN(1:nSNP,:,:) = effectsign;SIGN(:,:,:,i)=tmpSIGN;
           tmpP = zeros(nSNPperJOB,DB{db}.nSegments,nDB,nDB+1,'uint32');tmpP(1:nSNP,:,:,:) = uint32(factor*-log10(p));P(:,:,:,:,i)=tmpP; 
        end
        parfor_progress;
    end
    closeParForProgress(path,ID);toc;
    %% REDUCING DATA
    TEST = sum(POS,1);
    FinishedJOBs = find(TEST);
    nFinishedJOBs = length(FinishedJOBs);% when GWAS is fully finished this is not an issue , and all jobs will be finished
    POS = POS(:,FinishedJOBs);POS = POS(:);index = find(POS);POS = POS(index);% take zero positions, since these are empty data containers
    nSNP = length(POS);
    disp([num2str(nFinishedJOBs) ' Jobs Finished / ' num2str(nSNP) ' SNPs analyzed']);

    RS = RS(:,FinishedJOBs);RS = RS(:);RS = RS(index);
    CHR = CHR(:,FinishedJOBs);CHR = CHR(:);CHR = CHR(index);
    MAF = MAF(:,FinishedJOBs);MAF = MAF(:);MAF = MAF(index);
    HWE = HWE(:,FinishedJOBs);HWE = HWE(:);HWE = HWE(index);

    N = N(:,:,FinishedJOBs);
    N = permute(N,[1 3 2]);
    N = reshape(N,size(N,1)*size(N,2),size(N,3));
    N = N(index,:);
    
    A1 = A1(:,:,FinishedJOBs);
    A1 = permute(A1,[1 3 2]);
    A1 = reshape(A1,size(A1,1)*size(A1,2),size(A1,3));
    A1 = A1(index,:);
    
    A2 = A2(:,:,FinishedJOBs);
    A2 = permute(A2,[1 3 2]);
    A2 = reshape(A2,size(A2,1)*size(A2,2),size(A2,3));
    A2 = A2(index,:);

    CHI = CHI(:,:,:,FinishedJOBs);
    CHI = permute(CHI,[1 4 2 3]);
    CHI = reshape(CHI,size(CHI,1)*size(CHI,2),size(CHI,3),size(CHI,4));
    CHI = CHI(index,:,:);
    
    SIGN = SIGN(:,:,:,FinishedJOBs);
    SIGN = permute(SIGN,[1 4 2 3]);
    SIGN = reshape(SIGN,size(SIGN,1)*size(SIGN,2),size(SIGN,3),size(SIGN,4));
    SIGN = SIGN(index,:,:);

    P = P(:,:,:,:,FinishedJOBs);
    P = permute(P,[1 5 2 3 4]);
    P = reshape(P,size(P,1)*size(P,2),size(P,3),size(P,4),size(P,5));
    P = P(index,:,:,:);

    GWAS.CHR = CHR;
    GWAS.POS = POS;
    GWAS.RS = RS;
    GWAS.A1 = A1;
    GWAS.A2 = A2;
    GWAS.MAF = MAF;
    GWAS.HWE = HWE;
    GWAS.CHI = CHI;
    GWAS.SIGN = SIGN;
    GWAS.P = P;
    GWAS.N = N;
    clear CHR POS RS A1 A2 MAF HWE SNP CHI P N;
end