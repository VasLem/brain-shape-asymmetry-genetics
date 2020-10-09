function GWAS = LoadBRAINGWASmichp03(type)
    if nargin<1, type='info'; end
    switch lower(type)
        case 'info'
            full = false;
        case 'gwas'
            full = true;
    end
    disp(full)
    %% LOADING BRAIN GWAS
    disp('LOADING 2020 michp03 version');
    OUTPUTPATH = '/IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/GWASOUTPUT2020/';
    INPUTPATH = '/IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/SNPJOBSv500/';
    HS = [];
    HS.OUTPUTPATH = [OUTPUTPATH '/SWH_PA80_20200513/'];
    in = load('/IMAGEN/BRAIN/UKBIOBANK/PHENOTYPES/SWH/AUX_SWH_PA80_20200513.mat');
    HS.IID = in.DB.IID;
    HS.nSubjects = length(HS.IID);
    HS.Segments = in.PA.Segments;
    HS.nSegments = in.PA.nSegments;
    clear in;
    nSNPperJOB = 500;
    %% OPEN WORKERS IF YOU DON'T HAVE ANY YET
    try
        parpool('LocalSingle',30);
    catch
    end
    %% PREPARE DATA CONTAINERS
    files = dir([INPUTPATH '*.mat']);
    nJOBS = length(files);clear files;
    RS = cell(nSNPperJOB,nJOBS);
    CHR = zeros(nSNPperJOB,nJOBS,'uint8');
    POS = zeros(nSNPperJOB,nJOBS,'uint32');
    JOBID = zeros(nSNPperJOB,nJOBS,'uint16');
    A1 = cell(nSNPperJOB,nJOBS);
    A2 = cell(nSNPperJOB,nJOBS);
    MAF = zeros(nSNPperJOB,nJOBS,'uint8');
    if full
        N = zeros(nSNPperJOB,2,nJOBS,'uint16');
        PH = zeros(nSNPperJOB,HS.nSegments,nJOBS,'uint32');
        CHIH = zeros(nSNPperJOB,HS.nSegments,nJOBS,'single');
    end
    factor = 10000;
    %% LOAD RESULTS
    [path,ID] = setupParForProgress(nJOBS);tic;
    parfor i=1:nJOBS
        % i=1
        jobfile = dir([INPUTPATH 'JOB_' num2str(i) '_*']);
        if isempty(jobfile), continue; end
        in = load([INPUTPATH jobfile.name]);
        JOB = in.JOB;

        nSNP = length(in.JOB.POS);
        [ind12,ind21] = vlookupFast(HS.IID,JOB.IID);%
        JOB.SNP = JOB.SNP(ind12,:);% extract SNP data from overlapping subjects
        n = sum(JOB.SNP>=0,1);
        n = n(:);
        maf = zeros(nSNP,1);
        a1 = cell(nSNP,1);
        a2 = cell(nSNP,1);
        for s=1:1:nSNP
           % s=1;
           geno = JOB.SNP(:,s);  
           [geno,flip] = codetoM(geno);
           a1{s} = JOB.A1{s};
           a2{s} = JOB.A2{s};
           maf(s) = getMAF(geno);
           JOB.SNP(:,s) = geno;
        end
        JOBID(:,i) = uint16(i);
        if nSNP==nSNPperJOB
           RS(:,i) = in.JOB.RSID(:);
           A1(:,i) = a1(:);
           A2(:,i) = a2(:);
           POS(:,i) = uint32(in.JOB.POS(:));
           CHR(:,i) = in.JOB.CHRID(:);
           if full, N(:,i) = uint16(n);end
           MAF(:,i) = uint8(100*maf(:));
        else
           tmpRS = cell(nSNPperJOB,1);tmpRS(1:nSNP) = in.JOB.RSID;RS(:,i) = tmpRS(:);
           tmpA1 = cell(nSNPperJOB,1);tmpA1(1:nSNP) = a1;A1(:,i) = tmpA1(:);
           tmpA2 = cell(nSNPperJOB,1);tmpA2(1:nSNP) = a2;A2(:,i) = tmpA2(:);
           tmpPOS = zeros(nSNPperJOB,1);tmpPOS(1:nSNP) = in.JOB.POS;POS(:,i) = uint32(tmpPOS(:));
           tmpCHR = zeros(nSNPperJOB,1);tmpCHR(1:nSNP) = in.JOB.CHRID;CHR(:,i) = uint8(tmpCHR(:));
           tmpMAF = zeros(nSNPperJOB,1);tmpMAF(1:nSNP) = uint8(100*maf(:));MAF(:,i) = tmpMAF;
           if full, tmpN = zeros(nSNPperJOB,1);tmpN(1:nSNP) = uint16(n);N(:,i) = tmpN;end
        end
        if full
            in = load([HS.OUTPUTPATH jobfile.name]);
            if nSNP==nSNPperJOB
               pvalues = -log10(in.res.pvalues);
               CHI = in.res.CHI2;
            else
               pvalues = nan*zeros(nSNPperJOB,size(in.res.pvalues,2));
               pvalues(1:nSNP,:) = -log10(in.res.pvalues);
               CHI = nan*zeros(nSNPperJOB,size(in.res.CHI2,2));
               CHI(1:nSNP,:) = in.res.CHI2;
            end
            PH(:,:,i) = uint32(pvalues*factor);
            CHIH(:,:,i) = CHI;
        end
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
    if full,N = N(:,FinishedJOBs);N = N(:);N = N(index);end
    JOBID = JOBID(:,FinishedJOBs);JOBID = JOBID(:);JOBID = JOBID(index);

    if full
        PH = PH(:,:,FinishedJOBs);
        PH = permute(PH,[1 3 2]);
        PH = reshape(PH,size(PH,1)*size(PH,2),size(PH,3));
        PH = PH(index,:);

        CHIH = CHIH(:,:,FinishedJOBs);
        CHIH = permute(CHIH,[1 3 2]);
        CHIH = reshape(CHIH,size(CHIH,1)*size(CHIH,2),size(CHIH,3));
        CHIH = CHIH(index,:);
    end

    GWAS.CHR = CHR;
    GWAS.POS = POS;
    GWAS.RS = RS;
    GWAS.JOBID = JOBID;
    GWAS.A1 = A1;
    GWAS.A2 = A2;
    GWAS.MAF = MAF;
    clear CHR POS RS A1 A2 MAF;
    if full 
        GWAS.N = N;
        GWAS.P = PH;
        GWAS.CHI = CHIH;
        clear PH N CHIH;
    end
%% THE END
end



