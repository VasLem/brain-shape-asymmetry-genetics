function GWAS = LoadBRAINGWASScript
    %% LOADING BRAIN GWAS
    studypath = '/uz/data/avalok/mic/tmp/pclaes4/==MATLAB==/ActiveProjects/2020/UKB/DATA/';
    disp('LOADING 2020 LOCAL version');
    %OUTPUTPATH = '/uz/data/avalok/mic/tmp/SHARED/pclaes4/BRAIN/UKBIOBANK/GENOTYPES/GWASOUTPUT2020/';
    OUTPUTPATH = '/data/pclaes4/UKB/GWASOUTPUT2020/';
    %disp('LOADING pre2020 version');
    %OUTPUTPATH = '/uz/data/avalok/mic/tmp/SHARED/pclaes4/BRAIN/UKBIOBANK/GENOTYPES/GWASOUTPUT/';
    cd(studypath);
    GENOPATH = [studypath 'GENOTYPING/'];cd(GENOPATH);
    %INPUTPATH = '/uz/data/avalok/mic/tmp/SHARED/pclaes4/BRAIN/UKBIOBANK/GENOTYPES/SNPJOBSv500/';
    INPUTPATH = '/data/pclaes4/UKB/SNPJOBSv500/';
    nHS = 2;
    HS = cell(1,nHS);
    HS{1}.OUTPUTPATH = [OUTPUTPATH '/LH_PA80_20190928/'];
    in = load('/uz/data/avalok/mic/tmp/SHARED/pclaes4/BRAIN/UKBIOBANK/PHENOTYPES/LH/AUX_LH_PA80_20190928.mat');
    HS{1}.IID = in.DB.IID;
    HS{1}.nSubjects = length(HS{1}.IID);
    HS{1}.Segments = in.PA.Segments;
    HS{1}.nSegments = in.PA.nSegments;
    clear in;
    HS{2}.OUTPUTPATH = [OUTPUTPATH '/RH_PA80_20190928/'];
    in = load('/uz/data/avalok/mic/tmp/SHARED/pclaes4/BRAIN/UKBIOBANK/PHENOTYPES/RH/AUX_RH_PA80_20190928.mat');
    HS{2}.IID = in.DB.IID;
    HS{2}.nSubjects = length(HS{2}.IID);
    HS{2}.Segments = in.PA.Segments;
    HS{2}.nSegments = in.PA.nSegments;
    clear in;
    nSNPperJOB = 500;
    nSubjects = 19644;
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
    N = zeros(nSNPperJOB,2,nJOBS,'uint16');
    PLH = zeros(nSNPperJOB,HS{1}.nSegments,nJOBS,'uint32');
    PRH = zeros(nSNPperJOB,HS{2}.nSegments,nJOBS,'uint32');%
    %SNP = zeros(nSNPperJOB,nSubjects,nJOBS,'int8');
    factor = 10000;
    %% LOAD RESULTS
    [path,ID] = setupParForProgress(nJOBS);tic;
    parfor i=1:nJOBS
        % i=1
        jobfile = dir([INPUTPATH 'JOB_' num2str(i) '_*']);
        if isempty(jobfile), continue; end
        present = ones(1,2);
        files = cell(1,2);
        for j=1:1:2
            files{j} = dir([HS{j}.OUTPUTPATH jobfile.name]);
            if isempty(files{j}), present(j) = 0; end
        end
        if sum(present)<2, continue;end% present as output in both
        in = load([INPUTPATH jobfile.name]);
        JOB = in.JOB;

        nSNP = length(in.JOB.POS);
        [ind12,ind21] = vlookupFast(HS{1}.IID,JOB.IID);%
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
%            if flip
%                a1{s} = JOB.A2{s};
%                a2{s} = JOB.A1{s};
%            else
%                a1{s} = JOB.A1{s};
%                a2{s} = JOB.A2{s};
%            end
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
           N(:,i) = uint16(n);
           MAF(:,i) = uint8(100*maf(:));
           %SNP(:,:,i) = JOB.SNP';
        else
           tmpRS = cell(nSNPperJOB,1);tmpRS(1:nSNP) = in.JOB.RSID;RS(:,i) = tmpRS(:);
           tmpA1 = cell(nSNPperJOB,1);tmpA1(1:nSNP) = a1;A1(:,i) = tmpA1(:);
           tmpA2 = cell(nSNPperJOB,1);tmpA2(1:nSNP) = a2;A2(:,i) = tmpA2(:);
           tmpPOS = zeros(nSNPperJOB,1);tmpPOS(1:nSNP) = in.JOB.POS;POS(:,i) = uint32(tmpPOS(:));
           tmpCHR = zeros(nSNPperJOB,1);tmpCHR(1:nSNP) = in.JOB.CHRID;CHR(:,i) = uint8(tmpCHR(:));
           tmpMAF = zeros(nSNPperJOB,1);tmpMAF(1:nSNP) = uint8(100*maf(:));MAF(:,i) = tmpMAF;
           tmpN = zeros(nSNPperJOB,1);tmpN(1:nSNP) = uint16(n);N(:,i) = tmpN;
           %tmpSNP = zers(nSNPperJOB,nSubjects,'int8');tmpSNP(1:nSNP,:) = JOB.SNP';SNP(:,:,i) = tmpSNP;
        end
        for j=1:2
            % j=1;
            in = load([HS{j}.OUTPUTPATH jobfile.name]);
            if nSNP==nSNPperJOB
               pvalues = -log10(in.res.pvalues);
            else
               pvalues = nan*zeros(nSNPperJOB,size(in.res.pvalues,2));
               pvalues(1:nSNP,:) = -log10(in.res.pvalues);
            end
            switch j
                case 1
                    PLH(:,:,i) = uint32(pvalues*factor);
                case 2
                    PRH(:,:,i) = uint32(pvalues*factor);
            end
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
    N = N(:,FinishedJOBs);N = N(:);N = N(index);
    JOBID = JOBID(:,FinishedJOBs);JOBID = JOBID(:);JOBID = JOBID(index);

    PLH = PLH(:,:,FinishedJOBs);
    PLH = permute(PLH,[1 3 2]);
    PLH = reshape(PLH,size(PLH,1)*size(PLH,2),size(PLH,3));
    PLH = PLH(index,:);
    %PLH = single(PLH)/factor;

    PRH = PRH(:,:,FinishedJOBs);
    PRH = permute(PRH,[1 3 2]);
    PRH = reshape(PRH,size(PRH,1)*size(PRH,2),size(PRH,3));
    PRH = PRH(index,:);
    %PRH = single(PRH)/factor;

    GWAS.CHR = CHR;
    GWAS.POS = POS;
    GWAS.RS = RS;
    GWAS.JOBID = JOBID;
    GWAS.A1 = A1;
    GWAS.A2 = A2;
    GWAS.MAF = MAF;
    GWAS.N = N;
    GWAS.PLH = PLH;
    GWAS.PRH = PRH;
    clear CHR POS RS PLH PRH A1 A2 MAF N;
%% THE END
end



