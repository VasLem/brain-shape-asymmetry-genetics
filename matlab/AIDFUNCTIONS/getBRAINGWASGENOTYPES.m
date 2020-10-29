function OUT = getBRAINGWASGENOTYPES(SEL,IID)
    %% LOADING BRAIN GWAS
    INPUTPATH = '/IMAGEN/BRAIN/UKBIOBANK/GENOTYPES/SNPJOBSv500/';
    nSNPperJOB = 500;
    nSubjects = 19644;
    jobIDS = unique(SEL.JOBID);
    jobIDS = setdiff(jobIDS,0);
    nJOBS = length(jobIDS);
    %% OPEN WORKERS IF YOU DON'T HAVE ANY YET
    try
        parpool('LocalSingle',30);
    catch
    end
    %% PREPING DATA CONTAINERS
    RS = cell(nSNPperJOB,nJOBS);
    CHR = zeros(nSNPperJOB,nJOBS,'uint8');
    POS = zeros(nSNPperJOB,nJOBS,'uint32');
    SNP = zeros(nSNPperJOB,nSubjects,nJOBS,'int8');
    %% LOAD RESULTS
    [path,ID] = setupParForProgress(nJOBS);tic;
    parfor i=1:nJOBS
        %i=1;
        jobfile = dir([INPUTPATH 'JOB_' num2str(jobIDS(i)) '_*']);
        in = load([INPUTPATH jobfile.name]);
        JOB = in.JOB;
        nSNP = length(in.JOB.POS);
        [ind12,~] = vlookupFast(IID,JOB.IID);%
        JOB.SNP = JOB.SNP(ind12,:);
        if nSNP==nSNPperJOB
           SNP(:,:,i) = JOB.SNP';
           POS(:,i) = uint32(in.JOB.POS(:));
           CHR(:,i) = in.JOB.CHRID(:);
           RS(:,i) = in.JOB.RSID(:);
        else
           tmpRS = cell(nSNPperJOB,1);tmpRS(1:nSNP) = in.JOB.RSID;RS(:,i) = tmpRS(:);
           tmpPOS = zeros(nSNPperJOB,1);tmpPOS(1:nSNP) = in.JOB.POS;POS(:,i) = uint32(tmpPOS(:));
           tmpCHR = zeros(nSNPperJOB,1);tmpCHR(1:nSNP) = in.JOB.CHRID;CHR(:,i) = uint8(tmpCHR(:)); 
           tmpSNP = zeros(nSNPperJOB,nSubjects,'int8');tmpSNP(1:nSNP,:) = JOB.SNP';SNP(:,:,i) = tmpSNP;
        end
        parfor_progress;
    end
    closeParForProgress(path,ID);toc;
    %% REDUCING DATA
    disp('RESHAPING DATA');
    TEST = sum(POS,1);
    FinishedJOBs = find(TEST);
    nFinishedJOBs = length(FinishedJOBs);
    POS = POS(:,FinishedJOBs);POS = POS(:);index = find(POS);POS = POS(index);
    nSNP = length(POS);
    disp([num2str(nFinishedJOBs) ' Jobs LOADED / ' num2str(nSNP) ' SNPs']);
    RS = RS(:,FinishedJOBs);RS = RS(:);RS = RS(index);
    CHR = CHR(:,FinishedJOBs);CHR = CHR(:);CHR = CHR(index);
    SNP = SNP(:,:,FinishedJOBs);
    SNP = permute(SNP,[1 3 2]);
    SNP = reshape(SNP,size(SNP,1)*size(SNP,2),size(SNP,3));  
    %% MATCH BACK TO SELECTION AND GENERATE OUTPUT
    disp('MATCHING');
    [ind12,ind21] =  vlookupFast([uint32(SEL.CHR) SEL.POS],[uint32(CHR) POS],true);
    %[ind12,ind21] =  vlookupFast(SEL.RS,RS);
    disp(['Found ' num2str(length(ind12)) ' out of ' num2str(length(SEL.POS)) ' SNPs']);
    OUT = -1*ones(length(SEL.POS),size(SNP,2),'int8');
    OUT(ind21,:) = SNP(ind12,:);
    clear RS CHR POS SNP;
    disp('DONE');
%% THE END
end



