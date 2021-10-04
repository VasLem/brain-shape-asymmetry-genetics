function out = runGWASJOB_VSC_2020(JobFile,AuxFile,OutFile)
    %maxNumCompThreads(1);
    out = 0;
    try
    %% LOADING DATA
        disp('LOADING DATA');
        in = load(AuxFile);
        PA = in.PA;% Phenotype data
        DB = in.DB;% Covariate data
        clear in;
        in = load(JobFile);% genotype data
        JOB = in.JOB;
        clear in;
        CHRID = JOB.CHRID(1);
        disp('LOADING DONE');
    %% RUNNING JOB    
       disp(['STARTING: ' JobFile]);
       res = runJOB(JOB,DB,PA,'all');
       res.JOBNR = JobFile;res.CHRID = CHRID; %#ok<*PFBNS>
       disp(['DONE: ' JobFile]);
    %% SAVING
       disp('SAVING RESULTS');
       save(OutFile,'res');
       disp('SAVING DONE');
    catch
        disp('ERROR DIFFERENT FROM FALSE SNP');
        out = 1;
    end
    %% DONE
end

function out = runJOB(JOB,DB,PA,saving)
         
         if nargin<4, saving = 'all'; end % how much of the data is saved into out
         warning off;
         [ind12,ind21] = int_vlookupFast(DB.IID,JOB.IID);% match subjects having both genotypes and phenotypes    
         nSNP = length(JOB.POS);
         GT = JOB.SNP(ind12,:);
         %GT(GT==-1) = nan;
         COV = DB.COV(ind21,:);
         %JOB = rmfield(JOB,'SNP');
         DepVar = PA.DepVar;
         nSegments = PA.nSegments;
         PA = cell(1,nSegments);
         for i=1:1:nSegments
            PA{i} = double(DepVar{i}(ind21,:)); 
         end
         pval = nan*zeros(nSNP,nSegments);
         cor = nan*zeros(nSNP,nSegments);
         F = nan*zeros(nSNP,nSegments);
         CHI = nan*zeros(nSNP,nSegments);
         pCHI = nan*zeros(nSNP,nSegments);
         df1 = nan*zeros(nSNP,nSegments);
         df2 = nan*zeros(nSNP,nSegments);
         status = ones(1,nSNP,'uint8');  
         for s=1:1:nSNP% for each SNP
             disp(num2str(s));
             try
                % s=1;
                snp = double(GT(:,s));snp(snp==-1) = nan;
                index = find(~isnan(snp));% select snp and missing values
                snp = int_getResiduals(COV(index,:),snp(index));% correct ifo covariates
                forpval = nan*zeros(1,nSegments);
                forcor = nan*zeros(1,nSegments);
                forF = nan*zeros(1,nSegments);
                forCHI = nan*zeros(1,nSegments);
                forpCHI = nan*zeros(1,nSegments);
                fordf1 = nan*zeros(1,nSegments);
                fordf2 = nan*zeros(1,nSegments);
                for i=1:1:nSegments
                    % i=1;
                    [~,~,forcor(i),~,~,STATS] = canoncorr(snp,PA{i}(index,:));
                    forpval(i) = STATS.pF(end);
                    fordf1(i) = STATS.df1;
                    fordf2(i) = STATS.df2;
                    forCHI(i) = STATS.chisq(end);
                    forpCHI(i) = STATS.pChisq(end);
                    forF(i) = STATS.F(end);
                end
                pval(s,:) = forpval;
                cor(s,:) = forcor;
                F(s,:) = forF;
                CHI(s,:) = forCHI;
                pCHI(s,:) = forpCHI;
                df1(s,:) = fordf1;
                df2(s,:) = fordf2;  
             catch
                 disp(['SNP ' num2str(s) ' ERROR']);
                 status(s) = 2;
             end
         end
         warning on;
         out.pvalues = pval;
         out.status = status;
         switch saving
             case 'all'
                 out.COR = uint16(cor*10000);
                 out.F = single(F);
                 out.CHI2 = single(CHI);
                 out.pCHI2 = pCHI;
                 out.df1 = df1;
                 out.df2 = df2;
             otherwise
                 return;
         end
end

function out = int_getResiduals(X,Y)
                [~,~,~,~,M] = plsregress(X,Y,min(size(X,2),size(Y,2)));
                Y_est = [ones(size(X,1),1) X]*M;
                out = Y-Y_est;
end

function [ind12,ind21,raw] = int_vlookupFast(Array1,Array2,rows)
         if nargin<3, rows = false; end
         if rows
            [~,ind12] = ismember(Array1,Array2,'rows');
         else
            [~,ind12] = ismember(Array1,Array2);
         end
         raw = ind12;
         ind21 = find(ind12);
         ind12 = ind12(ind21);
end