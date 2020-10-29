function out = runCCAJOB(JOB,DB,PA,saving)
         
         if nargin<4, saving = 'all'; end % how much of the data is saved into out
         warning off;
         %[ind12,ind21] = int_vlookupFast(DB.IID,JOB.IID);% match subjects having both genotypes and phenotypes    
         [ind12,ind21] = vlookupFast(DB.IID,JOB.IID);% match subjects having both genotypes and phenotypes    
         nSNP = size(JOB.SNP,2);
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
         A = nan*zeros(nSNP,nSegments);
         F = nan*zeros(nSNP,nSegments);
         CHI = nan*zeros(nSNP,nSegments);
         pCHI = nan*zeros(nSNP,nSegments);
         df1 = nan*zeros(nSNP,nSegments);
         df2 = nan*zeros(nSNP,nSegments);
         B = cell(nSNP,nSegments);
         status = ones(1,nSNP,'uint8');  
         for s=1:1:nSNP% for each SNP
             %s=1;
             %disp(num2str(s));
             try
                % s=1;
                snp = double(GT(:,s));snp(snp==-1) = nan;
                index = find(~isnan(snp));% select snp and missing values
                %snp = int_getResiduals(COV(index,:),snp(index));% correct ifo covariates
                snp = getResiduals(COV(index,:),snp(index));% correct ifo covariates
                forpval = nan*zeros(1,nSegments);
                forcor = nan*zeros(1,nSegments);
                forA = nan*zeros(1,nSegments);
                forF = nan*zeros(1,nSegments);
                forCHI = nan*zeros(1,nSegments);
                forpCHI = nan*zeros(1,nSegments);
                fordf1 = nan*zeros(1,nSegments);
                fordf2 = nan*zeros(1,nSegments);
                forB = cell(1,nSegments);
                parfor i=1:nSegments
                    % i=1;
                    [forA(i),forB{i},forcor(i),~,~,STATS] = canoncorr(snp,PA{i}(index,:));                    
                    forB{i} = sign(forA(i))*forB{i};% correct for a possible flip;
                    forpval(i) = STATS.pF(end);
                    fordf1(i) = STATS.df1;
                    fordf2(i) = STATS.df2;
                    forCHI(i) = STATS.chisq(end);
                    forpCHI(i) = STATS.pChisq(end);
                    forF(i) = STATS.F(end);
                end
                pval(s,:) = forpval;
                cor(s,:) = forcor;
                A(s,:) = forA;
                F(s,:) = forF;
                CHI(s,:) = forCHI;
                pCHI(s,:) = forpCHI;
                df1(s,:) = fordf1;
                df2(s,:) = fordf2; 
                B(s,:) = forB;
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
                 out.A = A;
                 out.F = single(F);
                 out.CHI2 = single(CHI);
                 out.pCHI2 = pCHI;
                 out.df1 = df1;
                 out.df2 = df2;
                 out.B = B;
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


                    
%                     [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(snp,PA{i}(index,:),1);
%                     [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(snp,PA{i}(index,:),1);
%                     
%                     [out,outd] = angle(YL,forB{i})
%                     
%                     [out,outd] = angle(BETA(2,:)',forB{i})
%                     
%                     [XL,YL,XS,YS,BETA,PCTVAR,MSE] = plsregress(PA{i}(index,:),snp,1);
%                     
%                      [out,outd] = angle(XL,forB{i})