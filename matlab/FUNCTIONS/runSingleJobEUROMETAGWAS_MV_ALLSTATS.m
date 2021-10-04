function out = runSingleJobEUROMETAGWAS_MV_ALLSTATS(JOB,DB,PA,HI,type,saving,precision,progress,PAtransformation)
         if nargin<5, type = 'cca'; end % type is a forward CCA combined with RIP scoring, or reverse regression
         if nargin<6, saving = 'all'; end % how much of the data is saved into out
         if nargin<7, precision = 'double'; end% what is the precision of the resulting input and output data
         if nargin<8, progress = false; end
         if nargin<9, PAtransformation = 'none';end
         out.Type = type;
         warning off;
         nDB = length(DB);
         nSNP = length(JOB{1}.POS);
         status = ones(1,nSNP,'uint8');  
         % retrieving and preparing DATA to USE for each database
         DATA = cell(1,nDB);
         for i=1:1:nDB
             [ind12,ind21] = vlookupFast(DB{i}.IID,JOB{i}.IID);% match subjects having both genotypes and phenotypes
             DATA{i}.PA = reduceCellData(HI,PA,DB{i}.PAIND(ind21));% extract shape features from overlapping subjects
             NOTUSE = strcmp(JOB{i}.RSID,'NA');% FLAG SNPS THAT ARE PRESENT
             SNP = double(JOB{i}.SNP(ind12,:));% extract SNP data from overlapping subjects
             SNP(SNP==-1) = nan;
             COV = DB{i}.COV(ind21,:);% extract covariates from overlapping subjects
             DATA{i}.ORIGSNP = SNP;
             for s=1:1:nSNP % Correction of SNP variant ifo covariates
                 if NOTUSE(s), continue; end
                 snp = SNP(:,s);index = find(~isnan(snp));% select snp and missing values
                 if isempty(index),NOTUSE(s) = 1; continue; end
                 try
                    new = getResiduals(COV(index,:),snp(index));% correct ifo covariates
                    %new = snp(index);
                 catch % something else went wrong
                    NOTUSE(s) = 1;continue;
                 end
                 SNP(index,s) = new;
             end
             DATA{i}.SNP = SNP;
             DATA{i}.NOTUSE = NOTUSE;
             DATA{i}.A1 = JOB{i}.A1;
             DATA{i}.A2 = JOB{i}.A2;
         end
         %clear PA DB JOB COV SNP;% free memory
         % RUN ALL ANALYSIS
         switch lower(type)
             case 'cca'
                 pvalues = nan*ones(nSNP,HI.nLC,nDB,nDB,precision);% precision initiatilization
                 pfastvalues = nan*zeros(nSNP,HI.nLC,nDB,precision);% precision initiatilization
                 cor = nan*zeros(nSNP,HI.nLC,nDB,precision);% precision initiatilization
                 FULLSTATS = nan*zeros(nSNP,8,HI.nLC,nDB,precision);% precision initiatilization
                 B = cell(nSNP,HI.nLC,nDB);
                 if progress, [path,ID] = setupParForProgress(nSNP);end
                 parfor s=1:1:nSNP% for each SNP
                     warning off;
                     forpvalues = nan*ones(HI.nLC,nDB,nDB,precision);% precision initiatilization
                     forpfastvalues = nan*zeros(HI.nLC,nDB,precision);% precision initiatilization
                     forcor = nan*zeros(HI.nLC,nDB,precision);% precision initiatilization
                     forFULLSTATS = nan*zeros(8,HI.nLC,nDB,precision);% precision initiatilization
                     forB = cell(HI.nLC,nDB);
                     % s=3756;
                     %disp(num2str(s));
                     try
                         for dd=1:1:nDB % for each dataset as discovery
                             %dd=1
                             if DATA{dd}.NOTUSE(s), continue; end %#ok<*PFBNS>
                             pa = DATA{dd}.PA;% select shape features
                             X = DATA{dd}.SNP(:,s);% select snp data
                             index = find(~isnan(X));% select non missing SNP data
                             pa = reduceCellData(HI,pa,index);X = X(index);% select overlap 
                             da1 = DATA{dd}.A1{s};da2 = DATA{dd}.A2{s};
                             for m=1:1:HI.nLC% for each facial module
                                 % m=1;
                                 %Y = getCellData(HI,pa,m);% select facial module
                                 Y = tranformY(getCellData(HI,pa,m),PAtransformation);% select facial module
                                 [A,forB{m,dd},forcor(m,dd),~,~,STATS] = canoncorr(X,Y);% run cannocical correlation analysis, keep direction and pF
                                 forB{m,dd} = sign(A)*forB{m,dd};% correct for a possible flip; 
                                 forpvalues(m,dd,dd) = STATS.pF(end);% enter on diagonal in matrix
                                 forFULLSTATS(1,m,dd) = STATS.chisq(end);
                                 forFULLSTATS(2,m,dd) = STATS.pChisq(end);
                                 forFULLSTATS(3,m,dd) = STATS.F(end);
                                 forFULLSTATS(4,m,dd) = STATS.pF(end);
                                 forFULLSTATS(5,m,dd) = STATS.df1(end);
                                 forFULLSTATS(6,m,dd) = STATS.df2(end);
                                 forFULLSTATS(7,m,dd) = STATS.Wilks(end);
                                 forFULLSTATS(8,m,dd) = A;
                             end % module loop
                             for rd=1:1:nDB % for each other dataset as replication
                                 % rd= 2
                                 if dd==rd, continue;end% same dataset as discovery
                                 if DATA{rd}.NOTUSE(s), continue; end % SNP not present in particular dataset
                                 ra1 = DATA{rd}.A1{s};ra2 = DATA{rd}.A2{s};
                                 flip = 2;
                                 if strcmp(ra1,da1)||strcmp(ra2,da2), flip = 0; end
                                 if strcmp(ra2,da1)||strcmp(ra1,da2), flip = 1; end
                                 if flip==2, continue; end% Do not execute, because the Alleles do not match
                                 pa = DATA{rd}.PA;% select shape features
                                 X = DATA{rd}.SNP(:,s);% select snp data
                                 index = find(~isnan(X));% select non missing SNP data
                                 pa = reduceCellData(HI,pa,index);X = X(index);% select overlap 
                                 for m=1:1:HI.nLC% for each facial module
                                     % m=13;
                                     %Y = getCellData(HI,pa,m);% select facial module
                                     Y = tranformY(getCellData(HI,pa,m),PAtransformation);% select facial module
                                     %RIP = PPMMV3.getRIP(Y,B{m}');% don't forget to remove the PPMMV3. adn the transpose on B
                                     RIP = getRIP(Y,forB{m,dd});
                                     if flip==1, RIP = -1*RIP; end
                                     try
                                         STATS = regstats(RIP',X,'linear','tstat');% perform univariate regression on RIP value
                                         onesided = tcdf(-1*STATS.tstat.t(end),STATS.tstat.dfe);% return the one-sided uppertail p-value instead of the two sided
                                         forpvalues(m,dd,rd) = onesided;% enter on the off-diagonal elements in matrix
                                         %forpvalues(m,dd,rd) = STATS.tstat.pval(end);% enter on the off-diagonal elements in matrix
                                     catch % something went wrong
                                     end
                                 end % module loop
                             end % Replication db loop
                         end % Discovery db loop
                         % Getting meta pvalues
                         for db=1:1:nDB
                             %db=1;
                             for m=1:1:HI.nLC
                                 %m=1;
                                 P = squeeze(forpvalues(m,db,:));
                                 index = find(~isnan(P));
                                 if isempty(index), continue;end
                                 %forpfastvalues(m,db) = pfast(P(index));
                                 forpfastvalues(m,db) = stouffer(P(index));
                             end
                         end
                     catch
                         status(s) = uint8(2);
                     end
                     pvalues(s,:,:,:) = forpvalues;
                     pfastvalues(s,:,:) = forpfastvalues;% precision initiatilization
                     cor(s,:,:) = forcor;% precision initiatilization
                     FULLSTATS(s,:,:,:) = forFULLSTATS;
                     B(s,:,:) = forB;
                     if progress, parfor_progress;end
                 end % SNP loop
                 if progress, closeParForProgress(path,ID);end
                 % linking to output
                 out.status = status;
                 out.pvalues = pvalues;
                 out.pfastvalues= pfastvalues;
                 switch lower(saving)
                     case 'all'
                         out.COR = uint16(cor*10000);
                         %reduced storage
%                          for e=1:1:numel(B)
%                             B{e} = single(B{e}); 
%                          end
%                          out.B = B;
                         out.FULLSTATS = FULLSTATS;
                     case 'reduced'
                         return;
                 end
             case 'revreg'
                 % META RESULTS
                 META_B = cell(nSNP,HI.nLC);
                 META_Qe = nan*zeros(nSNP,HI.nLC,precision);% precision initiatilization
                 META_pQe = nan*zeros(nSNP,HI.nLC,precision);% precision initiatilization
                 META_Qb = nan*zeros(nSNP,HI.nLC,precision);% precision initiatilization
                 META_pQb = nan*zeros(nSNP,HI.nLC,precision);% precision initiatilization
                 META_mse = nan*zeros(nSNP,HI.nLC,precision);% precision initiatilization
                 pvalues = nan*zeros(nSNP,HI.nLC,nDB,precision);% precision initiatilization
                 pfastvalues = nan*zeros(nSNP,HI.nLC,precision);% precision initiatilization
                 
                 DBpresent = zeros(nSNP,nDB);
                 
                 if progress, [path,ID] = setupParForProgress(nSNP);end
                 parfor s=1:1:nSNP% for each SNP
                     %s=1;
                     warning off;
                     B = cell(HI.nLC,nDB);
                     covB = cell(HI.nLC,nDB);
                     forpvalues = nan*zeros(HI.nLC,nDB,precision);% precision initiatilization
                     forpfastvalues = nan*zeros(1,HI.nLC,precision);% precision initiatilization
                     beta = cell(1,HI.nLC);
                     Qe = nan*zeros(1,HI.nLC,precision);% precision initiatilization
                     pQe = nan*zeros(1,HI.nLC,precision);% precision initiatilization
                     Qb = nan*zeros(1,HI.nLC,precision);% precision initiatilization
                     pQb = nan*zeros(1,HI.nLC,precision);% precision initiatilization
                     mse = nan*zeros(1,HI.nLC,precision);% precision initiatilization
                     present = ones(1,nDB);
                     % s=1;
                     %disp(num2str(s));
                     try
                         % select reference allele
                         refdd = 0;
                         for dd=1:1:nDB
                             if ~(DATA{dd}.NOTUSE(s)), refdd = dd; break; end
                         end
                         refa1 = DATA{refdd}.A1{s};refa2 = DATA{refdd}.A2{s};
                         for dd=1:1:nDB % for each dataset as discovery
                             %dd=1
                             if DATA{dd}.NOTUSE(s), present(dd) = 0; continue;  end %#ok<*PFBNS>
                             dda1 = DATA{dd}.A1{s};dda2 = DATA{dd}.A2{s};
                             flip = 2;
                             if strcmp(dda1,refa1)||strcmp(dda2,refa2), flip = 0; end
                             if strcmp(dda2,refa1)||strcmp(dda1,refa2), flip = 1; end
                             if flip==2, present(dd) = 0; continue; end% I STILL CANNOT INCORPORATE PROPERLY                           
                             pa = DATA{dd}.PA;% select shape features
                             X = DATA{dd}.SNP(:,s);% select snp data
                             index = find(~isnan(X));% select non missing SNP data
                             pa = reduceCellData(HI,pa,index);X = X(index);% select overlap 
                             for m=1:1:HI.nLC% for each facial module
                                 % m=1;
                                 %Y = getCellData(HI,pa,m);% select facial module
                                 Y = tranformY(getCellData(HI,pa,m),PAtransformation);% select facial module
                                 revreg = regstats(X,Y,'linear',{'beta' 'covb' 'fstat'});% run reverse regression (independent = shape, dependent SNP)                                          
                                 forpvalues(m,dd) = revreg.fstat.pval;
                                 covB{m,dd} = revreg.covb;
                                 B{m,dd} = revreg.beta;
                             end % module loop
                         end % database loop   
                         % Getting pfast meta pvalues
                         for m=1:1:HI.nLC% META ANALYSIS OVER EACH MODULE
                             %m=1;
                             index = find(present);
                             if sum(index)>1
%                                 pa = DATA{index(1)}.PA;% select shape features
%                                 Y = getCellData(HI,pa,m);% select facial module 
%                                 nPred = size(Y,2)+1;
                                [beta{m},Qe(m),pQe(m),Qb(m),pQb(m),mse(m)] = metaReverseRegression(B(m,index),covB(m,index));
                             elseif sum(index)==1
                                 beta{m} = B(m,index);
                                 pQb = forpvalues(m,index);
                             else
                                 %do nothing;
                             end
                             P = squeeze(forpvalues(m,:)');
                             index = find(~isnan(P));
                             if isempty(index), continue;end
                             %forpfastvalues(m) = pfast(P(index));
                             forpfastvalues(m) = stouffer(P(index));
                         end % module loop
                         % Getting regression of slopes meta pvalues
                     catch
                         status(s) = uint8(2);
                     end       
                     META_B(s,:) = beta;
                     META_Qe(s,:) = Qe;
                     META_pQe(s,:) = pQe;
                     META_Qb(s,:) = Qb;
                     META_pQb(s,:) = pQb;
                     META_mse(s,:) = mse;
                     pvalues(s,:,:) = forpvalues;
                     pfastvalues(s,:) = forpfastvalues;
                     DBpresent(s,:) = present;
                     if progress, parfor_progress;end
                 end % SNP loop
                 if progress, closeParForProgress(path,ID);end
                 % linking to output
                 out.status = status;
                 out.META_pQb = META_pQb;
                 out.META_pQe = META_pQe;
                 switch lower(saving)
                     case 'all'
                         out.META_B = META_B;
                         out.META_Qb = META_Qb;
                         %out.META_pQb = META_pQb;
                         out.META_Qe = META_Qe;
                         %out.META_pQe = META_pQe;
                         out.META_mse = META_mse;
                         out.pvalues = pvalues;
                         out.pfastvalues = pfastvalues;
                         out.DBpresent = DBpresent;
                     case 'reduced'
                         return;
                 end
         end
         warning on;
end

function out = getResiduals(X,Y)
                [~,~,~,~,M] = plsregress(X,Y,min(size(X,2),size(Y,2)));
                Y_est = [ones(size(X,1),1) X]*M;
                out = Y-Y_est;
end

function [out] = getRIP(in,M)
            out = dot(in',repmat(M/norm(M),1,size(in,1)));
end

function out = tranformY(Y,transformation)
         switch lower(transformation)
             case 'none'
                 out = Y;
             case 'abs'% probably usefull for asymmetry
                 out = abs(Y);
             case 'square'
                 out = Y.^2;
             otherwise
                 out = Y;
         end
end

% function out = runSingleJobEUROMETAGWAS_dev(JOB,DB,PA,HI,type,saving,precision)
%          if nargin<5, type = 'cca'; end % type is a forward CCA combined with RIP scoring, or reverse regression
%          if nargin<6, saving = 'all'; end % how much of the data is saved into out
%          if nargin<7, precision = 'double'; end% what is the precision of the resulting input and output data
%          out.Type = type;
%          warning off;
%          nDB = length(DB);
%          nSNP = length(JOB{1}.POS);
%          status = ones(1,nSNP,'uint8');  
%          % retrieving and preparing DATA to USE for each database
%          DATA = cell(1,nDB);
%          for i=1:1:nDB
%              [ind12,ind21] = vlookupFast(DB{i}.IID,JOB{i}.IID);% mach subjects having both genotypes and phenotypes
%              DATA{i}.PA = reduceCellData(HI,PA,DB{i}.PAIND(ind21));% extract shape features from overlapping subjects
%              NOTUSE = strcmp(JOB{i}.RSID,'NA');% FLAG SNPS THAT ARE PRESENT
%              SNP = double(JOB{i}.SNP(ind12,:));% extract SNP data from overlapping subjects
%              SNP(SNP==-1) = nan;
%              COV = DB{i}.COV(ind21,:);% extract covariates from overlapping subjects
%              DATA{i}.ORIGSNP = SNP;
%              for s=1:1:nSNP % Correction of SNP variant ifo covariates
%                  if NOTUSE(s), continue; end
%                  snp = SNP(:,s);index = find(~isnan(snp));% select snp and missing values
%                  if isempty(index),NOTUSE(s) = 1; continue; end
%                  try
%                     new = getResiduals(COV(index,:),snp(index));% correct ifo covariates
%                  catch % something else went wrong
%                     NOTUSE(s) = 1;continue;
%                  end
%                  SNP(index,s) = new;
%              end
%              DATA{i}.SNP = SNP;
%              DATA{i}.NOTUSE = NOTUSE;
%              DATA{i}.A1 = JOB{i}.A1;
%              DATA{i}.A2 = JOB{i}.A2;
%          end
%          clear PA DB JOB COV SNP;% free memory
%          % RUN ALL ANALYSIS
%          
%          switch lower(type)
%              case 'cca'
%                  pvalues = nan*ones(nSNP,HI.nLC,nDB,nDB,precision);% precision initiatilization
%                  pfastvalues = nan*zeros(nSNP,HI.nLC,nDB,precision);% precision initiatilization
%                  cor = nan*zeros(nSNP,HI.nLC,nDB,precision);% precision initiatilization
%                  B = cell(nSNP,HI.nLC,nDB);
%                  parfor s=1:1:nSNP% for each SNP
%                      forpvalues = nan*ones(HI.nLC,nDB,nDB,precision);% precision initiatilization
%                      forpfastvalues = nan*zeros(HI.nLC,nDB,precision);% precision initiatilization
%                      forcor = nan*zeros(HI.nLC,nDB,precision);% precision initiatilization
%                      forB = cell(HI.nLC,nDB);
%                      % s=1;
%                      disp(num2str(s));
%                      try
%                          for dd=1:1:nDB % for each dataset as discovery
%                              %dd=1
%                              if DATA{dd}.NOTUSE(s), continue; end
%                              pa = DATA{dd}.PA;% select shape features
%                              X = DATA{dd}.SNP(:,s);% select snp data
%                              index = find(~isnan(X));% select non missing SNP data
%                              pa = reduceCellData(HI,pa,index);X = X(index);% select overlap 
%                              da1 = DATA{dd}.A1{s};da2 = DATA{dd}.A2{s};
%                              for m=1:1:HI.nLC% for each facial module
%                                  % m=1;
%                                  Y = getCellData(HI,pa,m);% select facial module
%                                  [~,forB{s,m,dd},forcor(s,m,dd),~,~,STATS] = canoncorr(X,Y);% run cannocical correlation analysis, keep direction and pF
%                                  pvalues(s,m,dd,dd) = STATS.pF(end);% enter on diagonal in matrix
%                              end % module loop
%                              for rd=1:1:nDB % for each other dataset as replication
%                                  % rd= 3
%                                  if dd==rd, continue;end% same dataset as discovery
%                                  if DATA{rd}.NOTUSE(s), continue; end % SNP not present in particular dataset
%                                  ra1 = DATA{rd}.A1{s};ra2 = DATA{rd}.A2{s};
%                                  flip = 2;
%                                  if strcmp(ra1,da1)||strcmp(ra2,da2), flip = 0; end
%                                  if strcmp(ra2,da1)||strcmp(ra1,da2), flip = 1; end
%                                  if flip==2, continue; end% Do not execute, because the Alleles do not match
%                                  pa = DATA{rd}.PA;% select shape features
%                                  X = DATA{rd}.SNP(:,s);% select snp data
%                                  index = find(~isnan(X));% select non missing SNP data
%                                  pa = reduceCellData(HI,pa,index);X = X(index);% select overlap 
%                                  for m=1:1:HI.nLC% for each facial module
%                                      % m=1;
%                                      Y = getCellData(HI,pa,m);% select facial module
%                                      %RIP = PPMMV3.getRIP(Y,B{m}');% don't forget to remove the PPMMV3. adn the transpose on B
%                                      RIP = getRIP(Y,B{s,m,dd});
%                                      if flip==1, RIP = -1*RIP; end
%                                      try
%                                          STATS = regstats(RIP',X,'linear','tstat');% perform univariate regression on RIP value
%                                          pvalues(s,m,dd,rd) = STATS.tstat.pval(end);% enter on the off-diagonal elements in matrix
%                                      catch % something went wrong
%                                      end
%                                  end % module loop
%                              end % Replication db loop
%                          end % Discovery db loop
%                          % Getting meta pvalues
%                          for db=1:1:nDB
%                              %db=1;
%                              for m=1:1:HI.nLC
%                                  %m=1;
%                                  P = squeeze(pvalues(s,m,db,:));
%                                  index = find(~isnan(P));
%                                  if isempty(index), continue;end
%                                  pfastvalues(s,m,db) = pfast(P(index));
%                              end
%                          end
%                      catch
%                          status(s) = uint8(2);
%                      end
%                  end % SNP loop
%              case 'reverse regression'
%                  %                                  % reverse regression
% %                                  revreg = regstats(X,Y,'linear',{'beta' 'covb' 'mse' 'fstat'});
% 
% 
%                  
%                  
%          end
%          % linking to output
%          out.status = status;
%          out.pvalues = pvalues;
%          out.pfastvalues= pfastvalues;
%          warning on;
% end
% 
% function out = getResiduals(X,Y)
%                 [~,~,~,~,M] = plsregress(X,Y,min(size(X,2),size(Y,2)));
%                 Y_est = [ones(size(X,1),1) X]*M;
%                 out = Y-Y_est;
% end
% 
% function [out] = getRIP(in,M)
%             out = dot(in',repmat(M/norm(M),1,size(in,1)));
% end