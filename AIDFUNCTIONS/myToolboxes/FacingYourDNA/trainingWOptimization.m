function out = trainingWOptimization(objSNP,objBase,GAoptions,nT)
         % INTITIALIZING
           snps = objSNP.SNPs;
           nSNP = length(snps);
           faces = objBase.DepVar;
           GT = getValues(objSNP,'GT');
           [n,~] = size(faces);
           COV = objBase.Cov;
           GB = objBase.GB;
           nCOV = objBase.nrCov;
           nGB = objBase.nrGB;
           BF = 0;
           nW = 1+nGB+nSNP;
         % TEST SET  
           Tind = randsample(n,nT,false);
           faces = faces(Tind,:);
           COV = COV(Tind,:);
           GB = GB(Tind,:);
           GT = GT(Tind,:);
           n = nT;
           input = eye(n,n);
           g = ind2logical(input(:));
         % RIPS
           % SNPS
           GRIP = nan*zeros(n,nSNP);
           for i=1:1:nSNP
               rip = SNP.updateRIP(faces,snps{i}.ST4MOD.M);
               rip = SNP.normalizeRIP(rip,snps{i}.ST4MOD.M,snps{i}.ST4MOD.Var);
               GRIP(:,i) = rip;
           end
           % BASE
           BRIP = updateRIP(faces,objBase.ST4MOD.M)';
           BRIP = BaseContainer.normalizeRIP(BRIP,objBase.ST4MOD.M,objBase.ST4MOD.Var);
           %SRIP = BRIP(:,1);GBRIP = BRIP(:,nCOV+1:end);
         % MATCHES
           SMatches = nan*zeros(1,n,n);
           GBMatches = nan*zeros(nGB,n,n);
           GMatches = nan*zeros(nSNP,n,n);
           wG = nan*zeros(nSNP,n,n);
           wGB = repmat(objBase.ST4MOD.RIPStat(5:end)',1,n,n);
           wS = 1*ones(1,n,n);
           parfor t=1:1:n
               [~,covD] = CovX2RIP(objBase,COV(t,:),3,300);
               [~,gbD] = GBX2RIP(objBase,GB(t,:),2,100);
               tmp = matchCovRIP(objBase,BRIP(:,1:nCOV),covD);
               SMatches(1,:,t) = tmp(1,:)';
               GBMatches(:,:,t) = matchGBRIP(objBase,BRIP(:,nCOV+1:end),gbD);
               tmpG = nan*zeros(nSNP,n);
               tmpwG = nan*zeros(nSNP,n);
               for s=1:1:nSNP
                    tmpG(s,:) = matchRIP(snps{s},GRIP(:,s)',GT(t,s),BF);
                    tmpwG(s,:) = GTTypicality(snps{s},GT(t,s))*ones(1,n);
               end
               GMatches(:,:,t) = tmpG;
               wG(:,:,t) = tmpwG;
           end
           %SMatches = wS.*-log(SMatches);
           %GBMatches = wGB.*-log(GBMatches);
           %GMatches = wG.*-log(GMatches);
           SMatches = -log(SMatches);
           GBMatches = -log(GBMatches);
           GMatches = -log(GMatches);
           Matches = [SMatches;GBMatches;GMatches];
           wM = [wS;wGB;wG];
         % UNWEIGHTED
           [out.R,out.EER] = optWScoreAndEval(Matches,wM,ones(nW,1),n,g);
         % WEIGTH OPTIMIZATION
           [out.BestW,out.AvgW] = GAW(Matches,wM,n,g,GAoptions);
           [out.BestWR,out.BestWEER] = optWScoreAndEval(Matches,wM,out.BestW,n,g);
           [out.AvgWR,out.AvgWEER] = optWScoreAndEval(Matches,wM,out.AvgW,n,g);
end