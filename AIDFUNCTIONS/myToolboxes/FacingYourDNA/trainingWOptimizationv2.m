function out = trainingWOptimizationv2(objSNP,objBase,testGT,testS,testGB,freq,nT,nR,GAoptions,SNPOnly)
           if nargin<10, SNPOnly = false; end;
         % INTITIALIZING
           snps = objSNP.SNPs;
           nSNP = length(snps);
           faces = objBase.DepVar;
           n = size(faces,1);
           GT = getValues(objSNP,'GT');
           COV = objBase.Cov;
           GB = objBase.GB;
           nCOV = objBase.nrCov;
           nGB = objBase.nrGB;
           BF = 0;
           nW = 1+nGB+nSNP; 
           DGT = GTDistance(testGT,GT,freq)';
           if ~isempty(testS)&&~SNPOnly
              DS = 1-(repmat(testS,n,1)==COV(:,1))';
           else
              DS = zeros(1,n);
           end
           if ~isempty(testGB)&&~SNPOnly
              DGB = sqrt(sum((repmat(testGB,n,1)-GB).^2,2))';
           else
              DGB = zeros(1,n);
           end
           D = DGT+DS+DGB;
           [~,sortind] = sort(D);
         % TEST SET  
           Tind = sortind(1:nT);
           tfaces = faces(Tind,:);
           tCOV = COV(Tind,:);
           tGB = GB(Tind,:);
           tGT = GT(Tind,:);
         % Opposite SET
           Oind = sortind(end-nT+1:end);
           ofaces = faces(Oind,:);
         % Random SET %nR
           rfaces = randomFaces(objBase,nR);
         % RIPS
           % SNPS
           tGRIP = nan*zeros(nT,nSNP);
           oGRIP = nan*zeros(nT,nSNP);
           rGRIP = nan*zeros(nR,nSNP);
           for i=1:1:nSNP
               rip = SNP.updateRIP(tfaces,snps{i}.ST4MOD.M);
               rip = SNP.normalizeRIP(rip,snps{i}.ST4MOD.M,snps{i}.ST4MOD.Var);
               tGRIP(:,i) = rip;
               rip = SNP.updateRIP(rfaces,snps{i}.ST4MOD.M);
               rip = SNP.normalizeRIP(rip,snps{i}.ST4MOD.M,snps{i}.ST4MOD.Var);
               rGRIP(:,i) = rip;
               rip = SNP.updateRIP(ofaces,snps{i}.ST4MOD.M);
               rip = SNP.normalizeRIP(rip,snps{i}.ST4MOD.M,snps{i}.ST4MOD.Var);
               oGRIP(:,i) = rip;
           end
           % BASE
           tBRIP = updateRIP(tfaces,objBase.ST4MOD.M)';
           tBRIP = BaseContainer.normalizeRIP(tBRIP,objBase.ST4MOD.M,objBase.ST4MOD.Var);
           oBRIP = updateRIP(ofaces,objBase.ST4MOD.M)';
           oBRIP = BaseContainer.normalizeRIP(oBRIP,objBase.ST4MOD.M,objBase.ST4MOD.Var);
           rBRIP = updateRIP(rfaces,objBase.ST4MOD.M)';
           rBRIP = BaseContainer.normalizeRIP(rBRIP,objBase.ST4MOD.M,objBase.ST4MOD.Var);
         % MATCHES
           tSMatches = nan*zeros(1,nT,nT);
           tGBMatches = nan*zeros(nGB,nT,nT);
           tGMatches = nan*zeros(nSNP,nT,nT);
           
           oSMatches = nan*zeros(1,nT,nT);
           oGBMatches = nan*zeros(nGB,nT,nT);
           oGMatches = nan*zeros(nSNP,nT,nT);
           
           rSMatches = nan*zeros(1,nT,nR);
           rGBMatches = nan*zeros(nGB,nT,nR);
           rGMatches = nan*zeros(nSNP,nT,nR);
           
           wG = nan*zeros(nSNP,nT);
           wGB = objBase.ST4MOD.RIPStat(5:end)';
           wS = 1;
           for t=1:1:nT
               [~,covD] = CovX2RIP(objBase,tCOV(t,:),3,300);
               [~,gbD] = GBX2RIP(objBase,tGB(t,:),2,100);
               tmp = matchCovRIP(objBase,tBRIP(:,1:nCOV),covD);
               tSMatches(1,t,:) = tmp(1,:)';
               tGBMatches(:,t,:) = matchGBRIP(objBase,tBRIP(:,nCOV+1:end),gbD);
               tmp = matchCovRIP(objBase,oBRIP(:,1:nCOV),covD);
               oSMatches(1,t,:) = tmp(1,:)';
               oGBMatches(:,t,:) = matchGBRIP(objBase,tBRIP(:,nCOV+1:end),gbD);
               tmp = matchCovRIP(objBase,rBRIP(:,1:nCOV),covD);
               rSMatches(1,t,:) = tmp(1,:);%rCOVMatches(:,:,t) = tmp;
               rGBMatches(:,t,:) = matchGBRIP(objBase,rBRIP(:,nCOV+1:end),gbD);
               tmptG = nan*zeros(nSNP,nT);
               tmpoG = nan*zeros(nSNP,nT);
               tmprG = nan*zeros(nSNP,nR);
               for s=1:1:nSNP
                    tmptG(s,:) = matchRIP(snps{s},tGRIP(:,s)',tGT(t,s),BF);
                    tmpoG(s,:) = matchRIP(snps{s},oGRIP(:,s)',tGT(t,s),BF);
                    tmprG(s,:) = matchRIP(snps{s},rGRIP(:,s)',GT(t,s),BF);
                    wG(s,t) = GTTypicality(snps{s},GT(t,s));
               end
               tGMatches(:,t,:) = tmptG;
               oGMatches(:,t,:) = tmpoG;
               rGMatches(:,t,:) = tmprG;
           end
           tSMatches = -log(tSMatches);tGBMatches = -log(tGBMatches);tGMatches = -log(tGMatches);
           oSMatches = -log(tSMatches);oGBMatches = -log(oGBMatches);oGMatches = -log(oGMatches);
           rSMatches = -log(rSMatches);rGBMatches = -log(rGBMatches);rGMatches = -log(rGMatches);  
           if SNPOnly
              tMatches = tGMatches;
              oMatches = oGMatches;
              rMatches = rGMatches;  
              wM = wG;
           else
              tMatches = [tSMatches;tGBMatches;tGMatches];
              oMatches = [oSMatches;oGBMatches;oGMatches];
              rMatches = [rSMatches;rGBMatches;rGMatches];
              wM = [wS*ones(1,nT);repmat(wGB,1,nT);wG];
           end
           [TtMatches,FtMatches] = seperateMatches(tMatches);
           FcMatches = permute([permute(oMatches,[1 3 2]), permute(FtMatches,[1 3 2])],[1 3 2]);
           out.NR = zeros(1,4);out.NEER = zeros(1,4);
           [out.NR(1),out.NEER(1)] = optWScoreAndEvalTF(TtMatches,rMatches,wM,ones(size(wM,1),1));
           [out.NR(2),out.NEER(2)] = optWScoreAndEvalTF(TtMatches,oMatches,wM,ones(size(wM,1),1));
           [out.NR(3),out.NEER(3)] = optWScoreAndEvalTF(TtMatches,FtMatches,wM,ones(size(wM,1),1));
           [out.NR(4),out.NEER(4)] = optWScoreAndEvalTF(TtMatches,FcMatches,wM,ones(size(wM,1),1));
          % optimize
           [bestW,avgW] = GAWvTF(TtMatches,FcMatches,wM,GAoptions);
           out.WR = zeros(1,4);out.WEER = zeros(1,4);
           [out.WR(1),out.WEER(1)] = optWScoreAndEvalTF(TtMatches,rMatches,wM,avgW);
           [out.WR(2),out.WEER(2)] = optWScoreAndEvalTF(TtMatches,oMatches,wM,avgW);
           [out.WR(3),out.WEER(3)] = optWScoreAndEvalTF(TtMatches,FtMatches,wM,avgW);
           [out.WR(4),out.WEER(4)] = optWScoreAndEvalTF(TtMatches,FcMatches,wM,avgW);
           out.AvgW = avgW;
           out.BestW = bestW;  
end
           
%            [bestW,avgW] = GAWvTF(TtMatches,rMatches,wM,GAoptions);
%            [bestW,avgW] = GAWvTF(TtMatches,oMatches,wM,GAoptions);
%            [bestW,avgW] = GAWvTF(TtMatches,FtMatches,wM,GAoptions);
%            
