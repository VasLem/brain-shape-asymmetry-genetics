function [out] = TRTrueMatchingWeightsFULL(objSNP,objBase)
          % INTITIALIZING
           snps = objSNP.SNPs;
           nSNP = length(snps);
           faces = objBase.DepVar;
           n = size(faces,1);
           tGT = getValues(objSNP,'GT');
           tCOV = objBase.Cov;
           tGB = objBase.GB;
           nCOV = objBase.nrCov;
           nGB = objBase.nrGB;
           BF = 0;
           tGRIP = nan*zeros(n,nSNP);
           for i=1:1:nSNP
               rip = SNP.updateRIP(faces,snps{i}.ST4MOD.M);
               rip = SNP.normalizeRIP(rip,snps{i}.ST4MOD.M,snps{i}.ST4MOD.Var);
               tGRIP(:,i) = rip;
           end
           % BASE
           tBRIP = updateRIP(faces,objBase.ST4MOD.M)';
           tBRIP = BaseContainer.normalizeRIP(tBRIP,objBase.ST4MOD.M,objBase.ST4MOD.Var);
         % MATCHES
           TSMatches = nan*zeros(nCOV,n);
           TGBMatches = nan*zeros(nGB,n);
           TGMatches = nan*zeros(nSNP,n);
           FSMatches = nan*zeros(nCOV,n);
           FGBMatches = nan*zeros(nGB,n);
           FGMatches = nan*zeros(nSNP,n);
           for t=1:1:n
               %disp(num2str(t));
               [~,covD] = CovX2RIP(objBase,tCOV(t,:),3,300);
               [~,gbD] = GBX2RIP(objBase,tGB(t,:),3,300);
               [~,tmp] = matchCovRIP(objBase,tBRIP(t,1:nCOV),covD);
               TSMatches(:,t) = tmp;
               [~,tmp] = inversematchCovRIP(objBase,tBRIP(t,1:nCOV),covD);
               FSMatches(:,t) = tmp;
               [~,TGBMatches(:,t)] = matchGBRIP(objBase,tBRIP(t,nCOV+1:end),gbD);
               [~,FGBMatches(:,t)] = inversematchGBRIP(objBase,tBRIP(t,nCOV+1:end),gbD);
               Ttmp = nan*zeros(nSNP,1);
               Ftmp = nan*zeros(nSNP,1);
               for s=1:1:nSNP
                    Ttmp(s) = matchRIP(snps{s},tGRIP(t,s)',tGT(t,s),BF);
                    Ftmp(s) = inversematchRIP(snps{s},tGRIP(t,s)',tGT(t,s),BF);
               end
               TGMatches(:,t) = Ttmp;
               FGMatches(:,t) = Ftmp;
           end
           TSMatches = -log(TSMatches);TGBMatches = -log(TGBMatches);TGMatches = -log(TGMatches);
           TMatches = [TSMatches;TGBMatches;TGMatches];
           AvgT = nanmean(TMatches,2);out.AvgT = AvgT;
           StdT = nanstd(TMatches,0,2);out.StdT = StdT;
           WAvgT = AvgT-min(AvgT);
           WAvgT = WAvgT/max(WAvgT);
           WAvgT = 1-WAvgT;out.WAvgT = WAvgT;
           %MedT = nanmedian(TMatches,2);out.MedT = MedT;
           %MadT = nanmad(TMatches,0,2);out.MadT = MadT;
           %WMedT = MedT-min(MedT);
           %WMedT = WMedT/max(WMedT);
           %WMedT = 1-WMedT;out.WMedT = WMedT;
           FSMatches = -log(FSMatches);FGBMatches = -log(FGBMatches);FGMatches = -log(FGMatches);
           FMatches = [FSMatches;FGBMatches;FGMatches];
           AvgF = nanmean(FMatches,2);out.AvgF = AvgF;
           StdF = nanstd(FMatches,0,2);out.StdF = StdF;
           D = abs(AvgT-AvgF);out.D = D;
           WD = D-min(D);
           WD = WD/max(WD);out.WD = WD;
           S = sqrt(((n-1).*StdT.^2+(n-1).*StdF.^2)./(n+n-2));out.S = S;
           CD = D./S;out.CD = CD;
           WCD = CD-min(CD);
           WCD = WCD/max(WCD);out.WCD = WCD;       
end
           
%            [bestW,avgW] = GAWvTF(TtMatches,rMatches,wM,GAoptions);
%            [bestW,avgW] = GAWvTF(TtMatches,oMatches,wM,GAoptions);
%            [bestW,avgW] = GAWvTF(TtMatches,FtMatches,wM,GAoptions);
%            
