classdef BIOMETRICDNAMOD < superClassSNP
    properties
        SNPs = {};
        Base = {};
    end
    properties (Dependent = true)
        nSNP;
        nUniqSNP;
        nBase;
    end
    properties
        Match = 'QD';
        Classify = 'SOFT';
        W = 'P';
        HTest = true;
        RIPNorm = true;
    end
    methods
        function obj = BIOMETRICDNAMOD(varargin)
            obj = obj@superClassSNP(varargin{:});         
        end
    end
    methods % GENERAL GETTING SETTING
        function out = get.nSNP(obj)
            out = length(obj.SNPs);
        end
        function out = get.nBase(obj)
            out = length(obj.Base);
        end
        function out = get.nUniqSNP(obj)
           if obj.nSNP==0, out = 0; return; end
           rs = getRS(obj);
           out = length(unique(rs));
        end  
    end
    methods % GENERAL INTERFACING
        function out = getRS(obj,index)
           if obj.nSNP==0, out = []; return;end
           if nargin < 2, index = 1:obj.nSNP; end
           out = cell(1,length(index));
           for i=1:1:length(index)
              if isempty(obj.SNPs{index(i)}), continue; end 
              if ~strcmp(obj.SNPs{index(i)}.Type,'USNP'), continue; end
              out{i} = obj.SNPs{index(i)}.RS;
           end
        end
        function out = getValues(obj,prop)
           if obj.nSNP==0, out = []; return;end
           counter = 0;
           for i=1:1:obj.nSNP % look for first valid SNP
               if ~isempty(obj.SNPs{i})&&strcmp(obj.SNPs{i}.Type,'USNP')
                   tmp = eval(['obj.SNPs{i}.' prop]);
                   break;
               end
               counter = counter + 1;
           end
           if counter==obj.nSNP, out = []; return; end% there were no valid SNPs
           out = nan*zeros(size(tmp,1),size(tmp,2),obj.nSNP);
           for j=counter+1:1:obj.nSNP
              if isempty(obj.SNPs{j}), continue; end 
              if ~strcmp(obj.SNPs{j}.Type,'USNP'), continue; end
              out(:,:,j) = eval(['obj.SNPs{j}.' prop]);
           end
           out = squeeze(out);
        end
    end
    methods % SNP ANALYSIS
        function [MCOR] = illustrateSNPCOR(obj)
           MCOR = nan*zeros(obj.nSNP,obj.nSNP);
           %SNPCOR = nan*zeros(obj.nSNP,obj.nSNP);
           for i=1:1:obj.nSNP
               for j=1:1:obj.nSNP
                   MCOR(i,j) = angle(obj.SNPs{i}.MOD.M',obj.SNPs{j}.MOD.M');
                   %tmpc = corrcoef(obj.SNPs{i}.GT,obj.SNPs{j}.GT);
                   %SNPCOR(i,j) = tmpc(1,2);
               end
           end
           figure;imagesc(abs(MCOR));set(gca,'clim',[0 1]);
           %figure;imagesc(abs(SNPCOR));set(gca,'clim',[0 1]); 
        end
        function illustrateSNP(obj,index,ShapeSpace)
           disp(['RS: ' obj.SNPs{index}.RS ' ' obj.SNPs{index}.Label]);
           disp(['A: ' num2str(obj.SNPs{index}.TestA.AvgA) ' pA: ' num2str(obj.SNPs{index}.TestA.pA)]);
           disp(['R2: ' num2str(obj.SNPs{index}.TestR2.R2) ' pR2: ' num2str(obj.SNPs{index}.TestR2.pR2)]);
           disp(['EER: ' num2str(obj.SNPs{index}.TestBio.Soft.EER);]);
           mod = obj.SNPs{index}.MOD;
           figure;plot(mod.GT,mod.RIP,'b.');
           C = mod.M.*ShapeSpace.EigStd';
           vec = ShapeSpace.EigVec*C';
           vec = reshape(vec,3,length(vec)/3);
           vec = sqrt(sum(vec.^2,1));
           scan = clone(ShapeSpace.Average);
           scan.Value = vec;
           scan.ColorMode = 'Indexed';
           viewer(scan);
        end
        function out = biometricsSNP(obj,TestDepVar,TestGT,short,vis)
            % Initialize
             nrT = size(TestDepVar,1);
            % extracting RIP values
             RIP = nan*zeros(nrT,obj.nSNP);
             wSNPA = nan*zeros(obj.nSNP,1);
             wSNPP = nan*zeros(obj.nSNP,1);
             for i=1:1:obj.nSNP
                 obj.SNPs{i}.HTest = obj.HTest;
                 obj.SNPs{i}.RIPNorm = obj.RIPNorm;
                 tmp = USNP.updateRIP(TestDepVar,getMODM(obj.SNPs{i},obj.SNPs{i}.MOD))';
                 if obj.SNPs{i}.RIPNorm, tmp = USNP.normalizeRIP(tmp,getMODM(obj.SNPs{i},obj.SNPs{i}.MOD),obj.SNPs{i}.MOD.Var);end
                 RIP(:,i) = tmp;
                 wSNPA(i) = max(obj.SNPs{i}.A,0);
                 wSNPP(i) = -log(obj.SNPs{i}.pA);
             end
             out.wSNPA = wSNPA;out.wSNPP = wSNPP;
            % extracting matching scores 
             TMatches = nan*zeros(obj.nSNP,nrT);
             FMatches = nan*zeros(obj.nSNP,nrT,nrT-1);
             wSNP = nan*zeros(obj.nSNP,nrT);
             parfor t=1:nrT
                %disp(num2str(t));
                % t=1;
                tInd = setdiff(1:nrT,t);
                tmpFMatches = nan*zeros(obj.nSNP,nrT-1);
                tmpTMatches = nan*zeros(obj.nSNP,1);
                tmpwSNP = nan*zeros(obj.nSNP,1);
                for i=1:1:obj.nSNP
                    % i=1;
                    [Distr,freq] = GT2RIP(obj.SNPs{i},TestGT(t,i));% DOUBLE CHECK
                    tmpwSNP(i,1) = 1-freq;
                    m = getMatch(obj,Distr,RIP(:,i)); %#ok<*PFBNS>
                    tmpFMatches(i,:) = m(tInd);
                    tmpTMatches(i,:) = m(t);
                end
                TMatches(:,t) = tmpTMatches;
                FMatches(:,t,:) = tmpFMatches;
                wSNP(:,t) = tmpwSNP;
             end
             out.wSNPF = wSNP;
             switch obj.Classify
               case 'SOFT'
                   TScores = match2score(obj,TMatches);out.TScores = TScores;
                   FScores = match2score(obj,FMatches);out.FScores = FScores;
               case 'HARD'
                   TScores = BASE.class2score(match2class(obj,TMatches));out.TScores = TScores;
                   FScores = BASE.class2score(match2class(obj,FMatches));out.FScores = FScores;
             end
             if short, return; end
             nr = 4;
             EER = nan*zeros(1,nr);Ranks = nan*zeros(1,nr,nrT);
             X = nan*zeros(nr,(nrT*nrT)-1);Y = nan*zeros(nr,(nrT*nrT)-1);
             r = 1;% UNWEIGHTED
             [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BIOMETRICDNAMOD.getBiometrics(TScores,FScores,ones(obj.nSNP,nrT),ones(obj.nSNP,1));
             r = 2;% ANGLE WEIGHTED
             [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BIOMETRICDNAMOD.getBiometrics(TScores,FScores,ones(obj.nSNP,nrT),wSNPA);
             r = 3;% FREQUENCY WEIGHTED
             [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BIOMETRICDNAMOD.getBiometrics(TScores,FScores,wSNP,ones(obj.nSNP,1));
             r = 4;% ANGLE & FREQUENCY WEIGHTED
             [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BIOMETRICDNAMOD.getBiometrics(TScores,FScores,wSNP,wSNPA);
             Ranks(1,:,:) = (Ranks(1,:,:)./nrT).*100;
             CumRanks = zeros(1,nr,100);
             for cr=1:1:100
                 CumRanks(:,:,cr) = (sum(Ranks<=cr,3)/nrT).*100;
             end
             out.EER = EER;
             out.CumRanks = CumRanks;
             if vis
               str = {'b-.' 'b-' 'g-.' 'g-'};
               figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
               xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
               title('Identification SNP');
               plot(1:100,1:100,'k--');
               for i=1:1:nr
                   plot(1:100,squeeze(CumRanks(1,i,:)),str{i},'LineWidth',1);
               end
               figure;hold on;title('Verification SNP');
               xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
               plot(0:0.01:1,0:0.01:1,'k--');
               plot(0:0.01:1,1:-0.01:0,'k-');
               for i=1:1:nr
                   plot(X(i,:),Y(i,:),str{i},'LineWidth',1);
               end
             end
        end
    end
    methods % BASE ANALYSIS
        function out = biometricsBASE(obj,TestDepVar,TestCOV,TestGB,COVind,GBind,short,vis,AngleTable)
                 L1CL1TestDepVar = BASEMOD.getL1CL1DepVar(TestDepVar);
                 nrT = size(L1CL1TestDepVar,1);
                 nrMOD = obj.nBase;nrCOV = length(COVind);nrGB = length(GBind);
                 out.TScoresCOV = nan*zeros(nrCOV*nrMOD,nrT);
                 out.TScoresGB = nan*zeros(nrGB*nrMOD,nrT);
                 out.FScoresCOV = nan*zeros(nrCOV*nrMOD,nrT,nrT-1);
                 out.FScoresGB = nan*zeros(nrGB*nrMOD,nrT,nrT-1);
                 out.wCOVA = nan*zeros(nrCOV*nrMOD,1);
                 out.wCOVP = nan*zeros(nrCOV*nrMOD,1);
                 out.wCOVPR = nan*zeros(nrCOV*nrMOD,1);
                 out.wCOVF = nan*zeros(nrCOV*nrMOD,nrT);
                 out.wGBA = nan*zeros(nrGB*nrMOD,1);
                 out.wGBP = nan*zeros(nrGB*nrMOD,1);
                 out.wGBPR = nan*zeros(nrGB*nrMOD,1);
                 out.wGBF = nan*zeros(nrGB*nrMOD,nrT);
                 COVcounter = 1;
                 GBcounter = 1;
                 for i=1:1:obj.nBase
                     % i = 1;
                     disp(num2str(i));
                     obj.Base{i}.Match = obj.Match;
                     obj.Base{i}.Classify = obj.Classify;
                     obj.Base{i}.HTest = obj.HTest;
                     obj.Base{i}.RIPNorm = obj.RIPNorm;
                     obj.Base{i}.W = obj.W;
                     subres = biometrics(obj.Base{i},TestDepVar,TestCOV,TestGB,COVind,true,false,AngleTable);
                     out.TScoresCOV(COVcounter:(COVcounter+nrCOV-1),:) = subres.TScoresCOV(COVind,:);
                     out.FScoresCOV(COVcounter:(COVcounter+nrCOV-1),:,:) = subres.FScoresCOV(COVind,:,:);
                     out.wCOVA(COVcounter:(COVcounter+nrCOV-1),:) = subres.wCOVA(COVind,:);
                     out.wCOVP(COVcounter:(COVcounter+nrCOV-1),:) = subres.wCOVP(COVind,:);
                     out.wCOVPR(COVcounter:(COVcounter+nrCOV-1),:) = subres.wCOVPR(COVind,:);
                     out.wCOVF(COVcounter:(COVcounter+nrCOV-1),:) = subres.wCOVF(COVind,:);
                     
                     out.TScoresGB(GBcounter:(GBcounter+nrGB-1),:) = subres.TScoresGB(GBind,:);
                     out.FScoresGB(GBcounter:(GBcounter+nrGB-1),:,:) = subres.FScoresGB(GBind,:,:);
                     out.wGBA(GBcounter:(GBcounter+nrGB-1),:) = subres.wGBA(GBind,:);
                     out.wGBP(GBcounter:(GBcounter+nrGB-1),:) = subres.wGBP(GBind,:);
                     out.wGBPR(GBcounter:(GBcounter+nrGB-1),:) = subres.wGBPR(GBind,:);
                     out.wGBF(GBcounter:(GBcounter+nrGB-1),:) = subres.wGBF(GBind,:);
                     
                     COVcounter = COVcounter+nrCOV;
                     GBcounter = GBcounter+nrGB;
                 end
                 if short; return; end
                 switch obj.W
                    case 'A'
                        wCOV = out.wCOVA;
                        wGB = out.wGBA;
                    case 'P'
                        wCOV = out.wCOVP;
                        wGB = out.wGBP;
                    case 'PR'
                        wCOV = out.wCOVPR;
                        wGB = out.wGBPR;
                 end
                 nr = 6;
                 EER = nan*zeros(1,nr);Ranks = nan*zeros(1,nr,nrT);AUC = nan*zeros(1,nr);
                 X = nan*zeros(nr,(nrT*nrT)-1);Y = nan*zeros(nr,(nrT*nrT)-1);
                 r = 0;
                 r=r+1;% COVARIATES UNWEIGHTED
                 [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics(out.TScoresCOV,out.FScoresCOV,out.wCOVF,ones(size(wCOV)));
                 %out1 = BASEMOD.evalSortedResults(ST,SF,TestDepVar,T);
                 r=r+1;% COVARIATES WEIGHTED
                 [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics(out.TScoresCOV,out.FScoresCOV,out.wCOVF,wCOV);
                 out1 = BIOMETRICDNAMOD.evalSortedResults(ST,SF,L1CL1TestDepVar,T);
                 r=r+1;% GENETIC BACKGROUND UNWEIGHTED
                 [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics(out.TScoresGB,out.FScoresGB,out.wGBF,ones(size(wGB)));
                 %out1 = BASEMOD.evalSortedResults(ST,SF,TestDepVar,T);
                 r=r+1;% GENETIC BACKGROUND WEIGHTED
                 [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics(out.TScoresGB,out.FScoresGB,out.wGBF,wGB);
                 out1 = BIOMETRICDNAMOD.evalSortedResults(ST,SF,L1CL1TestDepVar,T);
                 r=r+1;% COMBINED UNWEIGHTED
                 [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics([out.TScoresCOV;out.TScoresGB],[out.FScoresCOV; out.FScoresGB],[out.wCOVF;out.wGBF],ones(size([wCOV;wGB])));
                 %out1 = BASEMOD.evalSortedResults(ST,SF,TestDepVar,T);
                 r=r+1;% COMBINED WEIGHTED
                 [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics([out.TScoresCOV;out.TScoresGB],[out.FScoresCOV; out.FScoresGB],[out.wCOVF;out.wGBF],[wCOV;wGB]);
                 out1 = BIOMETRICDNAMOD.evalSortedResults(ST,SF,L1CL1TestDepVar,T);
                 Ranks(1,:,:) = (Ranks(1,:,:)./nrT).*100;
                 CumRanks = zeros(1,nr,100);
                 for cr=1:1:100
                     CumRanks(:,:,cr) = (sum(Ranks<=cr,3)/nrT).*100;
                 end
                 out.EER = EER;
                 out.AUC = AUC;
                 out.CumRanks = CumRanks;
                 if vis
                    str = {'g-.' 'g-' 'b-.' 'b-' 'k-.' 'k-'};
                    figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
                    xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
                    title('Identification BASEMOD');
                    plot(1:100,1:100,'k--');
                    for i=1:1:nr
                        plot(1:100,squeeze(CumRanks(1,i,:)),str{i},'LineWidth',1);
                    end
                    figure;hold on;title('Verification BASEMOD');
                    xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
                    plot(0:0.01:1,0:0.01:1,'k--');
                    plot(0:0.01:1,1:-0.01:0,'k-');
                    for i=1:1:nr
                        plot(X(i,:),Y(i,:),str{i},'LineWidth',1);
                    end
                 end                  
        end
    end
    methods % BIOMETRICS
       function out = getMatch(obj,distr,rip)
                if isempty(distr), out = nan*zeros(1,length(rip)); return; end
                pT = normpdf(rip,distr.Tmu,distr.Tsigma)';
                pF = normpdf(rip,distr.Fmu,distr.Fsigma)';
                switch obj.Match
                    case 'QD'
                        out = pT./pF;
                    case 'INLIER'
                        out = pT./(pT+pF);
                end  
       end
       function out = match2score(obj,matches)
          switch obj.Match
              case 'QD'
                  out = -log(matches);
                  %out = 1./matches;
              case 'INLIER'
                  out = 1-matches;
          end
       end
       function out = match2class(obj,matches)
          switch obj.Match
              case 'QD'
                  T = 1;
              case 'INLIER'
                  T = 0.5;
          end
          out = BIOMETRICDNAMOD.matchthresholding(matches,T);
       end
       function out = biometrics(obj,TestDepVar,TestCOV,TestGB,TestGT,COVind,vis)
             disp('Running Biometrics BASE');
             BioBase = biometricsBASE(obj,TestDepVar,TestCOV,TestGB,COVind,true,false);
             disp('Running Biometrics SNPs');
             BioSNP = biometricsSNP(obj,TestDepVar,TestGT,true,false);
           
             %COVind = 1;
             wF = [BioBase.wCOVF(COVind,:);BioBase.wGBF;BioSNP.wSNPF];
             wA = [BioBase.wCOVA(COVind);BioBase.wGBA;BioSNP.wSNPA];
             TScores = [BioBase.TScoresCOV(COVind,:);BioBase.TScoresGB;BioSNP.TScores];
             FScores = [BioBase.FScoresCOV(COVind,:,:);BioBase.FScoresGB;BioSNP.FScores];
           
           
             nr = 4;nrT = size(TestDepVar,1);
             EER = nan*zeros(1,nr);Ranks = nan*zeros(1,nr,nrT);AUC = nan*zeros(1,nr);
             X = nan*zeros(nr,(nrT*nrT)-1);Y = nan*zeros(nr,(nrT*nrT)-1);
             r = 1;% SNPs
             [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics(BioSNP.TScores,BioSNP.FScores,BioSNP.wSNPF,BioSNP.wSNPA);
             out1 = BIOMETRICDNAMOD.evalSortedResults(ST,SF,TestDepVar,T);
             r = r+1;% BASE GB
             [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics(BioBase.TScoresGB,BioBase.FScoresGB,ones(obj.Base.nGB,nrT),BioBase.wGBA);
             out1 = BIOMETRICDNAMOD.evalSortedResults(ST,SF,TestDepVar,T);
             r = r+1;% BASE COV
             [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics(BioBase.TScoresCOV(COVind,:),BioBase.FScoresCOV(COVind,:,:),ones(length(COVind),nrT),BioBase.wCOVA(COVind));
             out1 = BIOMETRICDNAMOD.evalSortedResults(ST,SF,TestDepVar,T);
%              r = 4;% BASE
%              [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BIOMETRICDNAMOD.getBiometrics([BioBase.TScoresCOV(COVind,:); BioBase.TScoresGB],[BioBase.FScoresCOV(COVind,:,:); BioBase.FScoresGB],...
%                                                          ones(obj.Base.nGB+length(COVind),nrT),[BioBase.wCOVA(COVind);BioBase.wGBA]);
%              r = 5;% BASE GB + SNPS
%              [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BIOMETRICDNAMOD.getBiometrics([BioBase.TScoresGB;BioSNP.TScores],[BioBase.FScoresGB;BioSNP.FScores],...
%                                                          [ones(obj.Base.nGB,nrT);BioSNP.wSNPF],[BioBase.wGBA;BioSNP.wSNPA]);                                                                                            
%              r = r+1;% UNWEIGHTED
%              [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),T] = BIOMETRICDNAMOD.getBiometrics(TScores,FScores,wF,ones(size(wF,1),1));
%              out = BIOMETRICDNAMOD.evalSortedResults(ST,SF,DepVar,T);
             r = r+1;% ANGLE WEIGHTED
             [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICDNAMOD.getBiometrics(TScores,FScores,wF,wA);
             out1 = BIOMETRICDNAMOD.evalSortedResults(ST,SF,TestDepVar,T);
             
             clear out;
             Ranks(1,:,:) = (Ranks(1,:,:)./nrT).*100;
             CumRanks = zeros(1,nr,100);
             for cr=1:1:100
                 CumRanks(:,:,cr) = (sum(Ranks<=cr,3)/nrT).*100;
             end
             out.EER = EER;
             out.CumRanks = CumRanks;
             out.AUC = AUC;
             if vis
               str = {'g-' 'b-' 'r-' 'k-.' 'k-'};
               figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
               xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
               title('Identification COMB');
               plot(1:100,1:100,'k--');
               for i=1:1:nr
                   plot(1:100,squeeze(CumRanks(1,i,:)),str{i},'LineWidth',1);
               end
               figure;hold on;title('Verification COMB');
               xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
               plot(0:0.01:1,0:0.01:1,'k--');
               plot(0:0.01:1,1:-0.01:0,'k-');
               for i=1:1:nr
                   plot(X(i,:),Y(i,:),str{i},'LineWidth',1);
               end
             end
       end
    end
    methods (Static = true)
       function out = matchthresholding(matches,T)
                out = matches>=T; 
       end
       function out = class2score(classes)
                out = 1-classes; 
       end
       function [R,EER,ST,SF,X,Y,AUC,T] = getBiometrics(TScores,FScores,wM,W)
                 nF = size(FScores,3);
                 nT = size(TScores,2);
                 tmpW = repmat(W,1,nT).*wM;
                 ST = BASE.accScores(TScores,tmpW)';
                 SF = BASE.accScores(FScores,repmat(tmpW,1,1,nF));
                 [R,EER,X,Y,AUC,T] = BIOMETRICDNAMOD.evalScores(ST,SF,nT,nF); 
       end
       function Score = accScores(Scores,tmpW)
                Score = squeeze(nansum(tmpW.*Scores,1)./nansum(tmpW,1));
       end
       function [R,EER,x,y,AUC,T] = evalScores(ST,SF,nT,nF)
                 % Rank analysis
                 R = (sum(SF<repmat(ST,1,nF),2)+1);%/(nF+1);
                 % EER analysis
                 g = [ones(1,nT), zeros(1,nT*nF)];
                 [tmp,order] = sort([ST;SF(:)]);
                 g = g(order);
                 true_neg = g == 0;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
                 true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
                 y = tpf(1:end-1)+dy./2;
                 x = fpf(1:end-1)+dx./2;
                 yn = 1-y;d = abs(x-yn);
                 [~,indmin] = min(d);
                 EER = ((x(indmin)+yn(indmin))/2);
                 AUC = sum(dx.*y);
                 T = tmp(indmin);
       end
       function [A,D] = evalSortedResults(ST,SF,DepVar,T)
                % DepVar = TestDepVar;
                nrT = size(DepVar,1);
                D = nan*zeros(nrT,nrT-1);
                A = nan*zeros(nrT,nrT-1);
                AT = nan*zeros(nrT,2);
                for t=1:1:nrT
                    %t=1;
                    tind = setdiff(1:nrT,t);
                    forD = pdist2(DepVar(t,:),DepVar(tind,:),'euclidean');
                    %forA = -1*(pdist2(DepVar(t,:),DepVar(tind,:),'cosine')-1);
                    forA = (pdist2(DepVar(t,:),DepVar(tind,:),'cosine'));
                    [~,ind] = sort(SF(t,:),'ascend');
                    D(t,:) = forD(ind);
                    A(t,:) = forA(ind);
                    AT(t,1) = mean(forA(SF(t,:)<=T));
                    AT(t,2) = mean(forA(SF(t,:)>=T));   
                end
                D = mean(D,1);
                A = mean(A,1);
                
                %figure;plot(1:nrT-1,D,'b-');
                figure;plot(((1:nrT-1)/(nrT-1))*100,A,'b-');hold on;
                set(gca,'ylim',[0.95 1.05],'xlim',[0 100]);
                plot(((1:nrT-1)/(nrT-1))*100,mean(A)*ones(1,nrT-1),'k-');
                
                p = polyfit(1:nrT-1,A,3);
                plot(((1:nrT-1)/(nrT-1))*100,polyval(p,1:nrT-1),'r-','lineWidth',2);
                
                
                figure;plot(((1:nrT-1)/(nrT-1))*100,D,'b-');hold on;
%                 set(gca,'ylim',[10 12],'xlim',[0 100]);
%                 plot(((1:nrT-1)/(nrT-1))*100,mean(D)*ones(1,nrT-1),'k-');
%                 
%                 p = polyfit(1:nrT-1,D,4);
%                 plot(((1:nrT-1)/(nrT-1))*100,polyval(p,1:nrT-1),'r-','lineWidth',2);
                
                %figure;boxplot(A);set(gca,'ylim',[0.9 1.1]);
       end
    end
end