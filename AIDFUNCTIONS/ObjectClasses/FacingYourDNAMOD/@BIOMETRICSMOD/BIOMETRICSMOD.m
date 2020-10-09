classdef BIOMETRICSMOD < superClassLight
    properties
       PPM = {};
    end
    properties (Dependent = true)
       nPPM;
       nSEX;
       nAGE;
       nW;
       nH;
       nCOV;
       nGB;
       nUnGB;
       nSNP;
       nUnSNP;
    end
    properties % MATCHING
       Match = 'QD';
       Classify = 'SOFT';
       AngleWeights = 'PR';
       BalOCC = true;
       FreqWeights = true;
       Kappa = 6;
       HTest = true;
       RIPNorm = true;
    end
    properties (Hidden = true)
       IDs;
       QDTMatches;
       QDFMatches;
       INTMatches;
       INFMatches;
       wF;
       wA;
       wP;
       wPR;
    end
    properties (Hidden = true, Dependent = true)
       SEXInd;
       AGEInd;
       WInd;
       HInd;
       COVInd;
       GBInd;
       SNPInd;
    end
    methods % CONSTRUCTOR
        function obj = BIOMETRICSMOD(varargin)
            obj = obj@superClassLight(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.nPPM(obj)
           if isempty(obj.PPM), out = 0; return; end
           out = length(obj.PPM);
        end
        function out = get.SEXInd(obj)
            tmp = SEXMOD;
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nSEX(obj)
            out = length(obj.SEXInd); 
        end
        function out = get.AGEInd(obj)
            tmp = COVMOD;tmp.VarName = 'AGE';
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nAGE(obj)
            out = length(obj.AGEInd); 
        end
        function out = get.WInd(obj)
            tmp = COVMOD;tmp.VarName = 'W';
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nW(obj)
            out = length(obj.WInd); 
        end
        function out = get.HInd(obj)
            tmp = COVMOD;tmp.VarName = 'H';
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nH(obj)
            out = length(obj.HInd); 
        end
        function out = get.COVInd(obj)
            tmp = GBMOD;
            out = find(obj.IDs<tmp.ID);
        end
        function out = get.nCOV(obj)
            out = length(obj.COVInd); 
        end
        function out = get.GBInd(obj)
            tmp = GBMOD;
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nGB(obj)
            out = length(obj.GBInd); 
        end
        function out = get.SNPInd(obj)
            tmp = SNPMOD;
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nSNP(obj)
            out = length(obj.SNPInd); 
        end
        function obj = set.Match(obj,in)
            obj.Match = in;
            setMatch(obj,in);
        end
        function obj = set.Classify(obj,in)
            obj.Classify = in;
            setClassify(obj,in);
        end
        function obj = set.Kappa(obj,in)
           obj.Kappa = in;
           setValues(obj,'Kappa',in);
        end
        function obj = set.HTest(obj,in)
            obj.HTest = in;
            setValues(obj,'HTest',in);
        end
        function obj = set.RIPNorm(obj,in)
            obj.RIPNorm;
            setValues(obj,'RIPNorm',in);
        end
        function out = get.nUnSNP(obj)
                 out = length(getUniqueSNPs(obj)); 
        end
        function out = get.nUnGB(obj)
                 out = length(getUniqueGBs(obj)); 
        end
    end
    methods % INTERFACING
        function setMatch(obj,in)
            if obj.nPPM==0, return;end
            for i=1:1:obj.nPPM
                if isempty(obj.PPM{i}), continue;end
                obj.PPM{i}.Match = in;
            end
        end
        function setClassify(obj,in)
            if obj.nPPM==0, return;end
            for i=1:1:obj.nPPM
                if isempty(obj.PPM{i}), continue;end
                obj.PPM{i}.Classify = in;
            end
        end
        function setValues(obj,prop,in)
            if obj.nPPM==0, return;end
            for i=1:1:obj.nPPM
                if isempty(obj.PPM{i}), continue;end
                if ~isprop(obj.PPM{i},prop), continue; end
                eval(['obj.PPM{i}.' prop ' = ' num2str(in) ';']);
            end
        end
        function out = getValues(obj,prop)
           if obj.nPPM==0, out = []; return;end
           counter = 0;
           for i=1:1:obj.nPPM % look for first valid SNP
               if ~isempty(obj.PPM{i})&&isprop(obj.PPM{i},prop)
                   tmp = eval(['obj.PPM{i}.' prop]);
                   break;
               end
               counter = counter + 1;
           end
           if counter==obj.nPPM, out = []; return; end% there were no valid SNPs
           out = nan*zeros(size(tmp,1),size(tmp,2),obj.nPPM);
           for j=counter+1:1:obj.nPPM
              if isempty(obj.PPM{j}), continue; end
              if ~isprop(obj.PPM{j},prop), continue; end
              out(:,:,j) = eval(['obj.PPM{j}.' prop]);
           end
           out = squeeze(out);
        end
        function out = getVarNames(obj,index)
           if obj.nPPM==0, out = []; return;end
           if nargin < 2, index = 1:obj.nPPM; end
           out = cell(1,length(index));
           for i=1:1:length(index)
              if isempty(obj.PPM{index(i)}), continue; end 
              out{i} = obj.PPM{index(i)}.VarName;
           end
        end
        function out = getUniqueSNPs(obj)
             if obj.nSNP==0, out = []; return;end
             out = unique(getVarNames(obj,obj.SNPInd));
        end
        function out = getUniqueGBs(obj)
             if obj.nGB==0, out = []; return;end
             out = unique(getValues(obj,'GBID'));
             out = out(~isnan(out));
        end
        function out = getIDs(obj,index)
           if obj.nPPM==0, out = []; return;end 
           if nargin < 2, index = 1:obj.nPPM; end
           out = nan*zeros(1,length(index));
           for i=1:1:length(index)
              if isempty(obj.PPM{index(i)}), continue; end 
              out(i) = obj.PPM{index(i)}.ID;
           end 
        end
        function storeIDs(obj)
           obj.IDs = getIDs(obj); 
        end
        function reducePPM(obj,index)
           obj.PPM = obj.PPM(index);
           storeIDs(obj);
        end
        function out = getTestX(obj,COV,GB,GT,RS)
            nrT = size(COV,1);
            out = nan*zeros(nrT,obj.nPPM);
            for i=1:1:obj.nPPM
               out(:,i) = BIOMETRICSMOD.ppmTestX(obj.PPM{i},COV,GB,GT,RS,nrT);
            end
        end
        function out = match2score(obj,matches)
          switch obj.Match
              case 'QD'
                  out = -log(matches);
              case 'INLIER'
                  out = 1-matches;
          end
        end
        function out = match2class(obj,matches)
          switch obj.Match
              case 'QD'
                  T = 1;
              case 'INLIER'
                  %T = 0.5;
                  T = 0.5*ones(obj.nPPM,1);
                  %T = getValues(obj,'TH')';
                  %T = T(:,2);
          end
          out = BIOMETRICSMOD.matchthresholding(matches,T);
        end
        function out = prunePPM(obj,nr)
            names = getVarNames(obj);
            pAR = getValues(obj,'pAR');
            unnames = unique(names);
            nrunnames = length(unnames);
            out = nan*zeros(nrunnames,nr);
            for i=1:1:nrunnames
                ind = find(strcmp(unnames{i},names));
                [~,subind] = sort(pAR(ind),'ascend');
                if length(subind)>=nr
                   out(i,:) = ind(subind(1:nr));
                else
                   out(i,1:length(subind)) = ind(subind); 
                end
            end
            out = out(:);
            out = out(find(~isnan(out)));
        end
        function [occ,ppmIDs] = getPPMOcc(obj)
                 if obj.nPPM==0, occ = []; ppmIDs = []; return;end
                 VAR = getVarNames(obj);
                 var = unique(VAR);
                 ppmIDs = nan*zeros(1,obj.nPPM);
                 for i=1:1:length(var)
                    ppmIDs(strcmp(var{i},VAR)) = i; 
                 end
                 occ = nan*zeros(1,obj.nPPM);
                 for i=1:1:obj.nPPM
                    occ(i) = sum(ppmIDs==ppmIDs(i));
                 end
        end
    end
    methods % BIOMETRICS
        function biometricMatchesWeights(obj,TestDepVar,TestCOV,TestGB,TestGT,TestRS,AngleTable)
               % INITIALIZING
                 L1CL1TestDepVar = BIOMETRICSMOD.getL1CL1DepVar(TestDepVar);
                 nrT = size(L1CL1TestDepVar,1);
                 %TestX = getTestX(obj,TestCOV,TestGB,TestGT,TestRS);
                 %QDTM = nan*zeros(obj.nPPM,nrT);
                 %QDFM = nan*zeros(obj.nPPM,nrT,nrT-1);
                 INTM = nan*zeros(obj.nPPM,nrT);
                 INFM = nan*zeros(obj.nPPM,nrT,nrT-1);
                 F = nan*zeros(obj.nPPM,nrT);
                 A = nan*zeros(obj.nPPM,1);
                 P = nan*zeros(obj.nPPM,1);
                 PR = nan*zeros(obj.nPPM,1);
               % RUNNING
                 f=waitbar(0,'MATCHING');
                 for i=1:1:obj.nPPM
                     %if mod(i,20)==0, disp([num2str(round((i/obj.nPPM)*100)) ' %']);end
                     if mod(i,20)==0, waitbar(i/obj.nPPM,f);drawnow; end
                     %[QDTM(i,:),QDFM(i,:,:),INTM(i,:),INFM(i,:,:),F(i,:),A(i),P(i),PR(i)] = ...
                     [~,~,INTM(i,:),INFM(i,:,:),F(i,:),A(i),P(i),PR(i)] = ...
                                                                                biometricMatchesWeights(obj.PPM{i},TestDepVar,...
                                                                                BIOMETRICSMOD.ppmTestX(obj.PPM{i},TestCOV,TestGB,TestGT,TestRS,nrT),...
                                                                                AngleTable);
                 end
                 delete(f);
                 %obj.QDTMatches = QDTM;
                 %obj.QDFMatches = QDFM;
                 obj.INTMatches = INTM;
                 obj.INFMatches = INFM;
                 obj.wF = F;
                 obj.wA = A;
                 obj.wP = P;
                 obj.wPR = PR;
        end
        function out = biometricPerformance(obj,TestDepVar,vis,index)
                 if nargin<4,index = {1:obj.nPPM};end
                 if isempty(index),index = {1:obj.nPPM};end
                 if ~iscell(index), index = {index};end
                 L1CL1TestDepVar = BIOMETRICSMOD.getL1CL1DepVar(TestDepVar);
                 nrT = size(L1CL1TestDepVar,1);
                 nG = length(index);
                 nr = nG*2;
                 EER = nan*zeros(1,nr);Ranks = nan*zeros(1,nr,nrT);AUC = nan*zeros(1,nr);
                 X = nan*zeros(nr,(nrT*nrT)-1);Y = nan*zeros(nr,(nrT*nrT)-1);
                 FM = cell(1,nG);
                 r = 0;
                 switch obj.AngleWeights
                     case 'A'
                         wa = obj.wA;
                     case 'P'
                         wa = obj.wP;
                     case 'PR'
                         wa = obj.wPR;
                     otherwise
                         wa = ones(size(obj.wA));
                 end
                 if obj.BalOCC
                    occ = getPPMOcc(obj)';
                    %wa = wa.*(1./occ);
                    wa = wa.*occ;
                 end
                 switch obj.FreqWeights
                     case true
                         wf = obj.wF;
                     case false
                         wf = ones(size(obj.wF));
                 end
                 switch obj.Match
                     case 'QD'
                         TMatches = obj.QDTMatches;
                         FMatches = obj.QDFMatches;
                     case 'INLIER'
                         TMatches = obj.INTMatches;
                         FMatches = obj.INFMatches;
                 end
                 switch obj.Classify
                   case 'SOFT'
                       TScores = match2score(obj,TMatches);
                       FScores = match2score(obj,FMatches);
                   case 'HARD'
                       TScores = BIOMETRICSMOD.class2score(match2class(obj,TMatches));
                       FScores = BIOMETRICSMOD.class2score(match2class(obj,FMatches));
                 end
                 for i=1:1:nG 
                     r=r+1;ind = index{i};
%                      [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:),AUC(:,r),~] = BIOMETRICSMOD.getBiometrics(TScores(ind,:),FScores(ind,:,:),wf(ind,:),ones(size(wa(ind))));
                     r=r+1;
                     [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICSMOD.getBiometrics(TScores(ind,:),FScores(ind,:,:),wf(ind,:),wa(ind));
                     [FM{i}.A,FM{i}.evalA,FM{i}.D,FM{i}.evalD,FM{i}.X] = BIOMETRICSMOD.evalSortedResults(ST,SF,L1CL1TestDepVar,T);
                 end
                 Ranks(1,:,:) = (Ranks(1,:,:)./nrT).*100;
                 CumRanks = zeros(1,nr,100);
                 for cr=1:1:100
                     CumRanks(:,:,cr) = (sum(Ranks<=cr,3)/nrT).*100;
                 end
                 out.EER = EER;
                 out.AUC = AUC;
                 out.CumRanks = CumRanks;
                 if vis
                    str = {'g-.' 'g-' 'b-.' 'b-' 'r-.' 'r-' 'k-.' 'k-' 'm-.' 'm-' 'c-.' 'c-'}; 
                    figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
                    xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
                    title('Identification');
                    plot(1:100,1:100,'k--');
                    for i=1:1:nr
                        plot(1:100,squeeze(CumRanks(1,i,:)),str{i},'LineWidth',2);
                    end
                    figure;hold on;title('Verification');
                    xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
                    plot(0:0.01:1,0:0.01:1,'k--');
                    plot(0:0.01:1,1:-0.01:0,'k-');
                    for i=1:1:nr
                        plot(X(i,:),Y(i,:),str{i},'LineWidth',2);
                    end
                    meanA = mean(FM{1}.A);
                    meanD = mean(FM{1}.D);
                    figure;set(gca,'ylim',[meanA-0.05 meanA+0.05],'xlim',[0 100]);hold on;grid on;
                    plot(((1:nrT-1)/(nrT-1))*100,meanA*ones(1,nrT-1),'k-');
                    title('False Matches Angle');
                    counter = 0;
                    for i=1:1:nG
                        counter = counter+2;
                        plot(FM{i}.X,FM{i}.evalA,str{counter},'lineWidth',2);
                    end
                    figure;set(gca,'ylim',[meanA-0.05 meanA+0.05],'xlim',[0 100]);hold on;grid on;
                    plot(((1:nrT-1)/(nrT-1))*100,meanA*ones(1,nrT-1),'k-');
                    title('False Matches Angle');
                    counter = 0;
                    for i=1:1:nG
                        counter = counter+2;
                        plot(FM{i}.X,FM{i}.A,str{counter});
                    end
                    figure;set(gca,'ylim',[meanD-2 meanD+2],'xlim',[0 100]);hold on;grid on;
                    plot(((1:nrT-1)/(nrT-1))*100,meanD*ones(1,nrT-1),'k-');
                    title('False Matches Distance');
                    counter = 0;
                    for i=1:1:nG
                        counter = counter+2;
                        plot(FM{i}.X,FM{i}.D,str{counter});
                    end
                    figure;set(gca,'ylim',[meanD-2 meanD+2],'xlim',[0 100]);hold on;grid on;
                    plot(((1:nrT-1)/(nrT-1))*100,meanD*ones(1,nrT-1),'k-');
                    title('False Matches Distance');
                    counter = 0;
                    for i=1:1:nG
                        counter = counter+2;
                        plot(FM{i}.X,FM{i}.evalD,str{counter},'lineWidth',2);
                    end      
                 end
        end
    end
    methods (Static = true)
        function out = ppmTestX(ppm,COV,GB,GT,RS,nrT)
            if isempty(ppm),out = nan*zeros(nrT,1);return;end
            switch ppm.ID
                case {1 2 3 4 5 6 7 8 9 10}
                    out = COV(:,ppm.ID);
                case 100
                    out = GB(:,ppm.GBID);
                case 1000
                    ind = find(strcmp(ppm.RS,RS));
                    if isempty(ind),out = nan*zeros(nrT,1);return;end
                    out = GT(:,ind);
                otherwise
                    out = nan*zeros(nrT,1);      
            end 
        end
        function [R,EER,ST,SF,X,Y,AUC,T] = getBiometrics(TScores,FScores,wM,W)
                 nF = size(FScores,3);
                 nT = size(TScores,2);
                 tmpW = repmat(W,1,nT).*wM;
                 ST = BIOMETRICSMOD.accScores(TScores,tmpW)';
                 SF = BIOMETRICSMOD.accScores(FScores,repmat(tmpW,1,1,nF));
                 [R,EER,X,Y,AUC,T] = BIOMETRICSMOD.evalScores(ST,SF,nT,nF); 
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
        function [A,evalA,D,evalD,X] = evalSortedResults(ST,SF,DepVar,T)
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
                X = ((1:nrT-1)/(nrT-1))*100;
                polyA = polyfit(1:nrT-1,A,3);
                evalA = polyval(polyA,1:nrT-1);
                polyD = polyfit(1:nrT-1,D,3);
                evalD = polyval(polyD,1:nrT-1);
                
%                 figure;plot(((1:nrT-1)/(nrT-1))*100,A,'b-');hold on;
%                 set(gca,'ylim',[0.94 1.06],'xlim',[0 100]);
%                 plot(((1:nrT-1)/(nrT-1))*100,mean(A)*ones(1,nrT-1),'k-');
%                 
%                 p = polyfit(1:nrT-1,A,3);
%                 plot(((1:nrT-1)/(nrT-1))*100,polyval(p,1:nrT-1),'r-','lineWidth',2);
%                 
%                 
%                 figure;plot(((1:nrT-1)/(nrT-1))*100,D,'b-');hold on;
%                 set(gca,'ylim',[10 13],'xlim',[0 100]);
% %                 set(gca,'ylim',[10 12],'xlim',[0 100]);
% %                 plot(((1:nrT-1)/(nrT-1))*100,mean(D)*ones(1,nrT-1),'k-');
% %                 
% %                 p = polyfit(1:nrT-1,D,4);
% %                 plot(((1:nrT-1)/(nrT-1))*100,polyval(p,1:nrT-1),'r-','lineWidth',2);
%                 
%                 %figure;boxplot(A);set(gca,'ylim',[0.9 1.1]);
        end
        function out = matchthresholding(matches,T)
            if size(matches,3)==1
                out = matches>=repmat(T,1,size(matches,2));
            else
                out = matches>=repmat(T,1,size(matches,2),size(matches,3));
            end
        end
        function out = class2score(classes)
                out = 1-classes; 
        end
        function out = getL1CL1DepVar(in)
           if ~iscell(in), out = in;return;end
           out = in{1}.DepVar{1};
       end
    end
end