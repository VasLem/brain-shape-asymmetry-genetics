classdef BIOMETRICSMV < superClassLight
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
       Classify = 'SOFT';
       HitWeights = true;
       FreqWeights = true;
       GroupNorm = true;
       SNPGroupNorm = true;
    end
    properties (Hidden = true)
       IDs;
       SOFTTMatches;
       SOFTFMatches;
       HARDTMatches;
       HARDFMatches;
       wF;
       wH;
       wL;
       SNPGroups;
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
        function obj = BIOMETRICSMV(varargin)
            obj = obj@superClassLight(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.nPPM(obj)
           if isempty(obj.PPM), out = 0; return; end
           out = length(obj.PPM);
        end
        function out = get.SEXInd(obj)
            tmp = SEXMV;
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nSEX(obj)
            out = length(obj.SEXInd); 
        end
        function out = get.AGEInd(obj)
            tmp = COVMV;tmp.VarName = 'AGE';
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nAGE(obj)
            out = length(obj.AGEInd); 
        end
        function out = get.WInd(obj)
            tmp = COVMV;tmp.VarName = 'W';
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nW(obj)
            out = length(obj.WInd); 
        end
        function out = get.HInd(obj)
            tmp = COVMV;tmp.VarName = 'H';
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nH(obj)
            out = length(obj.HInd); 
        end
        function out = get.COVInd(obj)
            tmp = GBMV;
            out = find(obj.IDs<tmp.ID);
        end
        function out = get.nCOV(obj)
            out = length(obj.COVInd); 
        end
        function out = get.GBInd(obj)
            tmp = GBMV;
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nGB(obj)
            out = length(obj.GBInd); 
        end
        function out = get.SNPInd(obj)
            tmp = SNPMV;
            out = find(obj.IDs==tmp.ID);
        end
        function out = get.nSNP(obj)
            out = length(obj.SNPInd); 
        end
        function obj = set.Classify(obj,in)
            obj.Classify = in;
            %setClassify(obj,in);
        end
        function out = get.nUnSNP(obj)
                 out = length(getUniqueSNPs(obj)); 
        end
        function out = get.nUnGB(obj)
                 out = length(getUniqueGBs(obj)); 
        end
    end
    methods % INTERFACING
        function setClassify(obj,in)
            return;
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
            nrT = size(GT,1);
            out = nan*zeros(nrT,obj.nPPM);
            for i=1:1:obj.nPPM
                if isempty(obj.PPM{i}), continue;end
               out(:,i) = getTestX(obj.PPM{i},COV,GB,GT,RS,nrT);
            end
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
        function out = getSNPGrouping(obj,T)
            val = getValues(obj,'GG');
            val = val(:,obj.SNPInd);
            val(val==0) = nan;
            n = size(val,2);
            %MI = zeros(n,n);
            NMI = zeros(n,n);
            parfor i=1:1:n
                v1 = val(:,i);
                ind1 = find(~isnan(v1));
                v1(v1==-1) = 2;
                %v1 = int8(v1);
                tmp = zeros(1,n);
                for j=1:1:n
                    v2 = val(:,j);
                    ind2 = find(~isnan(v2));
                     v2(v2==-1) = 2;
                     %v2 = int8(v2);
                     ind = intersect(ind1,ind2);
                     %MI(i,j) = mutualinfo(v1(ind),v2(ind));
                     [~,tmp(j)] = nmi(v1(ind),v2(ind));
                end
                NMI(i,:) = tmp;
            end
            %NMI = real(NMI);
            gg = zeros(1,n);
            list = 1:n;
            counter = 1;
            while ~isempty(list)
                new = list(1);
                mival = NMI(new,:);
                index = find(mival>=T);
                index = intersect(list,index);
                gg(index) = counter;
                list = setdiff(list,index);
                counter = counter +1;
            end
            out = gg;
        end
        function out = getSNPGroupingRS(obj)
            snpind = obj.SNPInd;
            VAR = getVarNames(obj);
            VAR = VAR(snpind);
            n = length(snpind);
            gg = zeros(1,n);
            list = 1:n;
            counter = 1;
            while ~isempty(list)
                new = VAR{list(1)};
                index = find(strcmp(new,VAR));
                index = intersect(list,index);
                gg(index) = counter;
                list = setdiff(list,index);
                counter = counter +1;
            end
            out = gg;
        end
        function fiterSNPPPM(obj,balth)
            % BADLY balanced SNP
            out = getBadBalancedSNP(obj,balth);
            if ~isempty(out), reducePPM(obj,setdiff(1:obj.nPPM,out));end
            % REDUNDANT SNP OCCURANCE
            out = getRedundantOccuranceSNP(obj);
            if ~isempty(out), reducePPM(obj,setdiff(1:obj.nPPM,out));end
            % SPLIT ADDITIVE MODELS
            %splitAdditiveSNP(obj);
            % BADLY balanced SNP, possibly induced by splitting additive
            % models
            %out = getBadBalancedSNP(obj,balth);
            %if ~isempty(out), reducePPM(obj,setdiff(1:obj.nPPM,out));end
        end
        function out = getBadBalancedSNP(obj,balth)
            out = [];
            if obj.nSNP==0, return; end
            snpind = obj.SNPInd;
            for i=1:obj.nSNP
                tmp = obj.PPM{snpind(i)}.GGGTBal;
                tmp = min(tmp);
                if tmp<=balth, out = [out snpind(i)]; end %#ok<AGROW>
            end
        end
        function out = getRedundantOccuranceSNP(obj)
            out = [];
            if obj.nSNP==0, return; end
            occ = getPPMOcc(obj);
            snpind = obj.SNPInd;
            VAR = getVarNames(obj);
            testind = getValues(obj,'TestInd');
            for i=1:1:obj.nSNP
                if occ(snpind(i))<2, continue; end
                if testind(snpind(i))>4, continue; end
                rs = VAR(snpind(i));
                ind = find(strcmp(rs,VAR));
                ind = setdiff(ind,snpind(i));
                for j=1:1:length(ind)
                    if testind(ind(j))>4, out = [out ind(j)];end %#ok<AGROW>
                end   
            end
            out = unique(out);
        end
        function splitAdditiveSNP(obj)
            if obj.nSNP==0, return; end
            keep = zeros(1,obj.nPPM);
            snpind = obj.SNPInd;
            VAR = getVarNames(obj);
            testind = getValues(obj,'TestInd');
            nr = sum(testind==4);
            add = cell(1,nr*3);
            counter = 0;
            f=waitbar(0,'SPLITTING');drawnow;
            for i=1:obj.nPPM
                % i=4
                if ~ismember(i,snpind), keep(i) = 1;continue; end % not a snp
                if ~(testind(i)==4), keep(i) = 1;continue; end
                rs = VAR(i);
                ind = find(strcmp(rs,VAR));
                ind = setdiff(ind,snpind(i));
                tests = setdiff(1:3,testind(ind));
                for j=1:1:length(tests)
                    tmp = clone(obj.PPM{i});
                    tmp.TestInd = tests(j);
                    counter = counter+1;
                    add{counter} = tmp;
                end
                waitbar(i/obj.nPPM,f);drawnow;
            end
            delete(f);
            add = add(1:counter);
            keep = obj.PPM(keep==1);
            obj.PPM = [keep add];
            storeIDs(obj);
        end
    end
    methods % BIOMETRICS
        function biometricMatchesWeights(obj,TestDepVar,TestCOV,TestGB,TestGT,TestRS)
               % INITIALIZING
                % TestGT = TestSNP;
                 L1CL1TestDepVar = BIOMETRICSMV.getL1CL1DepVar(TestDepVar);
                 nrT = size(L1CL1TestDepVar,1);
                 SOFTTM = nan*zeros(obj.nPPM,nrT);
                 SOFTFM = nan*zeros(obj.nPPM,nrT,nrT-1);
                 HARDTM = nan*zeros(obj.nPPM,nrT);
                 HARDFM = nan*zeros(obj.nPPM,nrT,nrT-1);
                 wF = nan*zeros(obj.nPPM,nrT); %#ok<*PROPLC>
                 wH = nan*zeros(obj.nPPM,1);
                 wL = nan*zeros(obj.nPPM,1);
               % RUNNING
                 f=waitbar(0,'MATCHING');drawnow;
                 ppm = obj.PPM;
                 %[path,ID] = setupParForProgress(length(ppm));
                 for i=1:1:length(ppm)
                     %if mod(i,1)==0, waitbar(i/obj.nPPM,f);drawnow; end
                     [SOFTTM(i,:),SOFTFM(i,:,:),HARDTM(i,:),HARDFM(i,:,:),wF(i,:),wH(i),wL(i)] = biometricMatchesWeights(ppm{i},TestDepVar,getTestX(ppm{i},TestCOV,TestGB,TestGT,TestRS,nrT));
                     %parfor_progress;
                     waitbar(i/obj.nPPM,f);drawnow;
                 end
                 %closeParForProgress(path,ID);
                 delete(f);
                 obj.SOFTTMatches = SOFTTM;
                 obj.SOFTFMatches = SOFTFM;
                 obj.HARDTMatches = HARDTM;
                 obj.HARDFMatches = HARDFM;
                 obj.wF = wF;
                 obj.wH = wH;
                 obj.wL = wL;
        end
        function out = biometricPerformance(obj,TestDepVar,vis,index)
                 if nargin<4,index = {1:obj.nPPM};end
                 if isempty(index),index = {1:obj.nPPM};end
                 if ~iscell(index), index = {index};end
                 L1CL1TestDepVar = BIOMETRICSMV.getL1CL1DepVar(TestDepVar);
                 nrT = size(L1CL1TestDepVar,1);
                 nG = length(index);
                 nr = nG*2;
                 EER = nan*zeros(1,nr);Ranks = nan*zeros(1,nr,nrT);AUC = nan*zeros(1,nr);
                 X = nan*zeros(nr,(nrT*nrT)-1);Y = nan*zeros(nr,(nrT*nrT)-1);
                 FM = cell(1,nG);
                 r = 0;
                 switch obj.HitWeights
                     case true
                         %wh = obj.wH;
                         wh = obj.wL;
                     case false
                         wh = ones(size(obj.wH));
                 end
                 switch obj.FreqWeights
                     case true
                         wf = obj.wF;
                     case false
                         wf = ones(size(obj.wF));
                         wf(isnan(obj.wF)) = nan;
                 end
                 wfh = repmat(wh,1,nrT).*wf;% combined frequency and hits weights
                 switch obj.GroupNorm
                     case true
                         sexind = obj.SEXInd;
                         ageind = obj.AGEInd;
                         hind = obj.HInd;
                         wind = obj.WInd;
                         gbind = obj.GBInd;
                         snpind = obj.SNPInd;
                         ggind = obj.SNPGroups;
                         for i=1:1:nrT
                             % SEX
                             if ~isempty(sexind)
                                wf(sexind,i) = BIOMETRICSMV.sum2one(wf(sexind,i)); 
                                wfh(sexind,i) = BIOMETRICSMV.sum2one(wfh(sexind,i)); 
                             end
                             % AGE
                             if ~isempty(ageind)
                                 wf(ageind,i) = BIOMETRICSMV.sum2one(wf(ageind,i)); 
                                 wfh(ageind,i) = BIOMETRICSMV.sum2one(wfh(ageind,i));
                             end
                             % HEIGHT 
                             if ~isempty(hind)
                                 wf(hind,i) = BIOMETRICSMV.sum2one(wf(hind,i)); 
                                 wfh(hind,i) = BIOMETRICSMV.sum2one(wfh(hind,i));
                             end
                             % WEIGHT 
                             if ~isempty(wind)
                                 wf(wind,i) = BIOMETRICSMV.sum2one(wf(wind,i)); 
                                 wfh(wind,i) = BIOMETRICSMV.sum2one(wfh(wind,i));
                             end
                             % GBBACKGROUND
                             if ~isempty(gbind)
                                 wf(gbind,i) = BIOMETRICSMV.sum2one(wf(gbind,i)); 
                                 wfh(gbind,i) = BIOMETRICSMV.sum2one(wfh(gbind,i));
                             end
                             % SNPs
                             if ~isempty(snpind)
                                 if isempty(ggind)
                                    wf(snpind,i) = BIOMETRICSMV.sum2one(wf(snpind,i)); 
                                    wfh(snpind,i) = BIOMETRICSMV.sum2one(wfh(snpind,i));
                                 else
                                    gg = unique(ggind); 
                                    ngg = length(gg);
                                    for j=1:1:ngg
                                        ind = find(ggind==j);
                                        avwf = nanmean(wf(snpind(ind),i));
                                        avwfh = nanmean(wfh(snpind(ind),i));
                                        wf(snpind(ind),i) = avwf*BIOMETRICSMV.sum2one(wf(snpind(ind),i)); 
                                        wfh(snpind(ind),i) = avwfh*BIOMETRICSMV.sum2one(wfh(snpind(ind),i));
                                    end
                                    wf(snpind,i) = BIOMETRICSMV.sum2one(wf(snpind,i)); 
                                    wfh(snpind,i) = BIOMETRICSMV.sum2one(wfh(snpind,i));
                                 end
                             end
                         end    
                     case false
                         % do nothing;
                         if obj.SNPGroupNorm
                             snpind = obj.SNPInd;
                             ggind = obj.SNPGroups;
                             for i=1:1:nrT
                                 gg = unique(ggind); 
                                 ngg = length(gg);
                                 for j=1:1:ngg
                                     ind = find(ggind==j);
                                     avwf = nanmean(wf(snpind(ind),i));
                                     avwfh = nanmean(wfh(snpind(ind),i));
                                     wf(snpind(ind),i) = avwf*BIOMETRICSMV.sum2one(wf(snpind(ind),i)); 
                                     wfh(snpind(ind),i) = avwfh*BIOMETRICSMV.sum2one(wfh(snpind(ind),i));
                                 end
                             end
                         end
                 end
                 switch obj.Classify
                   case 'SOFT'
                       TScores = obj.SOFTTMatches;
                       FScores = obj.SOFTFMatches;
                   case 'HARD'
                       TScores = obj.HARDTMatches;
                       FScores = obj.HARDFMatches;
                 end
                 for i=1:1:nG 
                     % STILL TO TAKE CARE OF TESTCASES CONTAINING ONLY NAN
                     % VALUES!!!
                     r=r+1;ind = index{i};
                     tmpwf = wf(ind,:);
%                      for c=1:1:size(tmpwf,2)
%                         tmpwf(:,c) = BIOMETRICSMV.sum2one(tmpwf(:,c)); 
%                      end
%                      for c=1:1:size(tmpwf,1)
%                         tmpwf(c,:) = BIOMETRICSMV.sum2one(tmpwf(c,:)); 
%                      end

                     [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:),AUC(:,r),~] = BIOMETRICSMV.getBiometrics(TScores(ind,:),FScores(ind,:,:),tmpwf);
                     r=r+1;
                     tmpwfh = wfh(ind,:);
%                      for c=1:1:size(tmpwf,2)
%                         tmpwfh(:,c) = BIOMETRICSMV.sum2one(tmpwfh(:,c)); 
%                      end
%                      for c=1:1:size(tmpwf,1)
%                         tmpwfh(c,:) = BIOMETRICSMV.sum2one(tmpwfh(c,:)); 
%                      end
                     [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICSMV.getBiometrics(TScores(ind,:),FScores(ind,:,:),tmpwfh);
                     [FM{i}.A,FM{i}.evalA,FM{i}.D,FM{i}.evalD,FM{i}.X] = BIOMETRICSMV.evalSortedResults(ST,SF,L1CL1TestDepVar,T);
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
                        if mod(i,2)==1, continue; end
                        plot(1:100,squeeze(CumRanks(1,i,:)),str{i},'LineWidth',2);
                    end
                    figure;hold on;title('Verification');
                    xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
                    plot(0:0.01:1,0:0.01:1,'k--');
                    plot(0:0.01:1,1:-0.01:0,'k-');
                    for i=1:1:nr
                        if mod(i,2)==1, continue; end
                        plot(X(i,:),Y(i,:),str{i},'LineWidth',2);
                    end
                    meanA = mean(FM{1}.A);
                    meanD = mean(FM{1}.D);
                    figure;set(gca,'ylim',[meanA-0.06 meanA+0.06],'xlim',[0 100]);hold on;grid on;
                    plot(((1:nrT-1)/(nrT-1))*100,meanA*ones(1,nrT-1),'k-');
                    title('False Matches Angle');
                    counter = 0;
                    for i=1:1:nG
                        counter = counter+2;
                        plot(FM{i}.X,FM{i}.evalA,str{counter},'lineWidth',2);
                    end
                    figure;set(gca,'ylim',[meanA-0.06 meanA+0.06],'xlim',[0 100]);hold on;grid on;
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
        function out = biometricPerformancev2(obj,TestDepVar,uindex,vis,index)
                 if nargin<4,index = {1:obj.nPPM};end
                 if isempty(index),index = {1:obj.nPPM};end
                 if ~iscell(index), index = {index};end
                 L1CL1TestDepVar = BIOMETRICSMV.getL1CL1DepVar(TestDepVar);
                 nrT = size(L1CL1TestDepVar,1);
                 nG = length(index);
                 nr = nG*2;
                 nrT2 = length(uindex);
                 EER = nan*zeros(1,nr);Ranks = nan*zeros(1,nr,length(uindex));AUC = nan*zeros(1,nr);
                 X = nan*zeros(nr,(nrT*nrT2)-1);Y = nan*zeros(nr,(nrT*nrT2)-1);
                 FM = cell(1,nG);
                 r = 0;
                 switch obj.HitWeights
                     case true
                         %wh = obj.wH;
                         wh = obj.wL;
                     case false
                         wh = ones(size(obj.wH));
                 end
                 switch obj.FreqWeights
                     case true
                         wf = obj.wF;
                     case false
                         wf = ones(size(obj.wF));
                         wf(isnan(obj.wF)) = nan;
                 end
                 wfh = repmat(wh,1,nrT).*wf;% combined frequency and hits weights
                 switch obj.GroupNorm
                     case true
                         sexind = obj.SEXInd;
                         ageind = obj.AGEInd;
                         hind = obj.HInd;
                         wind = obj.WInd;
                         gbind = obj.GBInd;
                         snpind = obj.SNPInd;
                         ggind = obj.SNPGroups;
                         for i=1:1:nrT
                             % SEX
                             if ~isempty(sexind)
                                wf(sexind,i) = BIOMETRICSMV.sum2one(wf(sexind,i)); 
                                wfh(sexind,i) = BIOMETRICSMV.sum2one(wfh(sexind,i)); 
                             end
                             % AGE
                             if ~isempty(ageind)
                                 wf(ageind,i) = BIOMETRICSMV.sum2one(wf(ageind,i)); 
                                 wfh(ageind,i) = BIOMETRICSMV.sum2one(wfh(ageind,i));
                             end
                             % HEIGHT 
                             if ~isempty(hind)
                                 wf(hind,i) = BIOMETRICSMV.sum2one(wf(hind,i)); 
                                 wfh(hind,i) = BIOMETRICSMV.sum2one(wfh(hind,i));
                             end
                             % WEIGHT 
                             if ~isempty(wind)
                                 wf(wind,i) = BIOMETRICSMV.sum2one(wf(wind,i)); 
                                 wfh(wind,i) = BIOMETRICSMV.sum2one(wfh(wind,i));
                             end
                             % GBBACKGROUND
                             if ~isempty(gbind)
                                 wf(gbind,i) = BIOMETRICSMV.sum2one(wf(gbind,i)); 
                                 wfh(gbind,i) = BIOMETRICSMV.sum2one(wfh(gbind,i));
                             end
                             % SNPs
                             if ~isempty(snpind)
                                 if isempty(ggind)
                                    wf(snpind,i) = BIOMETRICSMV.sum2one(wf(snpind,i)); 
                                    wfh(snpind,i) = BIOMETRICSMV.sum2one(wfh(snpind,i));
                                 else
                                    gg = unique(ggind); 
                                    ngg = length(gg);
                                    for j=1:1:ngg
                                        ind = find(ggind==j);
                                        avwf = nanmean(wf(snpind(ind),i));
                                        avwfh = nanmean(wfh(snpind(ind),i));
                                        wf(snpind(ind),i) = avwf*BIOMETRICSMV.sum2one(wf(snpind(ind),i)); 
                                        wfh(snpind(ind),i) = avwfh*BIOMETRICSMV.sum2one(wfh(snpind(ind),i));
                                    end
                                    wf(snpind,i) = BIOMETRICSMV.sum2one(wf(snpind,i)); 
                                    wfh(snpind,i) = BIOMETRICSMV.sum2one(wfh(snpind,i));
                                 end
                             end
                         end    
                     case false
                         % do nothing;
                         if obj.SNPGroupNorm
                             snpind = obj.SNPInd;
                             ggind = obj.SNPGroups;
                             for i=1:1:nrT
                                 gg = unique(ggind); 
                                 ngg = length(gg);
                                 for j=1:1:ngg
                                     ind = find(ggind==j);
                                     avwf = nanmean(wf(snpind(ind),i));
                                     avwfh = nanmean(wfh(snpind(ind),i));
                                     wf(snpind(ind),i) = avwf*BIOMETRICSMV.sum2one(wf(snpind(ind),i)); 
                                     wfh(snpind(ind),i) = avwfh*BIOMETRICSMV.sum2one(wfh(snpind(ind),i));
                                 end
                             end
                         end
                 end
                 switch obj.Classify
                   case 'SOFT'
                       TScores = obj.SOFTTMatches;
                       FScores = obj.SOFTFMatches;
                   case 'HARD'
                       TScores = obj.HARDTMatches;
                       FScores = obj.HARDFMatches;
                 end
                 TScores = TScores(:,uindex);
                 FScores = FScores(:,uindex,:);
                 wf = wf(:,uindex);
                 wfh = wfh(:,uindex);
                 for i=1:1:nG 
                     % STILL TO TAKE CARE OF TESTCASES CONTAINING ONLY NAN
                     % VALUES!!!
                     r=r+1;ind = index{i};
                     tmpwf = wf(ind,:);
                     [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:),AUC(:,r),~] = BIOMETRICSMV.getBiometrics(TScores(ind,:),FScores(ind,:,:),tmpwf);
                     r=r+1;
                     tmpwfh = wfh(ind,:);
                     [Ranks(:,r,:),EER(:,r),ST,SF,X(r,:),Y(r,:),AUC(:,r),T] = BIOMETRICSMV.getBiometrics(TScores(ind,:),FScores(ind,:,:),tmpwfh);
                     [FM{i}.A,FM{i}.evalA,FM{i}.D,FM{i}.evalD,FM{i}.X] = BIOMETRICSMV.evalSortedResultsv2(ST,SF,L1CL1TestDepVar,T,uindex);
                 end
                 Ranks(1,:,:) = (Ranks(1,:,:)./nrT).*100;
                 CumRanks = zeros(1,nr,100);
                 for cr=1:1:100
                     CumRanks(:,:,cr) = (sum(Ranks<=cr,3)/nrT2).*100;
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
                        if mod(i,2)==1, continue; end
                        plot(1:100,squeeze(CumRanks(1,i,:)),str{i},'LineWidth',2);
                    end
                    figure;hold on;title('Verification');
                    xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
                    plot(0:0.01:1,0:0.01:1,'k--');
                    plot(0:0.01:1,1:-0.01:0,'k-');
                    for i=1:1:nr
                        if mod(i,2)==1, continue; end
                        plot(X(i,:),Y(i,:),str{i},'LineWidth',2);
                    end
                    meanA = mean(FM{1}.A);
                    meanD = mean(FM{1}.D);
                    figure;set(gca,'ylim',[meanA-0.08 meanA+0.08],'xlim',[0 100]);hold on;grid on;
                    plot(((1:nrT-1)/(nrT-1))*100,meanA*ones(1,nrT-1),'k-');
                    title('False Matches Angle');
                    counter = 0;
                    for i=1:1:nG
                        counter = counter+2;
                        plot(FM{i}.X,FM{i}.evalA,str{counter},'lineWidth',2);
                    end
                    figure;set(gca,'ylim',[meanA-0.08 meanA+0.08],'xlim',[0 100]);hold on;grid on;
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
    methods % EVOMORPH
        function out = matchSolutions(obj,X,DepVar,index)
                 L1CL1DepVar = BIOMETRICSMV.getL1CL1DepVar(DepVar);
                 nrT = size(L1CL1DepVar,1);
                 SM = nan*zeros(obj.nPPM,nrT);
                 HM = nan*zeros(obj.nPPM,nrT);
                 wF = nan*zeros(obj.nPPM,1); %#ok<*PROPLC>
                 wH = nan*zeros(obj.nPPM,1);
                 wL = nan*zeros(obj.nPPM,1);
                 f = waitbar(0,'MATCHING');
                 for i=1:length(index)
                    [SM(index(i),:),HM(index(i),:),wF(index(i)),wH(index(i)),wL(index(i))] = matchSolutions(obj.PPM{index(i)},X(index(i)),DepVar,nrT); 
                    waitbar(i/length(index),f);
                 end
                 delete(f);
                 switch obj.HitWeights
                     case true
                         %wh = wH;
                         wh = wL;
                     case false
                         wh = ones(size(wH));
                 end
                 switch obj.FreqWeights
                     case true
                         wf = wF;
                     case false
                         wf = ones(size(wF));
                         wf(isnan(wF)) = nan;
                 end
                 wfh = wh.*wf;% combined frequency and hits weights
                 switch obj.GroupNorm
                     case true
                         sexind = obj.SEXInd;
                         ageind = obj.AGEInd;
                         hind = obj.HInd;
                         wind = obj.WInd;
                         gbind = obj.GBInd;
                         snpind = obj.SNPInd;
                         ggind = obj.SNPGroups;
                         for i=1:1:1
                             % SEX
                             if ~isempty(sexind)
                                wf(sexind,i) = BIOMETRICSMV.sum2one(wf(sexind,i)); 
                                wfh(sexind,i) = BIOMETRICSMV.sum2one(wfh(sexind,i)); 
                             end
                             % AGE
                             if ~isempty(ageind)
                                 wf(ageind,i) = BIOMETRICSMV.sum2one(wf(ageind,i)); 
                                 wfh(ageind,i) = BIOMETRICSMV.sum2one(wfh(ageind,i));
                             end
                             % HEIGHT 
                             if ~isempty(hind)
                                 wf(hind,i) = BIOMETRICSMV.sum2one(wf(hind,i)); 
                                 wfh(hind,i) = BIOMETRICSMV.sum2one(wfh(hind,i));
                             end
                             % WEIGHT 
                             if ~isempty(wind)
                                 wf(wind,i) = BIOMETRICSMV.sum2one(wf(wind,i)); 
                                 wfh(wind,i) = BIOMETRICSMV.sum2one(wfh(wind,i));
                             end
                             % GBBACKGROUND
                             if ~isempty(gbind)
                                 wf(gbind,i) = BIOMETRICSMV.sum2one(wf(gbind,i)); 
                                 wfh(gbind,i) = BIOMETRICSMV.sum2one(wfh(gbind,i));
                             end
                             % SNPs
                             if ~isempty(snpind)
                                 if isempty(ggind)
                                    wf(snpind,i) = BIOMETRICSMV.sum2one(wf(snpind,i)); 
                                    wfh(snpind,i) = BIOMETRICSMV.sum2one(wfh(snpind,i));
                                 else
                                    gg = unique(ggind); 
                                    ngg = length(gg);
                                    for j=1:1:ngg
                                        ind = find(ggind==j);
                                        avwf = nanmean(wf(snpind(ind),i));
                                        avwfh = nanmean(wfh(snpind(ind),i));
                                        wf(snpind(ind),i) = avwf*BIOMETRICSMV.sum2one(wf(snpind(ind),i)); 
                                        wfh(snpind(ind),i) = avwfh*BIOMETRICSMV.sum2one(wfh(snpind(ind),i));
                                    end
                                    wf(snpind,i) = BIOMETRICSMV.sum2one(wf(snpind,i)); 
                                    wfh(snpind,i) = BIOMETRICSMV.sum2one(wfh(snpind,i));
                                 end
                             end
                         end    
                     case false
                         % do nothing;
                         if obj.SNPGroupNorm
                             snpind = obj.SNPInd;
                             ggind = obj.SNPGroups;
                             for i=1:1:1
                                 gg = unique(ggind); 
                                 ngg = length(gg);
                                 for j=1:1:ngg
                                     ind = find(ggind==j);
                                     avwf = nanmean(wf(snpind(ind),i));
                                     avwfh = nanmean(wfh(snpind(ind),i));
                                     wf(snpind(ind),i) = avwf*BIOMETRICSMV.sum2one(wf(snpind(ind),i)); 
                                     wfh(snpind(ind),i) = avwfh*BIOMETRICSMV.sum2one(wfh(snpind(ind),i));
                                 end
                             end
                         end
                 end
                 wfh = BIOMETRICSMV.sum2one(wfh);
                 switch obj.Classify
                   case 'SOFT'
                       out = repmat(wfh,1,nrT).*SM.*HM;
                   case 'HARD'
                       out = repmat(wfh,1,nrT).*HM.*SM;
                 end
                 out = nansum(out,1);
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
        function [R,EER,ST,SF,X,Y,AUC,T] = getBiometrics(TScores,FScores,wM)
                 nF = size(FScores,3);
                 nT = size(TScores,2);
                 %tmpW = repmat(W,1,nT).*wM;
                 ST = BIOMETRICSMV.accScores(TScores,wM)';
                 SF = BIOMETRICSMV.accScores(FScores,repmat(wM,1,1,nF));
                 [R,EER,X,Y,AUC,T] = BIOMETRICSMV.evalScores(ST,SF,nT,nF); 
        end
        function Score = accScores(Scores,tmpW)
                Score = squeeze(nansum(tmpW.*Scores,1)./nansum(tmpW,1));
        end
        function [R,EER,x,y,AUC,T] = evalScores(ST,SF,nT,nF)
                 % Rank analysis
                 R = (sum(SF<=repmat(ST,1,nF),2)+1);%/(nF+1);
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
        function [A,evalA,D,evalD,X] = evalSortedResultsv2(ST,SF,DepVar,T,uindex)
            
                %[FM{i}.A,FM{i}.evalA,FM{i}.D,FM{i}.evalD,FM{i}.X] = BIOMETRICSMV.evalSortedResults(ST,SF,L1CL1TestDepVar,T);
                %DepVar = L1CL1TestDepVar;
                % DepVar = TestDepVar;
                nTot = size(DepVar,1);
                %nrT = length(uindex);
                [nrT,nrF] = size(SF);
                D = nan*zeros(nrT,nrF);
                A = nan*zeros(nrT,nrF);
                AT = nan*zeros(nrT,2);
                for t=1:1:nrT
                    %t=1;
                    tind = setdiff(1:nTot,uindex(t));
                    forD = pdist2(DepVar(uindex(t),:),DepVar(tind,:),'euclidean');
                    %forA = -1*(pdist2(DepVar(t,:),DepVar(tind,:),'cosine')-1);
                    forA = (pdist2(DepVar(uindex(t),:),DepVar(tind,:),'cosine'));
                    [~,ind] = sort(SF(t,:),'ascend');
                    D(t,:) = forD(ind);
                    A(t,:) = forA(ind);
                    AT(t,1) = mean(forA(SF(t,:)<=T));
                    AT(t,2) = mean(forA(SF(t,:)>=T));   
                end
                D = mean(D,1);
                A = mean(A,1);
                X = ((1:nrF)/(nrF-1))*100;
                polyA = polyfit(1:nrF,A,3);
                evalA = polyval(polyA,1:nrF);
                polyD = polyfit(1:nrF,D,3);
                evalD = polyval(polyD,1:nrF);
                
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
        function out = selectDepVar(DepVar,index)
                 out = cell(size(DepVar));   
                 for i=1:1:length(DepVar);
                     out{i}.DepVar = cell(size(DepVar{i}.DepVar));
                     for j=1:1:length(DepVar{i}.DepVar)
                        tmp = DepVar{i}.DepVar{j};
                        out{i}.DepVar{j} = tmp(index,:); 
                     end
                 end 
        end
        function out = sum2one(v)
            out = v/nansum(v);
        end
    end
end