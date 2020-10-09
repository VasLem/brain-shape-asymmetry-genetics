classdef CONTMV2 < CATMV2
    properties
        TargX;
        RangeX;
    end
    properties (Dependent = true)
    end
    properties (Hidden = true) 
        %ON = 100;
        KN = 100;
        MinSample = 10;
        Bins = 10;
        CORT = 0.1;
        SeqX = {};
        SeqTargX;
        nSeqX = 20;
        UseSeq = true;
        MinBal = 0.1;
    end
    properties(Hidden = true, Dependent = true)
        XREG;
        PosClass;
        ON;
    end
    methods % CONSTRUCTOR
        function obj = CONTMV2(varargin)
            obj = obj@CATMV2(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.XREG(obj)
            out = obj.X;
        end
        function out = get.PosClass(obj) %#ok<*MANU>
           out = 1; 
        end
        function out = get.ON(obj)
           out = obj.KN; 
        end
    end
    methods % PREDICTOR
        function trainPredictor(obj,X,TargetX,FS,K,t)
%            % PREPARATION
%              BKX = X;
%              [X,COST,W,ind,~,~] = prepPredictorTraining(obj,X,TargetX);
%              %W = ones(size(X));
%            % FEATURE SELECTION
%              featureSelection(obj,X(ind),BKX(ind),FS(ind,:));
%            % ACCUMULATE FEATURES
%              sFS = selectFeatures(obj,FS);
%            % TRAINING
%              fitClassifier(obj,sFS(ind,:),X(ind),COST,W(ind));
%             % LEARNING HOW TO NORMALIZE  
%              [~,PXC] = predict(obj,FS);
%              learnPXCNorm(obj,PXC);
              if nargin<6, t = 0; end
              switch obj.ClassifierCascade
                case false
                    trainSVMRBFSingle(obj,X,TargetX,FS,K,t);
                case true
                    trainSVMRBFCascade(obj,X,TargetX,FS,K,t);
              end
        end
        function trainSVMRBFSingle(obj,X,TargetX,FS,K,t)
             if nargin<5, K = 5;end
             if nargin<6, t = 0; end
           % PREPARATION
             obj.ClassifierType = 'SVMRBF';
             BKX = X;
             [X,COST,W,ind,~,~] = prepPredictorTraining(obj,BKX,TargetX);
             %W = ones(size(X));
             obj.PosClassBal = getPosClassBal(obj,X(ind));
             obj.nTraining = length(ind);
             obj.FSHits = obj.Hits;
           % FEATURE SELECTION
             FS = getFS(obj,FS);
             featureSelection(obj,X(ind),BKX(ind),FS(ind,:));
           % SVM LEARNING  
             if ~isobject(K), rng(1);K = cvpartition(length(ind),'KFold',K);end
             out = CATMV2.trainSVMRBF(FS(ind,obj.FSInd),X(ind),COST,W(ind),K,obj.OptSVMRBF,obj.PosClass,t);
             obj.Classifier = out.Classifier;
             obj.CVloss = out.cv;
             obj.pCV = out.pCV;
             obj.CVEER = out.EER;
             obj.CVG = out.G;
             obj.CVAUC = out.AUC;
             obj.OffsetPXC = out.OffsetPXC;
             obj.ScalePXC = out.ScalePXC;              
        end
        function trainSVMRBFCascade(obj,X,TargetX,FS,K,t)
             if nargin<5, K = 5;end
             if nargin<6, t = 0; end
           % PREPARATION
             obj.ClassifierType = 'SVMRBF';
             [X,COST,W,ind,~,~] = prepPredictorTraining(obj,X,TargetX);
             obj.PosClassBal = getPosClassBal(obj,X(ind));
             obj.nTraining = length(ind);
           % MODULE SELECTION
             obj.CascadeInd = find(obj.Hits);
             classifier = cell(1,obj.nCascade);
             cvloss = zeros(1,obj.nCascade);
             pcv = zeros(1,obj.nCascade);
             cveer = zeros(1,obj.nCascade);
             cvg = zeros(1,obj.nCascade);
             cvauc = zeros(1,obj.nCascade);
             offset = zeros(2,obj.nCascade);
             scale = zeros(2,obj.nCascade);
             if ~isobject(K),rng(1);K = cvpartition(length(ind),'KFold',K);end
             parfor i=1:obj.nCascade
                    forFS = getDepVar(obj,FS,obj.CascadeInd(i));
                    out = CATMV2.trainSVMRBF(forFS(ind,:),X(ind),COST,W(ind),K,obj.OptSVMRBF,obj.PosClass,t);
                    classifier{i} = out.Classifier;
                    cvloss(i) = out.cv;
                    cveer(i) = out.EER;
                    cvg(i) = out.G;
                    cvauc(i) = out.AUC;
                    offset(:,i) = out.OffsetPXC';
                    scale(:,i) = out.ScalePXC';
                    pcv(i) = out.pCV;
             end
             obj.Classifier = classifier;
             obj.CVloss = ones(1,obj.nLC);
             obj.CVloss(obj.CascadeInd) = cvloss;
             obj.CasCVloss = obj.CVloss;
             obj.pCV = ones(1,obj.nLC);
             obj.pCV(obj.CascadeInd) = pcv;
             obj.CVEER = ones(1,obj.nLC);
             obj.CVEER(obj.CascadeInd) = cveer;
             obj.CasCVEER = obj.CVEER;
             obj.CVG = cvg;
             obj.CVAUC = cvauc;
             obj.OffsetPXC = offset;
             obj.ScalePXC = scale;
        end
        function out = testPredictor(obj,X,FS,TrX,TrFS,display)
            % PREPARATION
              [X,ind] = prepPredictorTesting(obj,X);
              X = X(ind);
              FS = FS(ind,:);
            % OPPOSITE AXIS TESTING
              nT = length(ind); 
              TMatches = nan*zeros(nT,1);
              FMatches = nan*zeros(nT,obj.ON);
              TClass = nan*zeros(nT,1);
              FClass = nan*zeros(nT,obj.ON);
              nrF = obj.ON;
              [path,ID] = setupParForProgress(nT);
              parfor t=1:nT
                 diff = abs(X-X(t));
                 [~,Find] = sort(diff);
                 Find = Find(end-nrF+1:end);
                 forobj = clone(obj);
                 trainPredictor(forobj,TrX,X(t),TrFS);
                 [PX,PXC] = predict(forobj,FS);
                 [MS,MH] = Pred2Match(forobj,1,PX,PXC);
                 TMatches(t) = MS(t);
                 FMatches(t,:) = MS(Find);
                 TClass(t) = MH(t);
                 FClass(t,:) = MH(Find);
                 parfor_progress;
              end
              closeParForProgress(path,ID);
              [out.R1,out.R10,out.R20,out.CR] =  PPMMV2.getRANKSTAT(TMatches,FMatches);
              [out.EER,out.G,out.AUC,out.XX,out.YY,out.TH,out.Y] = PPMMV2.getEER(TMatches,FMatches(:));% CHECK ORDER
              [out.HPREC,out.HREC,out.HG,out.HTNF,out.HX,out.HY,out.HACC] = PPMMV2.getHardClass(1-TClass,1-FClass(:));          
            % VISUALIZATION
              if display
                 figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                 plot(out.XX,out.YY,'b-','LineWidth',1.5);
                 plot(out.HX,out.HY,'r.','MarkerSize',20);
                 title(['EER: ' num2str(out.EER) ' HG: ' num2str(out.HG)]);
                 figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
                 xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
                 title('Identification');
                 plot(1:100,1:100,'k--');
                 plot(1:100,out.CR,'g-','LineWidth',1.5);
              end
        end
        function trainPredictorSequence(obj,X,FS,K,t)
              if nargin < 5, t=0;end
              obj.SeqX = {};
              targets = min(X):((max(X)-min(X))/(obj.nSeqX-1)):max(X);
              obj.SeqTargX = targets;
              obj.Bins = round(obj.nSeqX/2);
              n = length(targets);
              tmp = cell(1,n);
              [path,ID] = setupParForProgress(n);
              parfor p=1:n
                 forobj = clone(obj);
                 trainPredictor(forobj,X,targets(p),FS,K,t);
                 tmp{p} = forobj;
                 parfor_progress;
              end
              closeParForProgress(path,ID);
              obj.SeqX = tmp;
        end
        function trainPredictorSequenceStaged(obj,X,FS,K,par,pcv,t)
              if nargin < 7, t=0;end
              if nargin < 6, pcv = 0.35; end
              if nargin < 5, par = 0.1; end
              obj.SeqX = {};
              targets = min(X):((max(X)-min(X))/(obj.nSeqX-1)):max(X);
              %obj.SeqTargX = targets;
              n = length(targets);
              tmp = cell(1,n);
              [path,ID] = setupParForProgress(n);
              good = zeros(1,n);
              parfor p=1:n
                 try
                     forobj = clone(obj);
                     forobj.pCrit = par;
                     forobj.HitType = 'pAR';
                     if forobj.nHits==0, continue; end
                     forobj.ClassifierCascade = true;
                     forobj.FSelection = 'NONE';
                     forobj.OptSVMRBF = true;
                     forobj.maxFS = 1000000; 
                     trainPredictor(forobj,X,targets(p),FS,K,t);
                     forobj.pCrit = pcv;
                     forobj.HitType = 'CasCVloss';
                     if forobj.nHits==0, continue; end
                     forobj.ClassifierCascade = false;
                     trainPredictor(forobj,X,targets(p),FS,K,t);
                     tmp{p} = forobj;
                     good(p) = 1;
                 catch
                     disp([num2str(p) '_' num2str(targets(p)) '_FAILED']);
                 end
                 parfor_progress;
              end
              closeParForProgress(path,ID);
              index = find(good);
              obj.SeqTargX = targets(index);
              obj.SeqX = tmp(index);
        end
        function out = testPredictorSequence(obj,X,FS,display)
            % PREPARATION
              [X,ind] = prepPredictorTesting(obj,X);
              X = X(ind);
              FS = PPMMV2.selectDepVar(FS,ind);%FS = FS(ind,:);
            % OPPOSITE AXIS TESTING
              nT = length(ind); 
              nrF = obj.ON;
              if nrF>=(nT-1), nrF = nT-1; end
              TMatches = nan*zeros(nT,1);
              FMatches = nan*zeros(nT,nrF);
              TClass = nan*zeros(nT,1);
              FClass = nan*zeros(nT,nrF);
              %[path,ID] = setupParForProgress(nT);
              for t=1:nT
                 diff = abs(X-X(t));
                 [~,Find] = sort(diff);
                 Find = Find(end-nrF+1:end);
                 [MS,MH] = match2Seq(obj,X(t),FS,nT);
                 TMatches(t) = MS(t);
                 FMatches(t,:) = MS(Find);
                 TClass(t) = MH(t);
                 FClass(t,:) = MH(Find);
                 %parfor_progress;
              end
              %closeParForProgress(path,ID);
              [out.R1,out.R10,out.R20,out.CR] =  PPMMV2.getRANKSTAT(TMatches,FMatches);
              [out.EER,out.G,out.AUC,out.XX,out.YY,out.TH,out.Y] = PPMMV2.getEER(TMatches,FMatches(:));% CHECK ORDER
              [out.HPREC,out.HREC,out.HG,out.HTNF,out.HX,out.HY,out.HACC] = PPMMV2.getHardClass(1-TClass,1-FClass(:));          
            % VISUALIZATION
              if display
                 figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                 plot(out.XX,out.YY,'b-','LineWidth',1.5);
                 plot(out.HX,out.HY,'r.','MarkerSize',20);
                 title(['EER: ' num2str(out.EER) ' HG: ' num2str(out.HG)]);
                 figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
                 xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
                 title('Identification');
                 plot(1:100,1:100,'k--');
                 plot(1:100,out.CR,'g-','LineWidth',1.5);
              end
        end
        function [MS,MH,ind] = match2Seq(obj,TargX,FS,nT)
            if isempty(obj.SeqX), error('Sequence not trained'); end      
            diff = abs(obj.SeqTargX-TargX);
            [~,ind] = min(diff);
            [PX,~,PXCN] = predict(obj.SeqX{ind},FS,nT);
            [MS,MH] = Pred2Match(obj.SeqX{ind},1,PX,PXCN);
        end
    end
    methods % BIOMETRICS
       function out = prepBiometricTest(obj,X)
           out = X; 
       end
       function [SOFTTMatches,SOFTFMatches,HARDTMatches,HARDFMatches,wF,wH,wL] = biometricMatchesWeights(obj,TestDepVar,TestX)
           % INITIALIZE 
            nrT = size(TestX,1);
            TestX = prepBiometricTest(obj,TestX);
           % WEIGHTS 
            wH = obj.nHits;
            wL = 2; % TO BE CORRECTED
           % FEATURES 
            %RIPF = getRIP(obj,TestDepVar);
            %PCF = getDepVar(obj,TestDepVar,0);
            %FS = [RIPF, PCF];
            FS = TestDepVar;
           % PREDICTIONS
            %[PX,PXC] = predict(obj,FS);
           % MATCHES 
            SOFTTMatches = nan*zeros(1,nrT);
            SOFTFMatches = nan*zeros(1,nrT,nrT-1);
            HARDTMatches = nan*zeros(1,nrT);
            HARDFMatches = nan*zeros(1,nrT,nrT-1);
            wF = nan*zeros(1,nrT);
            for t=1:nrT
                % t=2;
                tInd = setdiff(1:nrT,t);
                wF(t) = 1-getFreq(obj,TestX(t));
                if isnan(wF(t)), continue; end
                [SM,HM,sind] = match2Seq(obj,TestX(t),FS,nrT);
                wF(t) = obj.SeqX{sind}.CVloss;
                wF(t) = -log(wF(t)./(1-wF(t)));
                wF(t) = 1;
                if wF(t)<0, wF(t) = 0; end
                %wF(t) = wF(t)*obj.SeqX{sind}.PosClassBal;
                SOFTTMatches(:,t) = SM(t);
                SOFTFMatches(:,t,:) = SM(tInd);
                HARDTMatches(:,t) = HM(t);
                HARDFMatches(:,t,:) = HM(tInd);
            end
        end
    end
    methods % EVOMORPH
        function out = prepEvoMorphTest(obj,X)
           out = X; 
        end
        function [SM,HM,wF,wH,wL] = matchSolutions(obj,TestX,FS,nrT)
            TestX = prepEvoMorphTest(obj,TestX);
            if isnan(TestX), SM = nan*zeros(nrT,1);HM = nan*zeros(nrT,1);wF = nan;wH = nan;wL = nan;return;end
            switch obj.UseSeq
                case true
                    [SM,HM,ind] = match2Seq(obj,TestX,FS,nrT);
                    wL = max(-log(obj.SeqX{ind}.CVloss./(1-obj.SeqX{ind}.CVloss)),0);
                    %wL = wL*obj.PosClassBal;
                case false
                    if ~(TestX==obj.TargetX), error('test and target value CONTMV do not match'); end
                    [PX,~,PXCN] = predict(obj,FS);
                    [SM,HM] = Pred2Match(obj,1,PX,PXCN);
                    wL = max(-log(obj.CVloss./(1-obj.CVloss)),0);
                    %wL = wL*obj.PosClassBal;
            end
            wF = 1-getFreq(obj,TestX);
            wH = obj.nHits;
        end
    end
    methods % INTERFACING
        function [NX,COST,W,ind,Tind,Find] = prepPredictorTraining(obj,X,TargetX)
            %X = BKX;
            
            obj.TargX = TargetX;
            NX = nan*zeros(size(X));
            W = nan*zeros(size(X));
            ind = find(~isnan(X));
            
            range = max(X)-min(X);
            range = range/obj.Bins;
           % select closest examples
            dif = abs(TargetX-X(ind));
            Tind2 = find(dif<=range);
            bal = length(Tind2)/length(ind);
            while bal<obj.MinBal
                range = range*1.5;
                Tind2 = find(dif<=range);
                bal = length(Tind2)/length(ind);
            end
            obj.RangeX = range;
              
%             IB = CONTMV2.getTrOutliers(dif',0.5);
            %index = find(IB>=0.5);
            Tind = ind(Tind2);
            NX(Tind) = obj.PosClass;
            tmp = dif(Tind2);
            if max(tmp)==0
               W(Tind) = 1;
            else
               tmp = tmp-min(tmp);tmp = tmp/max(tmp);tmp = 1-tmp;
               %mval = min(tmp(tmp>0));tmp(tmp==0) = mval;
               W(Tind) = tmp;
            end
            W(Tind) = CONTMV2.sum2one(W(Tind));
           % opposite examples
%             IB = CONTMV2.getTrOutliers(dif',1);
%             Find2 = find(IB<0.5);
            Find2 = find(dif>range);
            Find = ind(Find2);
            NX(Find) = -1*obj.PosClass;
            tmp = dif(Find2);
            tmp = tmp-min(tmp);tmp = tmp/max(tmp);
            %mval = min(tmp(tmp>0));tmp(tmp==0) = mval;
            W(Find) = tmp;
            W(Find) = CONTMV2.sum2one(W(Find));
            %ind = union(Tind,Find);
            W(Tind) = W(Tind)*(1-(length(Tind)/length(ind)));
            W(Find) = W(Find)*(1-(length(Find)/length(ind)));
            %W(ind) = W(ind)+0.01;
            ind = find(W>0);     
            %W(ind) = W(ind)+0.000000000001;
           % Defining COST 
            COST.ClassNames = [obj.PosClass -1*obj.PosClass];
            COST.ClassificationCosts = zeros(2,2);
            %COST.ClassificationCosts(1,2) = 1;
            %COST.ClassificationCosts(2,1) = 1;
            COST.ClassificationCosts(1,2) = 1-(length(Tind)/length(ind));
            COST.ClassificationCosts(2,1) = 1-(length(Find)/length(ind));
            %W = nan*zeros(size(X));
            %W(Tind) = 1-(length(Tind)/length(ind));
            %W(Find) = 1-(length(Find)/length(ind));
        end
        function [X,ind] = prepPredictorTesting(obj,X)
            ind = find(~isnan(X));
        end
        function out = getVarInfo(obj,DepVar)
            [sortv,ind] = sort(obj.X(~isnan(obj.X)),'ascend');
           % lower 25%
            indl = round(0.25*length(ind));
            out.Mb = sortv(indl);
            Mindex = find(obj.X<=out.Mb);
            out.Mel = mean(obj.X(Mindex));
            out.Mrange = [min(obj.X(Mindex)) max(obj.X(Mindex))];
            out.MDepVar = mean(DepVar(Mindex,:));
           % upper 25%
            indu = round(0.75*length(ind));
            out.Pb = sortv(indu);
            Pindex = find(obj.X>=out.Pb);
            out.Pel = mean(obj.X(Pindex));
            out.Prange = [min(obj.X(Pindex)) max(obj.X(Pindex))];
            out.PDepVar = mean(DepVar(Pindex,:));
            out.Range = out.Pel-out.Mel;
        end
        function out = getFreq(obj,in) %#ok<*INUSL>
           if isnan(in), out = nan; return; end
           out = 0;
           %out = nan*zeros(size(in));
           %ind =  find(~isnan(in));
           %if isempty(in), return; end
           %[N,x] = hist(obj.XMOD,min(min(obj.XMOD),min(in)):(max(obj.XMOD)-min(obj.XMOD))/obj.Bins:max(max(obj.XMOD),max(in)));
           %out(ind) = interp1(x,N,in(ind))/length(obj.XMOD);
        end
        function out = getW(obj,in)
           out = nan*zeros(size(in));
           ind =  find(~isnan(in));
           if isempty(in), return; end
           [N,x] = hist(in(ind),obj.Bins);
           out(ind) = 1-interp1(x,N,in(ind))/length(ind); 
        end
        function postFeatureLearning(obj) %#ok<MANU>
           return;
        end
    end
    methods (Static = true)
        function W = getTrOutliers(v,kappa)
            AVG = 0;
            nrO = length(v);W = ones(1,nrO);
            for i=1:1:10
                STD = sqrt(sum(W.*((v).^2))/sum(W));
                IP = normpdf(v,AVG,STD);
                L = (1/(sqrt(2*pi)*STD))*exp(-0.5*kappa^2);
                W = IP./(IP+L);
            end
        end
        function out = sum2one(v)
            out = v/nansum(v);
        end
    end
end


% function [NX,COST,W,ind,Tind,Find] = prepPredictorTraining(obj,X,TargetX)
%             obj.TargX = TargetX;
%             NX = nan*zeros(size(X));
%             W = nan*zeros(size(X));
%             ind = find(~isnan(X));
%             
%             range = max(X)-min(X);
%             range = range/obj.Bins;
%             obj.RangeX = range;
%             %medX = median(X(ind));
%             %madX = mad(X(ind));
%             %range = mad(X(ind));
%            % select closest examples
%             dif = abs(TargetX-X(ind));
%             Tind2 = find(dif<=range);
%             Tind = ind(Tind2);
%             %[~,ind2] = sort(dif,'ascend');
%             %Tind2 = ind2(1:obj.KN);
%             %Tind = ind(Tind2);
%             NX(Tind) = obj.PosClass;
%             tmp = dif(Tind2);
%             if max(tmp)==0
%                W(Tind) = 1;
%             else
%                tmp = tmp-min(tmp);tmp = tmp/max(tmp);tmp = 1-tmp;
%                %mval = min(tmp(tmp>0));tmp(tmp==0) = mval;
%                W(Tind) = tmp;
%             end
%            % opposite examples
%             Find2 = find(dif>range);
%             %Find2 = ind2(end-obj.KN+1:end);
%             Find = ind(Find2);
%             NX(Find) = -1*obj.PosClass;
%             tmp = dif(Find2);
%             tmp = tmp-min(tmp);tmp = tmp/max(tmp);
%             %mval = min(tmp(tmp>0));tmp(tmp==0) = mval;
%             W(Find) = tmp;
%             ind = union(Tind,Find);
%             W(ind) = W(ind)+1;
%            % Defining COST 
%             COST.ClassNames = [obj.PosClass -1*obj.PosClass];
%             COST.ClassificationCosts = zeros(2,2);
%             %COST.ClassificationCosts(1,2) = 1;
%             %COST.ClassificationCosts(2,1) = 1;
%             COST.ClassificationCosts(1,2) = 1-(length(Tind)/length(ind));
%             COST.ClassificationCosts(2,1) = 1-(length(Find)/length(ind));
%             W = nan*zeros(size(X));
%             W(Tind) = 1-(length(Tind)/length(ind));
%             W(Find) = 1-(length(Find)/length(ind));
%         end

% function [MS,MH] = match2Seq(obj,TargX,FS)
%             if isempty(obj.SeqX), error('Sequence not trained'); end
%             n = length(obj.SeqX);
%             diff = zeros(1,length(obj.SeqX));
%             in = zeros(1,length(obj.SeqX));
%             nT = size(FS,1);
%             MS = zeros(nT,n);
%             MH = zeros(nT,n);
%             for i=1:n
%                 diff(i) = abs(obj.SeqX{i}.TargX-TargX);
%                 if diff(i)<=obj.SeqX{i}.RangeX
%                    in(i) = 1;
%                 else
%                    in(i) = -1;%continue;
%                 end
%                 [PX,PXC] = predict(obj.SeqX{i},FS);
%                 [MS(:,i),MH(:,i)] = Pred2Match(obj.SeqX{i},in(i),PX,PXC);
%             end
%             %good = find(in==1);
%             good = 1:n;
%             W = 1./diff(good);
%             sumW = sum(W);
%             MS = MS(:,good);
%             MS = sum(repmat(W,nT,1).*MS,2)/sumW;
%             MH = MH(:,good);
%             MH = sum(repmat(W,nT,1).*MH,2)/sumW;  
%         end
