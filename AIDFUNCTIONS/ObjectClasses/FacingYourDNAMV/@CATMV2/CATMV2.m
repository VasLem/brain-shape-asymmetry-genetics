classdef CATMV2 < PPMMV2
    properties
        ClassifierType = 'SVMRBF';
        ClassifierCascade = false;
        X;
        Beta = [];
    end
    properties (Dependent = true)
        Predictor;
    end
    properties (Hidden = true)
        Classifier;
        CascadeInd;
        FSHits;
        CVloss = 1;
        CasCVloss = 1;
        pCV = 1;
        CVEER = 1;
        CasCVEER = 1;
        CVG = 0;
        CVAUC = 0;
        OptSVMRBF = true;
        FST=single(0.45);
        REGFST = single(0.47);
        OffsetPXC;
        ScalePXC;
    end
    properties(Hidden = true, Dependent = true)
        treeTempl;
        predictortype;
        nCascade;
    end
    properties (Abstract = true)
        PosClass; 
    end
    methods % CONSTRUCTOR
        function obj = CATMV2(varargin)
            obj = obj@PPMMV2(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.Predictor(obj)
            out = obj.Classifier;
        end
        function out = get.treeTempl(obj)
            if obj.nFS==0, out = []; return; end
            out = templateTree('MaxNumSplits',obj.nrSplits,'Type','Classification'); 
        end
        function out = get.predictortype(obj) %#ok<MANU>
           out = 'classifier';
        end
        function obj = set.X(obj,in)
           % Conversion for memory saving 
           %ind = isnan(in); in = int8(in);in(ind) = -128;
           if strcmp(obj.Type,'SNPMV'), obj.X = single(in);return; end
           if strcmp(obj.Type,'SEXMV'), obj.X = single(in);return; end
           obj.X = in;
        end
        function out = get.X(obj)
            %ind = obj.X==-128;
            %out = single(obj.X);
            %out(ind) = nan;
            out = obj.X;
        end
        function out = get.nCascade(obj)
            out = length(obj.CascadeInd);
        end
        function out = get.FSHits(obj)
           if isempty(obj.FSHits), out = obj.Hits; return;end
           out = obj.FSHits;
        end
    end
    methods % PREDICTOR
        function trainPredictor(obj,X,FS,K,COV,t)
            if nargin<6, t = 0; end
            if nargin<5, COV = []; end
            if ~isempty(COV)
               index = PPMMV2.notNAN(COV);
               X = X(index);
               FS = PPMMV2.selectDepVar(FS,index);
               COV = COV(index,:);
               [FS,obj.Beta] = redDepVar(obj,FS,COV);
            else
               obj.Beta = [];
            end
            switch obj.ClassifierCascade
                case false
                    trainSVMRBFSingle(obj,X,FS,K,t);
                case true
                    trainSVMRBFCascade(obj,X,FS,K,t);
            end
        end
        function trainSVMRBFSingle(obj,X,FS,K,t)     
             if nargin<4, K = 5;end
             if nargin<5, t = 0; end
           % PREPARATION
             obj.ClassifierType = 'SVMRBF';
             BKX = X;
             [X,COST,W,ind,~,~] = prepPredictorTraining(obj,X);
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
        function trainSVMRBFCascade(obj,X,FS,K,t)
             if nargin<4, K = 5;end
             if nargin<5, t = 0; end
           % PREPARATION
             obj.ClassifierType = 'SVMRBF';
             [X,COST,W,ind,~,~] = prepPredictorTraining(obj,X);
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
        function out = testPredictor(obj,X,FS,COV,display)
            if nargin < 5, display = false;end
            if nargin < 4, COV = [];end
            if ~isempty(COV)
               index = PPMMV2.notNAN(COV);
               X = X(index);
               FS = PPMMV2.selectDepVar(FS,index);
               COV = COV(index,:);
               FS = redDepVarBeta(obj,FS,COV);
            end
            switch obj.ClassifierCascade
                case false
                    out = testSinglePredictor(obj,X,FS,display);
                case true
                    out = testCascadePredictor(obj,X,FS,display);
            end
        end
        function out = testSinglePredictor(obj,X,FS,display)
            % PREPARATION
              [~,~,Tind,Find,Tval] = prepPredictorTesting(obj,X);
            % PREDICTION
              [PX,~,PXCN] = predict(obj,FS,length(X));   
              switch Tval
                  case -1
                      PXCN = PXCN(:,2);
                  case 1
                      PXCN = PXCN(:,1);
              end
            % EVALUATION  
               res = PX==Tval;
               %[out.R1,out.R10,out.R20,out.CR] =  PPMMV2.getRANKSTAT(TMatches,FMatches);
               [out.EER,out.G,out.AUC,out.XX,out.YY,out.TH,out.Y,out.pAUC] = PPMMV2.getEER(PXCN(Tind),PXCN(Find));% CHECK ORDER
               [out.HPREC,out.HREC,out.HG,out.HTNF,out.HX,out.HY,out.HACC] = PPMMV2.getHardClass(res(Tind),res(Find));
            % VISUALIZATION
               if display
                  figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                  plot(0:0.1:1,1:-0.1:0,'k--');
                  plot(out.XX,out.YY,'b-','LineWidth',1.5);
                  plot(out.HX,out.HY,'r.','MarkerSize',20);
                  title(['EER: ' num2str(out.EER) ' HG: ' num2str(out.HG) ' CV: ' num2str(obj.CVloss) ' AUC: ' num2str(out.AUC)]);
                  %figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
                  %xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
                  %title('Identification');
                  %plot(1:100,1:100,'k--');
                  %plot(1:100,out.CR,'g-','LineWidth',1.5);
               end
        end
        function [PX,PXC,PXCN] = predict(obj,FS,nS)
            switch obj.ClassifierCascade
                case false
                    % ACCUMULATE FEATURES
                      FS = getFS(obj,FS);% Hit Based selection
                      FS = selectFeatures(obj,FS);% Feature Selection
                    % PREDICT  
                      [PX,PXC] = predict(obj.Predictor,FS);
                      PXCN = CATMV2.statNormPXC(PXC,obj.OffsetPXC,obj.ScalePXC);
                case true
                    %error('To be implemented'); 
                    allPX = zeros(nS,obj.nCascade);
                    allPXC = zeros(nS,2,obj.nCascade);
                    allPXCN = zeros(nS,2,obj.nCascade);
                    for i=1:obj.nCascade
                      forFS = getDepVar(obj,FS,obj.CascadeInd(i));
                      [allPX(:,i),allPXC(:,:,i)] = predict(obj.Classifier{i},forFS);
                      allPXCN(:,:,i) = CATMV2.statNormPXC(allPXC(:,:,i),obj.OffsetPXC(:,i)',obj.ScalePXC(:,i)');   
                    end
                    %wC = ones(1,obj.nCascade);% UNWEIGHTED COMBINATION!
                    wC = -log(obj.CVloss./(1-obj.CVloss));
                    %wC = -log(obj.CVEER./(1-obj.CVEER));
                    %[~,index] = sort(wC,'descend');
                    %n = round(0.5*length(wC));
                    %wC(index(end-n+1:end)) = 0;
                    wC(obj.CVloss>0.4) = 0;
                    %wC(index(1:n)) = 1;
                    % MAJORITY VOTING
                    PX = sum(repmat(wC,nS,1).*allPX,2)./sum(wC);
                    PX(PX<0) = -1;
                    PX(PX>=0) = 1;
                    PXC = nan*zeros(nS,2);
                    PXCN = nan*zeros(nS,2);
                    for i=1:2
                        PXC(:,i) = sum(repmat(wC,nS,1).*squeeze(allPXC(:,i,:)),2)./sum(wC);
                        PXCN(:,i) = sum(repmat(wC,nS,1).*squeeze(allPXCN(:,i,:)),2)./sum(wC);
                    end
            end
        end
        function [out] = testCascadePredictor(obj,X,FS,display,val)
            % PREPARATION
              [~,~,Tind,Find,Tval] = prepPredictorTesting(obj,X);
              nT = length(X);
              if nargin<5, 
                  val = obj.CVloss(obj.CascadeInd); 
                  %val = obj.pCV; 
                  lim = [0 0.5];
              else
                  val = val(obj.CascadeInd);
                  lim = [min(val) max(val)];
              end
              if display
                 figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'-','Color',[0.5 0.5 0.5]);
                 plot(0:0.1:1,1:-0.1:0,'-','Color',[0.5 0.5 0.5]);
                 colormap(gca,'jet');
                 legstr = cell(1,obj.nCascade+4);
                 legstr{1} = ' ';
                 legstr{2} = ' ';
                 map = colormap('jet');
                 [~,sortind] = sort(val,'ascend');
                 %lim = [min(val) max(val)];
                 set(gca,'clim',lim);
                 colind = floor(1:64/length(val):64);
                 tmpcolors = map(colind,:);
                 colors = tmpcolors;
                 for i=1:1:length(sortind)
                    colors(sortind(i),:) = tmpcolors(i,:);    
                 end
                 %colorbar;
              end
              
            % FOR EACH PREDICTOR IN THE CASCADE
              allPX = zeros(nT,obj.nCascade);
              allPXCN = zeros(nT,obj.nCascade);
              allRes = zeros(nT,obj.nCascade);
              out = cell(1,obj.nCascade+4);        
              for i=1:obj.nCascade
                  forFS = getDepVar(obj,FS,obj.CascadeInd(i));
                  [PX,PXCN] = predict(obj.Classifier{i},forFS);
                  PXCN = CATMV2.statNormPXC(PXCN,obj.OffsetPXC(:,i)',obj.ScalePXC(:,i)');
                  switch Tval
                      case -1
                          PXCN = PXCN(:,2);
                      case 1
                          PXCN = PXCN(:,1);
                  end
                  allPX(:,i) = PX;
                  allPXCN(:,i) = PXCN;
                  res = PX==Tval;
                  allRes(:,i) = 1-res;
                  [out{i}.EER,out{i}.G,out{i}.AUC,out{i}.XX,out{i}.YY,out{i}.TH,out{i}.Y] = PPMMV2.getEER(PXCN(Tind),PXCN(Find));% CHECK ORDER
                  [out{i}.HPREC,out{i}.HREC,out{i}.HG,out{i}.HTNF,out{i}.HX,out{i}.HY,out{i}.HACC] = PPMMV2.getHardClass(res(Tind),res(Find));
                  if display
                     %scatter(out{i}.XX,out{i}.YY,5,repmat(val(i),length(out{i}.XX),1));
                     plot(out{i}.XX,out{i}.YY,'-','LineWidth',1,'Color',colors(i,:));
                     %plot(out{i}.HX,out{i}.HY,'r.','MarkerSize',10);
                     legstr{i+2} = ['L' num2str(obj.Levels(obj.CascadeInd(i))) 'C' num2str(obj.Clusters(obj.CascadeInd(i))) ' ' num2str(val(i))];
                  end
              end
            % FULL EVALUATION UNWEIGHTED
              % DEFINE WEIGHTS
              wC = ones(1,obj.nCascade);
              % MAJORITY VOTING
              PX = sum(repmat(wC,nT,1).*allPX,2)./sum(wC);
              PX(PX<0) = -1;
              PX(PX>=0) = 1; 
              % ACCUMULATE HARD CLASSIFICATION
              i = i+1;
              PXN = sum(repmat(wC,nT,1).*allRes,2)./sum(wC);
              [out{i}.EER,out{i}.G,out{i}.AUC,out{i}.XX,out{i}.YY,out{i}.TH,out{i}.Y,out{i}.pAUC] = PPMMV2.getEER(PXN(Tind),PXN(Find));% CHECK ORDER 
              res = PX==Tval;
              [out{i}.HPREC,out{i}.HREC,out{i}.HG,out{i}.HTNF,out{i}.HX,out{i}.HY,out{i}.HACC] = PPMMV2.getHardClass(res(Tind),res(Find));
              if display
                 plot(out{i}.XX,out{i}.YY,'m-.','LineWidth',2);
                 %plot(cas{1}.HX,cas{1}.HY,'k.','MarkerSize',20);
                 legstr{i+2} = 'UW AHC';
              end
              % ACCUMULATE SOFT CLASSIFICATION
              i= i+1;
              PXN = sum(repmat(wC,nT,1).*allPXCN,2)./sum(wC);
              [out{i}.EER,out{i}.G,out{i}.AUC,out{i}.XX,out{i}.YY,out{i}.TH,out{i}.Y,out{i}.pAUC] = PPMMV2.getEER(PXN(Tind),PXN(Find));% CHECK ORDER 
              res = PX==Tval;
              [out{i}.HPREC,out{i}.HREC,out{i}.HG,out{i}.HTNF,out{i}.HX,out{i}.HY,out{i}.HACC] = PPMMV2.getHardClass(res(Tind),res(Find));
              if display
                 plot(out{i}.XX,out{i}.YY,'m-','LineWidth',2);
                 %plot(cas{1}.HX,cas{1}.HY,'k.','MarkerSize',20);
                 legstr{i+2} = 'UW ASC';
              end
              
            % FULL EVALUATION WEIGHTED
              % DEFINE WEIGHTS
              wC = -log(obj.CVloss(obj.CascadeInd)./(1-obj.CVloss(obj.CascadeInd)));
              wC(wC<0) = 0;
              %wC = -log10(obj.pCCA(obj.CascadeInd));
              %wC(obj.CVEER>0.45) = 0;
              %wC = -log(obj.pCV./(1-obj.pCV));
              %wC = zeros(1,obj.nCascade);
              %wC(obj.pCV<=0.3) = 1;
              %wC(obj.CVloss>0.35) = 0;
%               [~,index] = sort(wC,'descend');
%               n = round(0.5*length(wC));
%               wC(index(end-n+1:end)) = 0;
%               wC(index(1:n)) = 1;
              % MAJORITY VOTING
              PX = sum(repmat(wC,nT,1).*allPX,2)./sum(wC);
              PX(PX<0) = -1;
              PX(PX>=0) = 1; 
              % ACCUMULATE HARD CLASSIFICATION
              i = i+1;
              PXN = sum(repmat(wC,nT,1).*allRes,2)./sum(wC);
              [out{i}.EER,out{i}.G,out{i}.AUC,out{i}.XX,out{i}.YY,out{i}.TH,out{i}.Y,out{i}.pAUC] = PPMMV2.getEER(PXN(Tind),PXN(Find));% CHECK ORDER 
              res = PX==Tval;
              [out{i}.HPREC,out{i}.HREC,out{i}.HG,out{i}.HTNF,out{i}.HX,out{i}.HY,out{i}.HACC] = PPMMV2.getHardClass(res(Tind),res(Find));
              if display
                 plot(out{i}.XX,out{i}.YY,'k-.','LineWidth',2);
                 %plot(cas{1}.HX,cas{1}.HY,'k.','MarkerSize',20);
                 legstr{i+2} = 'W AHC';
              end
              % ACCUMULATE SOFT CLASSIFICATION
              i= i+1;
              PXN = sum(repmat(wC,nT,1).*allPXCN,2)./sum(wC);
              [out{i}.EER,out{i}.G,out{i}.AUC,out{i}.XX,out{i}.YY,out{i}.TH,out{i}.Y,out{i}.pAUC] = PPMMV2.getEER(PXN(Tind),PXN(Find));% CHECK ORDER 
              res = PX==Tval;
              [out{i}.HPREC,out{i}.HREC,out{i}.HG,out{i}.HTNF,out{i}.HX,out{i}.HY,out{i}.HACC] = PPMMV2.getHardClass(res(Tind),res(Find));
              if display
                 plot(out{i}.XX,out{i}.YY,'k-','LineWidth',2);
                 %plot(cas{1}.HX,cas{1}.HY,'k.','MarkerSize',20);
                 legstr{i+2} = 'W ASC';
                 legend(legstr,'Location','SouthEast');
                 title(['EER: ' num2str(out{i}.EER) ' HG: ' num2str(out{i}.HG) ' AUC: ' num2str(out{i}.AUC)]);
              end
        end
    end
    methods % BIOMETRICS
        function [SM,HM] = Pred2Match(obj,X,PX,PXC)
            if isscalar(X)
               if isnan(X), SM = nan*zeros(size(PX));HM = SM; return; end
               HM = (1-(repmat(X,size(PX))==PX));
               %PXC = normPXC(obj,PXC);
               switch X
                   case -1
                       SM = PXC(:,2);
                   case 1
                       SM = PXC(:,1);
               end
               return;
            end
            SM = nan*zeros(size(PXC,1),size(X,1));
            HM = nan*zeros(size(PXC,1),size(X,1));
            parfor i=1:1:size(X,1)
                [SM(:,i),HM(:,i)] = Pred2Match(obj,X(i),PX,PXC) 
            end
        end
        function [SOFTTMatches,SOFTFMatches,HARDTMatches,HARDFMatches,wF,wH,wL] = biometricMatchesWeights(obj,FS,TestX)
           % INITIALIZE 
            nrT = size(TestX,1);
            TestX = prepBiometricTest(obj,TestX);
           % WEIGHTS 
            %wH = obj.nHits;
            wH = obj.nHits;
            wL = -log(obj.CVloss./(1-obj.CVloss));
            if wL<0, wL = 0; end
            %wL = wL*obj.PosClassBal;
           % PREDICTIONS
            [PX,~,PXCN] = predict(obj,FS);
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
%                 if isnan(wF(t))
%                    TestX(t) = PX(t); 
%                    wF(t) = 1-getFreq(obj,TestX(t));
%                 end
                if isnan(wF(t)), continue; end
                [SM,HM] = Pred2Match(obj,TestX(t),PX,PXCN);
                SOFTTMatches(:,t) = SM(t);
                SOFTFMatches(:,t,:) = SM(tInd);
                HARDTMatches(:,t) = HM(t);
                HARDFMatches(:,t,:) = HM(tInd);
            end
        end
    end
    methods % EVOMORPH
        function [SM,HM,wF,wH,wL] = matchSolutions(obj,TestX,FS,nrT)
            TestX = prepEvoMorphTest(obj,TestX);
            if isnan(TestX), SM = nan*zeros(nrT,1);HM = nan*zeros(nrT,1);wF = nan;wH = nan;wL = nan;return;end
            [PX,~,PXCN] = predict(obj,FS);
            [SM,HM] = Pred2Match(obj,TestX,PX,PXCN);
            wF = 1-getFreq(obj,TestX);
            wH = obj.nHits;
            wL = max(-log(obj.CVloss./(1-obj.CVloss)),0);
            %wL = wL*obj.PosClassBal;
        end
    end
    methods % INTERFACING
        function learnPXCNorm(obj,PXC)
           obj.OffsetPXC  = min(PXC);
           obj.ScalePXC = max(PXC)-min(PXC);
        end
        function out = normPXC(obj,PXC)
           n = size(PXC,1);
           out =  (PXC-repmat(obj.OffsetPXC,n,1))./repmat(obj.ScalePXC,n,1);
           out(out<0) = 0;
           out(out>1) = 1;
        end
        function filterClassifierCascade(obj,th)
            
            
            
        end
    end
    methods (Static = true)
        function out = trainSVMRBF(FS,X,COST,W,K,opt,Tval,t)
            if nargin < 8, t = 0; end
            %FS = FS(ind,obj.FSInd);
            %X = X(ind);
            %W = W(ind);
            %Tval = obj.PosClass;
            switch opt
                case true
                    %COST.ClassificationCosts(1,2) = 1;
                    %COST.ClassificationCosts(2,1) = 1;
                otherwise 
            end
            COST.ClassificationCosts(1,2) = 1;
            COST.ClassificationCosts(2,1) = 1;
            C = 1;
            tmp = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
                          'KernelScale','auto','Cost',COST,'Weights',W,'Prior','empirical');
            sigma = tmp.KernelParameters.Scale;          
            if ~isobject(K), K = cvpartition(length(X),'KFold',K); end          
            switch opt
                case false
                    out.Classifier = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
                                             'KernelScale',sigma,'Cost',COST,'Weights',W,'Prior','empirical');                                       
                case true
                    init.scale = log(sigma);
                    init.c = log(1);
                    COST.ClassificationCosts(1,2) = 1;
                    COST.ClassificationCosts(2,1) = 1;
                    [out.cv,sigma,C] = CATMV2.optSVMRBF1(FS,X,COST,W,K,init);
                    out.Classifier = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
                                             'KernelScale',sigma,'Cost',COST,'Weights',W,'Prior','empirical');                     
            end
            % CROSS VALIDATION LOSS
             out.cv = kfoldLoss(fitcsvm(FS,X,'CVPartition',K,...
                               'KernelFunction','rbf','BoxConstraint',C,...
                               'Cost',COST,'Weights',W,'Prior','empirical',...
                               'KernelScale',sigma));
            % CROSS VALIDATION EER
             [out.EER,out.G,out.AUC] = CATMV2.kfoldEER(FS,X,K,COST,W,C,sigma,Tval,false);
             [~,PXC] = predict(out.Classifier,FS);
             out.OffsetPXC  = min(PXC);
             out.ScalePXC = max(PXC)-min(PXC);
             if t==0, out.pCV = 0;out.pEER = 0;out.pG = 0;out.pAUC = 0;return; end
             cvcount = zeros(1,t);
             %eercount = nan*zeros(1,t);
             %gcount = nan*zeros(1,t);
             %auccount = nan*zeros(1,t);
             parfor i=1:t
                rng(i);
                ind = randperm(length(X));
                forcv = kfoldLoss(fitcsvm(FS,X(ind),'CVPartition',K,...
                               'KernelFunction','rbf','BoxConstraint',C,...
                               'Cost',COST,'Weights',W(ind),'Prior','empirical',...
                               'KernelScale',sigma));
                %disp(num2str(forcv));
                cvcount(i) = forcv<=out.cv;
                %if mod(i,50)==0
                %   pCV = (sum(cvcount(1:i))+1)/(i+1); 
                %   acc = 10/i;
                %   if pCV>=acc, break;end
                %end  
             end
             %out.pCV = (sum(cvcount(1:i))+1)/(i+1); 
             out.pCV = (sum(cvcount(1:t))+1)/(t+1); 
        end
        function [fval,sigma,C] = optSVMRBF1(FS,X,COST,W,K,init)
            if nargin<5, K = 4; end
            if ~isobject(K), K = cvpartition(length(X),'KFold',K); end
            minfn = @(z)kfoldLoss(fitcsvm(FS,X,'CVPartition',K,'Standardize',false,...
                       'KernelFunction','rbf','BoxConstraint',exp(z(2)),...
                       'Cost',COST,'Weights',W,...
                       'KernelScale',exp(z(1)),'Prior','empirical'));
            options = psoptimset('UseParallel', true, 'CompletePoll', 'on', 'Vectorized', 'off','TolMesh',5e-5,'Display','off');
            [searchmin, fval] = patternsearch(minfn,[init.scale, init.c],[],[],[],[],-5*ones(2,1),5*ones(2,1),options);
            z = exp(searchmin);
            sigma = z(1);
            C = z(2);
        end
        function [fval,sigma,C,fst] = optParSVMRBFINIT(FS,X,fsinfo,COST,W,init,K)
            %X = X(ind);FS = FS(ind,:);W = W(ind);K = 4;fsinfo = obj.FSInfo;
            if nargin<6, K = 5; end
            if isscalar(K), K = cvpartition(length(X),'KFold',K); end       
            minfn = @(z)CATMV2.kfoldLossSVMRBF(FS,X,K,COST,W,fsinfo,z);
            options = psoptimset('UseParallel', true, 'CompletePoll', 'on', 'Vectorized', 'off','TolMesh',5e-5,'Display','iter');
            lB = [-5*ones(2,1);log(min(fsinfo,[],2))];
            uB = [5*ones(2,1);log(max(fsinfo,[],2))];
            [searchmin, fval] = patternsearch(minfn,[init.scale, init.c, init.fst],[],[],[],[],lB,uB,options);
            %[searchmin, fval] = patternsearch(minfn,randn(1,3),[],[],[],[],lB,uB,options);
            %[searchmin, fval] = patternsearch(minfn,[init.scale, init.c, init.fst],[],[],[],[],[],[],options);
            z = exp(searchmin);
            sigma = z(1);
            C = z(2);
            fst = z(3);
        end
        function [fval,sigma,C,fst] = optParSVMRBFRAND(FS,X,fsinfo,COST,W,init,K)
            %X = X(ind);FS = FS(ind,:);W = W(ind);K = 4;fsinfo = obj.FSInfo;
            if nargin<6, K = 5; end
            if isscalar(K), K = cvpartition(length(X),'KFold',K); end       
            minfn = @(z)CATMV2.kfoldLossSVMRBF(FS,X,K,COST,W,fsinfo,z);
            options = psoptimset('UseParallel', true, 'CompletePoll', 'on', 'Vectorized', 'off','TolMesh',5e-5,'Display','off');
            lB = [-5*ones(2,1);log(min(fsinfo,[],2))];
            uB = [5*ones(2,1);log(max(fsinfo,[],2))];
            t = 16;
            searchmin = nan*zeros(t,3);
            fval = nan*zeros(t);
            parfor i=1:t
                  init = zeros(1,3);
                  for j=1:1:3
                      init(j) = lB(j) + (uB(j)-lB(j)).*rand(1,1);
                  end
                  [searchmin(i,:), fval(i)] = patternsearch(minfn,init,[],[],[],[],lB,uB,options);
            end
            %[searchmin, fval] = patternsearch(minfn,randn(1,3),[],[],[],[],lB,uB,options);
            %[searchmin, fval] = patternsearch(minfn,[init.scale, init.c, init.fst],[],[],[],[],[],[],options);
            [~,ind] = min(fval);
            fval = fval(ind(1));
            searchmin = searchmin(ind(1),:);
            z = exp(searchmin);
            sigma = z(1);
            C = z(2);
            fst = z(3);
        end
        function loss = kfoldLossSVMRBF(FS,X,K,COST,W,fsinfo,z)
                 loss = kfoldLoss(fitcsvm(FS(:,fsinfo>=exp(z(3))),X,'CVPartition',K,...
                                  'KernelFunction','rbf','BoxConstraint',exp(z(2)),...
                                  'Cost',COST,'Weights',W,'Prior','empirical',...
                                  'KernelScale',exp(z(1))));
        end
        function [EER,G,AUC] = kfoldEER(fs,x,K,COST,w,C,sigma,Tval,display)
                 if nargin < 9, display = false; end
                 %fres = cell(1,K.NumTestSets)
                 eer = zeros(1,K.NumTestSets);
                 g = zeros(1,K.NumTestSets);
                 auc = zeros(1,K.NumTestSets);
                 %X = zeros(1,K.NumTestSets);
                 %Y = zeros(1,K.NumTestSets);
                 %fs = FS(ind,obj.FSInd);
                 %x = X(ind);
                 %w = W(ind);        
                 %C = 1;
                 %sigma = 10;
                 if display
                   figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                   plot(0:0.1:1,1:-0.1:0,'k--');
                 end
                 for i=1:K.NumTestSets
                     TrInd = K.training(i);
                     TestInd = K.test(i);
                     % BUILD CLASSIFIER
                     Classifier = fitcsvm(fs(TrInd,:),x(TrInd),'Standardize',false,'KernelFunction','rbf','BoxConstraint',C,...
                                          'KernelScale',sigma,'Cost',COST,'Weights',w(TrInd),'Prior','empirical');
                     % CLASSIFY                 
                     [PX,PXCN] = predict(Classifier,fs(TestInd,:));
                     Find = find(x(TestInd)==-1*Tval);
                     Tind = find(x(TestInd)==Tval);  
                     switch Tval
                          case -1
                              PXCN = PXCN(:,2);
                          case 1
                              PXCN = PXCN(:,1);
                     end
                     % EVALUATION  
                     %res = PX==Tval;
                     %[out.R1,out.R10,out.R20,out.CR] =  PPMMV2.getRANKSTAT(TMatches,FMatches);
                     %[fres{i}.EER,fres{i}.G,fres{i}.AUC] = PPMMV2.getEER(PXCN(Tind),PXCN(Find));% CHECK ORDER
                     [eer(i),g(i),auc(i),X,Y] = PPMMV2.getEER(PXCN(Tind),PXCN(Find));% CHECK ORDER
                     if display
                         plot(X,Y,'b-','LineWidth',1.5);
                     end
                     %[fres{i}.HPREC,fres{i}.HREC,fres{i}.HG,fres{i}.HTNF,fres{i}.HX,fres{i}.HY,fres{i}.HACC] = PPMMV2.getHardClass(res(Tind),res(Find));
                 end
%                  EER = mean(eer);
%                  G = mean(g);
%                  AUC = mean(auc);
                 EER = max(eer);
                 G = min(g);
                 AUC = min(auc);
        end
        function out = statNormPXC(PXC,OffsetPXC,ScalePXC)
            n = size(PXC,1);
            out =  (PXC-repmat(OffsetPXC,n,1))./repmat(ScalePXC,n,1);
            out(out<0) = 0;
            out(out>1) = 1; 
        end      
    end
end



% CVO = cvpartition(species,'KFold',10);
%       err = zeros(CVO.NumTestSets,1);
%       for i = 1:CVO.NumTestSets
%           trIdx = CVO.training(i);
%           teIdx = CVO.test(i);
%           ytest = classify(meas(teIdx,:),meas(trIdx,:),species(trIdx,:));
%           err(i) = sum(~strcmp(ytest,species(teIdx)));
%       end
%       cvErr = sum(err)/sum(CVO.TestSize);
    
%     function out = crossvalidate(obj,X,FS,K)
%            if nargin<4, K = 10;end
%            n = length(X);
%            F = crossvalind('Kfold',n,K);
%            numfs = 10:10:obj.maxFS;
%            nfs = length(numfs);
%            EER = nan*zeros(K,nfs);
%            HG = nan*zeros(K,nfs);
%            parfor k=1:K
%                Tind = find(F==k);
%                Trind = setdiff(1:n,Tind);
%                forobj = clone(obj);
%                forX = X(Trind); %#ok<PFBNS>
%                forFS = FS(Trind,:); %#ok<PFBNS>
%                forTX = X(Tind);
%                forTFS = FS(Tind,:);
%                % PREPARATION
%                forBKX = forX;
%                [forX,COST,W,ind,~,~] = prepPredictorTraining(forobj,forX);
%                % FEATURE SELECTION
%                featureSelection(forobj,forX(ind),forBKX(ind),forFS(ind,:));
%                fsind = forobj.FSInd;
%                tmpEER = nan*zeros(1,nfs);
%                tmpHG = nan*zeros(1,nfs);
%                for i=1:1:nfs
%                    % PREP FS INDEX
%                      forobj.FSInd = fsind(1:numfs(i)); %#ok<PFBNS>
%                    % ACCUMULATE FEATURES
%                      sFS = selectFeatures(forobj,forFS); 
%                    % TRAINING
%                      fitClassifier(forobj,sFS(ind,:),forX(ind),COST,W(ind));
%                    % LEARNING HOW TO NORMALIZE  
%                      [~,PXC] = predict(forobj,forFS);
%                      learnPXCNorm(forobj,PXC);
%                    % TESTING  
%                      test = testPredictor(forobj,forTX,forTFS,false);
%                      tmpEER(i) = test.EER;
%                      tmpHG(i)  = test.HG;
%                end
%                EER(k,:) = tmpEER;
%                HG(k,:) = tmpHG;
%            end
%            out.EER = EER;
%            out.HG = HG;
%            [minEER,ind] = min(mean(EER,1));
%            out.nFSEER = numfs(ind);
%            out.minEER = minEER;
%            out.maxEER = max(mean(EER,1));
%            out.minEERHG = mean(HG(:,ind));
%            tmp = EER(:,ind);
%            out.FailedEER = sum(tmp>=0.5);      
%            [maxHG,ind] = max(mean(HG,1));
%            out.nFSHG = numfs(ind);
%            out.maxHG = maxHG;
%            out.minHG = min(mean(HG,1));
%            out.maxHGEER = mean(EER(:,ind));
%            tmp = HG(:,ind);
%            out.FailedHG = sum(tmp<=0.5);
%         end


% function fitClassifier(obj,FS,X,COST,W)
%             rng(1);
%             switch obj.ClassifierType
%                 case 'SVM'
%                     obj.Classifier = fitcsvm(FS,X,'Standardize',true,'Cost',COST,'Weights',W);
%                     %obj.Classifier = fitPosterior(obj.Classifier);
%                 case 'SVMRBF'
%                     switch obj.OptSVMRBF
%                         case false  
%                             %COST.ClassificationCosts(1,2) = 1;
%                             %COST.ClassificationCosts(2,1) = 1;
%                             tmp = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',1,...
%                                                      'KernelScale','auto','Cost',COST,'Weights',W,'Prior','empirical');
%                             obj.Classifier = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',1,...
%                                                      'KernelScale',tmp.KernelParameters.Scale,'Cost',COST,'Weights',W,'Prior','empirical');
%                            cv = crossval(obj.Classifier,'Kfold',5);
%                            obj.CVloss = kfoldLoss(cv);
%                         case true
%                             COST.ClassificationCosts(1,2) = 1;
%                             COST.ClassificationCosts(2,1) = 1;
%                             %COST.ClassificationCosts
%                             %[obj.CVloss,sigma,C,nFS] = CATMV2.optSVMRBF(FS,X,COST,W,5);
%                             %obj.FSInd = obj.FSInd(1:nFS);
%                             %obj.Classifier = fitcsvm(FS(:,1:nFS),X,'KernelFunction','rbf',...
%                             %                         'KernelScale',sigma,'BoxConstraint',C,...
%                             %                         'Cost',COST,'Weights',W);
%                             [obj.CVloss,sigma,C] = CATMV2.optSVMRBF1(FS,X,COST,W,5);
%                             obj.Classifier = fitcsvm(FS,X,'KernelFunction','rbf','Standardize',false,...
%                                                      'KernelScale',sigma,'BoxConstraint',C,...
%                                                      'Cost',COST,'Weights',W);
%                             
%                     end
%                     %cv = crossval(obj.Classifier,'Kfold',5);
%                     %obj.CVloss = kfoldLoss(cv);
%                 case 'ROBUSTBOOST'
%                     obj.Classifier = fitensemble(FS,X,'RobustBoost',obj.nrTrees,obj.treeTempl,'Cost',COST,'Weights',W);
%                 case 'RUSBOOST'
%                     obj.Classifier = fitensemble(FS,X,'RUSBoost',obj.nrTrees,obj.treeTempl,'Cost',COST,'Weights',W);
%                  otherwise
%                     error('WRONG TYPE OF CLASSIFIER');
%             end
%         end