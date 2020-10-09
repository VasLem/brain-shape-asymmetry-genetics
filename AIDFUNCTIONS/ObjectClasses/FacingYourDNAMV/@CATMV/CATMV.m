classdef CATMV < PPMMV
    properties
        ClassifierType = 'SVMRBF';
        X;
    end
    properties (Dependent = true)
        Predictor;
    end
    properties (Hidden = true)
        Classifier;
        CVloss = 0;
        OptSVMRBF = true;
        FST=single(0.45);
        REGFST = single(0.47);
        OffsetPXC;
        ScalePXC;
    end
    properties(Hidden = true, Dependent = true)
        treeTempl;
        predictortype;
    end
    properties (Abstract = true)
        PosClass; 
    end
    methods % CONSTRUCTOR
        function obj = CATMV(varargin)
            obj = obj@PPMMV(varargin{:});         
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
    end
    methods % PREDICTOR
        function trainPredictor(obj,X,FS)
           % PREPARATION
             BKX = X;
             [X,COST,W,ind,~,~] = prepPredictorTraining(obj,X);
           % FEATURE SELECTION
             FS = getFS(obj,FS);
             featureSelection(obj,X(ind),BKX(ind),FS(ind,:));
           % ACCUMULATE FEATURES
             sFS = selectFeatures(obj,FS);
           % TRAINING
             fitClassifier(obj,sFS(ind,:),X(ind),COST,W(ind));
           % LEARNING HOW TO NORMALIZE  
             [~,PXC] = predict(obj,FS);
             learnPXCNorm(obj,PXC);
        end
        function trainSVMRBF(obj,X,FS,K)
             rng(1);
             if nargin<4, K = 5;end
           % PREPARATION
             obj.ClassifierType = 'SVMRBF';
             BKX = X;
             [X,COST,W,ind,~,~] = prepPredictorTraining(obj,X);
           % FEATURE SELECTION
             FS = getFS(obj,FS);
             startnFS = obj.maxFS;
             obj.maxFS = obj.maxFS*2;%size(FS,2);
             featureSelection(obj,X(ind),BKX(ind),FS(ind,:));
           % GET INITIAL PARAMETER SETTINGS
             fst = obj.FSInfo(:,startnFS);
             index = find(obj.FSInfo>=fst);
             init.fst = log(fst);
             init.c = log(1);
             %COST.ClassificationCosts(1,2) = 1;
             %COST.ClassificationCosts(2,1) = 1;
             tmp = fitcsvm(FS(ind,obj.FSInd(index)),X(ind),'Standardize',false,'KernelFunction','rbf','BoxConstraint',1,...
                          'KernelScale','auto','Cost',COST,'Weights',W(ind),'Prior','empirical');              
             init.scale = log(tmp.KernelParameters.Scale);
            % OPTIMIZE PARAMETERS
            switch obj.OptSVMRBF
                case false
                    obj.FSInd = obj.FSInd(index);obj.maxFS = length(obj.FSInd);
                    obj.Classifier = fitcsvm(FS(ind,obj.FSInd),X(ind),'Standardize',false,'KernelFunction','rbf','BoxConstraint',1,...
                                             'KernelScale',tmp.KernelParameters.Scale,'Cost',COST,'Weights',W(ind),'Prior','empirical');
                    cv = crossval(obj.Classifier,'Kfold',K);
                    obj.CVloss = kfoldLoss(cv);
                case true
                     COST.ClassificationCosts(1,2) = 1;
                     COST.ClassificationCosts(2,1) = 1;
                     %[obj.CVloss,sigma,C,fst] = CATMV.optParSVMRBFRAND(FS(ind,obj.FSInd),X(ind),obj.FSInfo,COST,W(ind),init,K);
                     [obj.CVloss,sigma,C,fst] = CATMV.optParSVMRBFINIT(FS(ind,obj.FSInd),X(ind),obj.FSInfo,COST,W(ind),init,K);
                     index = find(obj.FSInfo>=fst);
                     obj.FSInd = obj.FSInd(index);
                     obj.maxFS = length(obj.FSInd);
                     obj.Classifier = fitcsvm(FS(ind,obj.FSInd),X(ind),...
                                          'KernelFunction','rbf','BoxConstraint',C,...
                                          'Cost',COST,'Weights',W(ind),'Prior','empirical',...
                                          'KernelScale',sigma);
            end
            % LEARNING HOW TO NORMALIZE  
             [~,PXC] = predict(obj,FS);
             learnPXCNorm(obj,PXC);                   
        end
        function fitClassifier(obj,FS,X,COST,W)
            rng(1);
            switch obj.ClassifierType
                case 'SVM'
                    obj.Classifier = fitcsvm(FS,X,'Standardize',true,'Cost',COST,'Weights',W);
                    %obj.Classifier = fitPosterior(obj.Classifier);
                case 'SVMRBF'
                    switch obj.OptSVMRBF
                        case false  
                            %COST.ClassificationCosts(1,2) = 1;
                            %COST.ClassificationCosts(2,1) = 1;
                            tmp = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',1,...
                                                     'KernelScale','auto','Cost',COST,'Weights',W,'Prior','empirical');
                            obj.Classifier = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf','BoxConstraint',1,...
                                                     'KernelScale',tmp.KernelParameters.Scale,'Cost',COST,'Weights',W,'Prior','empirical');
                           cv = crossval(obj.Classifier,'Kfold',5);
                           obj.CVloss = kfoldLoss(cv);
                        case true
                            COST.ClassificationCosts(1,2) = 1;
                            COST.ClassificationCosts(2,1) = 1;
                            %COST.ClassificationCosts
                            %[obj.CVloss,sigma,C,nFS] = CATMV.optSVMRBF(FS,X,COST,W,5);
                            %obj.FSInd = obj.FSInd(1:nFS);
                            %obj.Classifier = fitcsvm(FS(:,1:nFS),X,'KernelFunction','rbf',...
                            %                         'KernelScale',sigma,'BoxConstraint',C,...
                            %                         'Cost',COST,'Weights',W);
                            [obj.CVloss,sigma,C] = CATMV.optSVMRBF1(FS,X,COST,W,5);
                            obj.Classifier = fitcsvm(FS,X,'KernelFunction','rbf','Standardize',false,...
                                                     'KernelScale',sigma,'BoxConstraint',C,...
                                                     'Cost',COST,'Weights',W);
                            
                    end
                    %cv = crossval(obj.Classifier,'Kfold',5);
                    %obj.CVloss = kfoldLoss(cv);
                case 'ROBUSTBOOST'
                    obj.Classifier = fitensemble(FS,X,'RobustBoost',obj.nrTrees,obj.treeTempl,'Cost',COST,'Weights',W);
                case 'RUSBOOST'
                    obj.Classifier = fitensemble(FS,X,'RUSBoost',obj.nrTrees,obj.treeTempl,'Cost',COST,'Weights',W);
                 otherwise
                    error('WRONG TYPE OF CLASSIFIER');
            end
        end
        function out = testPredictor(obj,X,FS,display)
            % PREPARATION
              [~,~,Tind,Find,Tval] = prepPredictorTesting(obj,X);
              FS = getFS(obj,FS);
            % PREDICTION
              [PX,PXCN] = predict(obj,FS);
              PXCN = normPXC(obj,PXCN);
              switch Tval
                  case -1
                      PXCN = PXCN(:,2);
                  case 1
                      PXCN = PXCN(:,1);
              end
            % EVALUATION  
               res = PX==Tval;
               %[out.R1,out.R10,out.R20,out.CR] =  PPMMV.getRANKSTAT(TMatches,FMatches);
               [out.EER,out.G,out.AUC,out.XX,out.YY,out.TH,out.Y] = PPMMV.getEER(PXCN(Tind),PXCN(Find));% CHECK ORDER
               [out.HPREC,out.HREC,out.HG,out.HTNF,out.HX,out.HY,out.HACC] = PPMMV.getHardClass(res(Tind),res(Find));
            % VISUALIZATION
               if display
                  figure;hold on;grid on;plot(0:0.1:1,0:0.1:1,'k-');
                  plot(0:0.1:1,1:-0.1:0,'k--');
                  plot(out.XX,out.YY,'b-','LineWidth',1.5);
                  plot(out.HX,out.HY,'r.','MarkerSize',20);
                  title(['EER: ' num2str(out.EER) ' HG: ' num2str(out.HG)]);
                  %figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
                  %xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
                  %title('Identification');
                  %plot(1:100,1:100,'k--');
                  %plot(1:100,out.CR,'g-','LineWidth',1.5);
               end
        end
        function [PX,PXC] = predict(obj,FS)
            % ACCUMULATE FEATURES
              FS = selectFeatures(obj,FS);
            % PREDICT  
              [PX,PXC] = predict(obj.Predictor,FS);
              %PXCN = PPMMV.normVAL(PXC);% NOT SURE THIS IS THE RIGHT THING TO DO IN THE BIOMETRICS CASE
              if iscell(PX), PX = cell2mat(PX); end
        end
    end
    methods % BIOMETRICS
        function [SM,HM] = Pred2Match(obj,X,PX,PXC)
            if isscalar(X)
               if isnan(X), SM = nan*zeros(size(PX));HM = SM; return; end
               HM = (1-(repmat(X,size(PX))==PX));
               PXC = normPXC(obj,PXC);
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
        function [SOFTTMatches,SOFTFMatches,HARDTMatches,HARDFMatches,wF,wH,wL] = biometricMatchesWeights(obj,TestDepVar,TestX)
           % INITIALIZE 
            nrT = size(TestX,1);
            TestX = prepBiometricTest(obj,TestX);
           % WEIGHTS 
            %wH = obj.nHits;
            wH = obj.wpAR;
            wL = (0.5-obj.CVloss);
            wL = max(wL,0);
            wL = wL/0.5;
           % FEATURES 
            RIPF = getRIP(obj,TestDepVar);
            PCF = getDepVar(obj,TestDepVar,0);
            %FS = [PCF, RIPF];
            if ~strcmp(obj.Type,'SNPMV')
               FS = [RIPF, PCF];
            else
               %FS = PCF;
               FS = getFS(obj,TestDepVar);
            end
           % PREDICTIONS
            [PX,PXC] = predict(obj,FS);
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
                [SM,HM] = Pred2Match(obj,TestX(t),PX,PXC);
                SOFTTMatches(:,t) = SM(t);
                SOFTFMatches(:,t,:) = SM(tInd);
                HARDTMatches(:,t) = HM(t);
                HARDFMatches(:,t,:) = HM(tInd);
            end
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
    end
    methods (Static = true)
        function [fval,sigma,C] = optSVMRBF1(FS,X,COST,W,K)
            if nargin<5, K = 4; end
            if isscalar(K), K = cvpartition(length(X),'KFold',K); end
            tmp = fitcsvm(FS,X,'Standardize',false,'KernelFunction','rbf',...
                          'KernelScale','auto','Cost',COST,'Weights',W);
            initscale = log(tmp.KernelParameters.Scale);
            initc = log(1);
            minfn = @(z)kfoldLoss(fitcsvm(FS,X,'CVPartition',K,'Standardize',false,...
                       'KernelFunction','rbf','BoxConstraint',exp(z(2)),...
                       'Cost',COST,'Weights',W,...
                       'KernelScale',exp(z(1))));
            options = psoptimset('UseParallel', true, 'CompletePoll', 'on', 'Vectorized', 'off','TolMesh',5e-5,'Display','iter');
            [searchmin, fval] = patternsearch(minfn,[initscale, initc],[],[],[],[],-5*ones(2,1),5*ones(2,1),options);
            %[searchmin, fval] = patternsearch(minfn,randn(1,2),[],[],[],[],-5*ones(2,1),5*ones(2,1),options);
            z = exp(searchmin);
            sigma = z(1);
            C = z(2);
        end
        function [fval,sigma,C,fst] = optParSVMRBFINIT(FS,X,fsinfo,COST,W,init,K)
            %X = X(ind);FS = FS(ind,:);W = W(ind);K = 4;fsinfo = obj.FSInfo;
            if nargin<6, K = 5; end
            if isscalar(K), K = cvpartition(length(X),'KFold',K); end       
            minfn = @(z)CATMV.kfoldLossSVMRBF(FS,X,K,COST,W,fsinfo,z);
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
            minfn = @(z)CATMV.kfoldLossSVMRBF(FS,X,K,COST,W,fsinfo,z);
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
    end
end

    
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