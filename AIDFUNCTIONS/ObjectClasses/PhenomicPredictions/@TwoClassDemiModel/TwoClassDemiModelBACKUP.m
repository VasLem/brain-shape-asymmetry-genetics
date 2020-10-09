classdef TwoClassDemiModel < superClass
    properties
        Tag = 'Sex';% Tag on which the demiModel operates, needed to link to input data in an excell file
        Level = 'Genomic';% Genomic or Genotype, what is the level of DNA that the demi model is active (used to discriminate demimodels used in the creation of the base face or not)
        D;% Data matrix, rows = samples, collumns = features
        G;% Group Vector, -1 males, 1 females
        C;% Conditioning variables if present
        Kfold = 10;
        CrossVal = 'DOB_SCV';% options are SCV, DOB_SCV
        NNDistance = 'Euclidean';% options are Euclidean or Mahalanobis
        Balance= 'none';% options are none ADASYN, SMOTE, SMOTE EER, RUSBoost, EUSBoost
        BalanceRuns = 10;% Number of times the data is rebalanced in a random process, resuling classifiers will be bagged
        BalanceTh = 75;% preset threshold for the maximum tolerated degree of class imbalance ratio
        BalanceB = 0.9;% desired level of balance after balancing, = 1 means fully balanced
        Classifier = 'PLSR RIP ROC';% options are PLSR (BRIM) RIP ROC, PLSR (BRIM) RIP NB,
        TrainingOut;% structure to store output information from training phase and classifier models
        DStd;% standard deviations (eigen values in PCA space e.g.)
        Ref = 'zeros';% options are zeros, average
        Status = 0;% Status can be 0 (initialization), 1 (crossvalidation) or 2 (full deploying)
        BaggedClassifiers;% A cell structure containing all bagged classifiers
        DataCleaning = 'none';% options are none , tomek or ENN
        DataTransformation = 'none';% Data transformation prior to everything options are none (basically using eucledian distances), Mahalanobis, LMNN (large margin nearest neighbor)
        TrainingDataTransformation = 'none';% Data transformation within the Training loops
    end
    properties (Dependent = true)
       nrS;%nr of Samples
       nrF;%nr of Features
       nrC;% nr of conditioning variables
       DW;% standarizing weights for features (to switch to Mahalanobis distance e.g.)
       MinorLab;% label for the minority group
       MajorLab;% label for the majority group
       nrMinor;% number of Minor samples
       nrMajor;% numberof Major samples
       BalanceDegree; % The degree of (im)balance
       DRef;% Reference feature values (to be used in computing RIP values e.g.)
    end
    properties (Hidden = true, Dependent = true)% variables to be hidden (algorithmic importance only)
        nrBal;% the actual number of balance runs to perform
        nrS1;% nr of samples in -1 class
        nrS2;% nr of samples in 1 class
        nrOuterRuns;% the number 
    end
    properties (Hidden = true)
        trainD;% data to train with when being rebalanced
        trainG;% G info to train with when being rebalanced
        trainC;% Covariate info to train with when being rebalanced
    end
    methods % Constructor
        function obj = TwoClassDemiModel(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
         function out = get.nrS(obj)
            if isempty(obj.D), out = 0; return; end
            out = size(obj.D,1);
         end
         function out = get.nrS1(obj)
             if isempty(obj.D)||isempty(obj.G), out = 0; return; end
             out = length(find(obj.G==-1));
         end
         function out = get.nrS2(obj)
             if isempty(obj.D)||isempty(obj.G), out = 0; return; end
             out = length(find(obj.G==1));
         end
         function out = get.MinorLab(obj)
             if obj.nrS1==obj.nrS2, out = -1; return; end
             if obj.nrS1<obj.nrS2, out = -1; return; end
             out = 1;
         end
         function out = get.MajorLab(obj)
             if obj.nrS1==obj.nrS2, out = 1; return; end
             if obj.nrS1<obj.nrS2, out = 1; return; end
             out = -1;
         end
         function out = get.nrMinor(obj)
              if isempty(obj.D)||isempty(obj.G), out = 0; return; end
              out = length(find(obj.G==obj.MinorLab));
         end
         function out = get.nrMajor(obj)
              if isempty(obj.D)||isempty(obj.G), out = 0; return; end
              out = length(find(obj.G==obj.MajorLab));
         end
         function out = get.nrF(obj)
            if isempty(obj.D), out = 0; return; end
            out = size(obj.D,2);
         end
         function out = get.nrC(obj)
            if isempty(obj.C), out = 0; return; end
            out = size(obj.C,2);
         end
         function out = get.DW(obj)
            if isempty(obj.D), out = []; return; end
            if isempty(obj.DStd), out = ones(1,obj.nrF); return; end
            out = obj.DStd;
         end
         function out = get.DRef(obj)
            if isempty(obj.D), out = []; return; end
            switch obj.Ref
                case 'zeros'
                    out = zeros(1,obj.nrF);
                case 'average'
                    out = mean(obj.D);
            end
         end
         function out = get.BalanceDegree(obj)
             out = obj.nrMinor/obj.nrMajor;
         end
    end
    methods % General Interface Functions
        function eliminateNAN(obj)
            index = (1:obj.nrS);
            [i,~] = find(isnan(obj.G));
            i = unique(i);
            index = setdiff(index,i);
            obj.D = obj.D(index,:);
            obj.G = obj.G(index,:);
            if obj.nrC == 0, return; end
            obj.C = obj.C(index,:);
        end
        function out = getFolds(obj,D,G)
           nrD = size(D,1);
           switch obj.CrossVal
               case 'SCV'
                   out = crossvalind('Kfold',nrD,obj.Kfold);
               case 'DOB_SCV'
                   
                   out = DOB_SCV(obj.Kfold,d,G);
           end
        end
        function transformData(obj)
           switch obj.DataTransformation
               case 'none'
                   return;
               case 'Mahalanobis'
                   obj.D = obj.D./repmat(obj.DW,size(obj.D,1),1);
           end
        end
        function crossvalidClassifier(obj)
           obj.Status = 0;
           % Remove samples with sex unknown
           eliminateNAN(obj);
           % splitting data into crossvalidation folds
           obj.Status = 1;
           F = getFolds(obj,obj.D,obj.G);
           EvalData = cell(obj.Kfold,1);
           parfor_progress(obj.Kfold);
           parfor k=1:obj.Kfold
               TestInd = find(F==k);
               TrInd = setdiff(1:obj.nrS,TestInd); %#ok<*PFBNS>
               out = trainClassifier(obj,TrInd,TestInd);
               EvalData{k} = out.Eval;
               parfor_progress;
           end
           parfor_progress(0);
           obj.TrainingOut = TwoClassDemiModel.extractEvaluations(EvalData);
           obj.Status = 2;
           %bagging(obj);
        end
    end
    methods % implementing classifiers
        function out = trainClassifier(obj,TrInd,TestInd)
            out.Classifiers = cell(obj.Kfold,obj.BalanceRuns);
            for b=1:1:obj.BalanceRuns
                switch obj.Balance
                    case 'none' % pass on the same (original) data every time
                        imbalanced(obj,TrInd);
                    case 'ADASYN'
                        ADASYN(obj,TrInd);
                end
                switch obj.Classifier
                       case {'PLSR RIP ROC' 'PLSR RIP NB'}
                           runout = trainPLSR_RIP(obj);
                end
                out.Classifiers(:,b) = runout(:);
            end
            if obj.Status == 1% in cross validation mode
               [out.Classes,out.Condidence] = TwoClassDemiModel.classifyWithModel(out.Classifiers,obj.D(TestInd,:),obj.Classifier,obj.DRef,obj.DW);
               out.Eval = getEvaluation(obj,out,TestInd);
            end
        end
        function out = trainPLSR_RIP(obj)
                % devide into folds to extract training parameters
                F = getFolds(obj,obj.trainD,obj.trainG);
                out = cell(obj.Kfold,1);
                % for every fold train a classifier
                for k2=1:1:obj.Kfold
                    TestInd = find(F==k2);
                    TrInd = setdiff(1:size(obj.trainD,1),TestInd);
                    B = obj.trainD(TrInd,:);
                    if obj.nrC == 0
                       A = obj.trainG(TrInd);
                    else
                       A = [obj.trainC(TrInd,:) obj.trainG(TrInd)];
                       [A,B] = eliminateNAN(A,B);% G does not have any nan values left, but C might have
                    end
                    % establish regression
                    [~,~,~,~,class.M] = plsregress(A,B,size(A,2));
                    class.TestRip = getRIP(obj.DRef',obj.DW',obj.trainD(TestInd,:)',class.M);
                    switch obj.Classifier
                        case 'PLSR RIP ROC'
                            % perform ROC analysis
                            [x,y,t,auc] = perfcurve(obj.trainG(TestInd),class.TestRip,1);
                            lower = -1;higher = 1;
                            if auc<0.5
                               [x,y,t,auc] = perfcurve(obj.trainG(TestInd),class.TestRip,-1);
                               lower = 1;higher = -1;
                            end
                            yn = 1-y;d = abs(x-yn);[~,indmin] = min(d);
                            class.eer = (x(indmin)+yn(indmin))/2;
                            class.T = t(indmin);
                            class.auc = auc;
                            class.lower = lower;
                            class.higher = higher;
                        case 'PLSR RIP NB'
                    end
                    out{k2,1} = class;
                end
        end
        function out = getEvaluation(obj,res,TestInd)
            % by definition, minority is the positive group, and the
            % majority the negative group
            posest = res.Classes(obj.G(TestInd)==obj.MinorLab);
            negest = res.Classes(obj.G(TestInd)==obj.MajorLab);
            TP = length(find(posest==obj.MinorLab));
            FN = length(posest)-TP;
            TN = length(find(negest==obj.MajorLab));
            FP = length(negest)-TN;
            out.OA = (TP+TN)/length(TestInd);
            out.Precision = TP/(TP+FP);
            out.Recall = TP/(TP+FN);
            out.G = sqrt((TP/(TP+FN))*(TN/(TN+FP)));
        end
        function out = classify(obj,D)
           out = TwoClassDemiModel.classifyWithModel(obj.BaggedClassifiers,D,obj.Classifier,obj.DRef,obj.DW);
        end
    end
    methods % implemented Balancing methods
        function imbalanced(obj,TrInd)
           obj.trainD = obj.D(TrInd,:);
           obj.trainG = obj.G(TrInd,:);
           if ~obj.nrC==0, obj.trainC = obj.C(TrInd,:);end 
        end
        function ADASYN(obj,TrInd)
            D = obj.D(TrInd,:);
            G = obj.G(TrInd,:);
            K = 10;% number of closest examples to be found within ADASYN
            if ~obj.nrC==0, C = obj.C(TrInd,:);end
            indm = find(G==obj.MinorLab);
            indM = find(G==obj.MajorLab);
            nrm = length(indm);
            nrM = length(indM); %#ok<*PROP>
            dcrit = nrm/nrM;
            if dcrit >= obj.BalanceTh
               obj.trainD = D; obj.trainG = G;
               if ~obj.nrC==0, obj.trainC = C; end
               return;
            end
            nrSyn = (nrM-nrm)*obj.BalanceB;
            switch obj.NNDistance
                   case 'Euclidean'
                        d = D; 
                   case 'Mahalanobis'
                        d = D./repmat(obj.DW,size(D,1),1); 
            end
            r = zeros(1,nrm);
            for i=1:1:nrm
                mex = indm(i);
                distances = sqrt(sum((repmat(d(mex,:),nrm+nrM,1)-d).^2,2));
                [~,index] = sort(distances,'ascend');
                index = index(2:K+1);
                r(i) = length(intersect(index,indM))/K;
            end
            r = r/sum(r);
            g = round(r*nrSyn);
            totg = sum(g);
            Dadd = zeros(totg,size(D,2));
            Gadd = obj.MinorLab*ones(totg,1);
            dm = d(indm,:);
            Dm = D(indm,:);
            if ~obj.nrC == 0, 
               Cm = C(indm,:);
               Cadd = zeros(totg,size(C,2));
            end
            counter = 0;
            for i=1:1:nrm
                if g(i)==0, continue; end
                distances = sqrt(sum((repmat(dm(i,:),nrm,1)-dm).^2,2));
                [~,index] = sort(distances,'ascend');
                index = index(2:K+1);
                for l=1:g(i)
                    counter = counter+1;
                    nb = randi(K,1,1);
                    alpha = rand(1);
                    Dadd(counter,:) = Dm(i,:) + (Dm(index(nb),:)-Dm(i,:))*alpha;
                    if ~obj.nrC==0, Cadd(counter,:) = Cm(i,:) + (Cm(index(nb),:)-Cm(i,:))*alpha;end
                end
            end
            D = [D; Dadd];
            G = [G; Gadd];
            if ~obj.nrC==0, C = [C; Cadd]; end
            ind = randperm(size(D,1));
            obj.trainD = D(ind,:);
            obj.trainG = G(ind,:);
            if ~obj.nrC==0, obj.trainC = C(ind,:); end
        end
    end
    methods (Static = true)
        function [Classes,Confidence] = classifyWithModel(model,D,method,DRef,DW)
            model = model(:);
            nrM = numel(model);
            nrT = size(D,1);
            out = nan*zeros(nrT,nrM);
            %parfor i=1:nrM
            for i=1:nrM
                class = model{i};
                switch method
                    case 'PLSR RIP ROC'% simple thresholding on rip values
                        rip = getRIP(DRef',DW',D',class.M);
                        est = class.lower*ones(1,nrT);
                        est(rip>class.T) = class.higher;
                end
                out(:,i) = est(:);
            end
            % majority voting
            C1votes = out==-1;
            C1votes = sum(C1votes,2)/nrM;
            C2votes = out==1;
            C2votes = sum(C2votes,2)/nrM;
            Classes = ones(nrT,1);
            index = C1votes>=C2votes;
            Classes(index) = -1;
            Confidence = C2votes;
            Confidence(index) = C1votes(index);
        end
        function out = extractEvaluations(eval)
           eval = eval(:);
           nrE = length(eval);
           out.OA = nan*zeros(1,nrE);
           out.Precision = nan*zeros(1,nrE);
           out.Recall = nan*zeros(1,nrE);
           out.G = nan*zeros(1,nrE);
           for i=1:1:nrE
              out.OA(i) = eval{i}.OA;
              out.Precision(i) = eval{i}.Precision;
              out.Recall(i) = eval{i}.Recall;
              out.G(i) = eval{i}.G;
           end
           out.AvgOA = mean(out.OA);
           out.StdOA = std(out.OA);
           out.AvgPrecision = mean(out.Precision);
           out.StdPrecision = std(out.Precision);
           out.AvgRecall = mean(out.Recall);
           out.StdRecall = std(out.Recall);
           out.AvgG = mean(out.G);
           out.StdG = std(out.G);
        end
    end
end