classdef SexDemiModel < superClass
    properties
        Tag = 'Sex';% Tag on which the demiModel operates, needed to link to input data in an excell file
        Level = 'Genomic';% Genomic or Genotype, what is the level of DNA that the demi model is active (used to discriminate demimodels used in the creation of the base face or not)
        D;% Data matrix, rows = samples, collumns = features
        G;% Group Vector, -1 males, 1 females
        C;% Conditioning variables if present
        Kfold = 10;
        CrossVal = 'DOB_SCV';% options are SCV, DOB_SCV
        NNDistance = 'Euclidean';% options are Euclidean or Mahalanobis
        Balance= 'ADASYN';% options are none ADASYN, SMOTE, SMOTE EER, RUSBoost, EUSBoost
        BalanceRuns = 10;% Number of times the data is rebalanced in a random process, resuling classifiers will be bagged
        Classifier = 'PLSR RIP ROC';% options are PLSR (BRIM) RIP ROC, PLSR (BRIM) RIP NB,
        TrainingOut;% structure to store output information from training phase and classifier models
        DStd;% standard deviations (eigen values in PCA space e.g.)
        Ref = 'zeros';% options are zeros, average
    end
    properties (Dependent = true)
       nrS;%nr of Samples
       nrF;%nr of Features
       nrC;% nr of conditioning variables
       DW;% standarizing weights for features (to switch to Mahalanobis distance e.g.)
       DRef;% Reference feature values (to be used in computing RIP values e.g.)
       BalanceDegree; % The degree of (im)balance
    end
    properties (Hidden = true, Dependent = true)% variables to be hidden (algorithmic importance only)
        nrBal;% the actual number of balance runs to perform
    end
    methods % Constructor
        function obj = SexDemiModel(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
         function out = get.nrS(obj)
            if isempty(obj.D), out = 0; return; end
            out = size(obj.D,1);
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
         function out = get.nrBal(obj)
             if strcmp(obj.Balance,'none');out = 1;return; end
             out = obj.BalanceRuns;
         end
    end
    methods % General Interface Functions
        function eliminateNAN(obj)
            index = (1:obj.nrSamples);
            [i,~] = find(isnan(obj.G));
            i = unique(i);
            index = setdiff(index,i);
            obj.D = obj.D(index,:);
            obj.G = obj.G(index,:);
            obj.C = obj.C(index,:);
        end
        function out = getFolds(obj)
           switch obj.CrossValid
               case 'SCV'
                   out = crossvalind('Kfold',obj.nrSamples,obj.Kfold);
               case 'DOB_SCV'
                   switch obj.NNDistance
                       case 'Euclidean'
                          d = obj.D; 
                       case 'Mahalanobis'
                          d = obj.D./repmat(obj.DW,obj.nrSamples,1); 
                   end
                   out = DOB_SCV(obj.Kfold,d,obj.G);
           end
        end
        function trainModel(obj)
           % Remove samples with sex unknown
           eliminateNAN(obj);
           % splitting data into crossvalidation folds
           F = getFolds(obj);
           FoldData = cell(K,obj.nrBal);
           parfor k=1:obj.Kfold
               TestInd = find(F==k);
               TrInd = setdiff(1:obj.nrSamples,TestInd); %#ok<*PFBNS>
               switch obj.Balance
                   case 'none'
                       out = Imbalanced(obj,TestInd,TrInd);
                   case 'ADASYN'
               end
               
           end
           %bagging(obj);
        end
        function bagging(obj)
            % bag the different classifiers over different folds
        end
    end
    methods % implemented Balancing methods
        function out = Imbalanced(obj,TestInd,TrInd)
            % do nothing in terms of balancing, just build the classifier
            out = buildClassifier(obj,TestInd,TrInd);
        end
    end
    methods % implementing classifiers
        function out = buildClassifier(obj,TestInd,TrInd)
            switch obj.Method
                   case {'PLSR RIP ROC' 'PLSR RIP NB'}
                       out = PLSR_RIP(obj,TestInd,TrInd);
            end
        end
        function out = getGmean(obj,model,TestInd)
            est = model.Test;
            ex = obj.G(TestInd);
            pos = find(ex==1);
            
            
            
        end
        function out = PLSR_RIP(obj,TestInd,TrInd)
            out = [];
            B = obj.D(TrInd,:);
            if obj.nrC == 0
               A = obj.G(TrInd);
            else
               A = [obj.C(TrInd,:) obj.G(TrInd)];
               [A,B] = eliminateNAN(A,B);
            end
            % establish regression
            [~,~,~,~,out.M] = plsregress(A,B,size(A,2));
            % estimate RIP values for testcases
            out.Rip = getRIP(obj.DRef,obj.DW,obj.D(TestInd,:)',out.M);
            switch obj.Method
                case 'PLSR RIP ROC'
                    [x,y,t,auc] = perfcurve(obj.G(TestInd),out.Rip,1);
                    if auc<0.5
                       [x,y,t,auc] = perfcurve(obj.G(TestInd),out.Rip,-1);
                    end
                    yn = 1-y;
                    d = abs(x-yn);
                    [~,indmin] = min(d);
                    out.eer = (x(indmin)+yn(indmin))/2;
                    out.T = t(indmin);
                    out.auc = auc;
                    out.Test = ones(1,length(TestInd));
                    out.Test(rip<=out.T) = -1;% perform classification to build test measures from
                case 'PLSR RIP NB'
            end
            out.Gmean = getGmean(obj,out,TestInd);
        end
    end
end