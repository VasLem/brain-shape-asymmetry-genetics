classdef TwoClassDemiModel < superClass
    properties
        Tag = 'Sex';% Tag on which the demiModel operates, needed to link to input data in an excell file
        Level = 'Genomic';% Genomic or Genotype, what is the level of DNA that the demi model is active (used to discriminate demimodels used in the creation of the base face or not)
        D;% Data matrix, rows = samples, collumns = features
        G;% Group Vector, -1 males, 1 females
        C;% Conditioning variables if present
        Kfold = 10;% cross validation folds
        Klmnn = 5;% number of neigbors in LMNN
        Kadasyn = 5;% number of neighbors ADASYN
        CrossVal = 'DOB_SCV';% options are SCV, DOB_SCV
        NNDistance = 'Euclidean';% options are Euclidean or Mahalanobis
        Balance= 'none';% options are none ADASYN, SMOTE, SMOTE EER, RUSBoost, EUSBoost
        BalanceRuns = 10;% Number of times the data is rebalanced in a random process, resuling classifiers will be bagged
        BalanceTh = 75;% preset threshold for the maximum tolerated degree of class imbalance ratio
        BalanceB = 0.9;% desired level of balance after balancing, = 1 means fully balanced
        Classifier = 'PLSR RIP ROC';% options are PLSR (BRIM) RIP ROC, PLSR (BRIM) RIP NB, KNN, Decision Tree, LSVM (linear support vector machine), KSVM (kernel support vector machine)
        TrainingOut;% structure to store output information from training phase and classifier models
        DStd;% standard deviations (eigen values in PCA space e.g.)
        Ref = 'zeros';% options are zeros, average
        Status = 0;% Status can be 0 (initialization), 1 (crossvalidation) or 2 (full deploying)
        BaggedClassifiers;% A cell structure containing all bagged classifiers
        DataCleaning = 'none';% options are none , tomek or ENN
        DataTransformation = 'none';% Data transformation prior to everything options are none (basically using eucledian distances), seuclidean (Mahalanobis), LMNN (large margin nearest neighbor)
        TrainingDataTransformation = 'none';% Data transformation within the Training loops, options are none, seuclidean or LMNN
        Verbose = true;
        NPManova = [];
    end
    properties (Dependent = true)
       nrS;%nr of Samples
       nrF;%nr of Features
       nrC;% nr of conditioning variables
       DW;% Data transformation
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
        L;% data transformation matrix
        transD;% data after transformation, hidden variable in order not to loose original data
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
            if isempty(obj.DStd), out = std(obj.D); return; end
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
        function getNPManova(obj,D,G)
            % this is a typical non parametric ANOVA test to figure out the
            % difficulty of the classification problem based on distances
            % between observations in feature space
            obj.NPManova = oneWayNPManova;
            indm = find(G==obj.MinorLab);
            indM = find(G==obj.MajorLab);
            D = [D(indm,:);D(indM,:)];
            obj.NPManova.D = squareform(pdist(D));
            obj.NPManova.n = [length(indm) length(indM)];
            perform(obj.NPManova);
        end
        function out = getFolds(obj,D,G)
           nrD = size(D,1);
           switch obj.CrossVal
               case 'SCV'
                   out = crossvalind('Kfold',nrD,obj.Kfold);
               case 'DOB_SCV'
                   out = DOB_SCV(obj.Kfold,D,G);
           end
        end
        function transformData(obj)
           switch obj.DataTransformation
               case 'none'
                   obj.transD = obj.D;
               case 'seuclidean'
                   obj.transD = obj.D./repmat(obj.DW,size(obj.D,1),1);
               case 'LMNN'
                   disp('Performing LMNN data transformation, this can take a while...');
                   D = obj.D./repmat(obj.DW,size(obj.D,1),1);% pre normalizing
                   tic;obj.L = lmnn2(D',obj.G',obj.Klmnn,eye(size(obj.D,2)),'quiet',1);
                   obj.transD = (obj.L*D')';
                   toc;
           end
        end
        function Info = crossvalidClassifier(obj)
           obj.Status = 0;
           % Remove samples with sex unknown
           eliminateNAN(obj);
           % transform Data
           transformData(obj);
           % perform NP-MANOVA
           getNPManova(obj,obj.transD,obj.G);
           % splitting data into crossvalidation folds
           obj.Status = 1;
           F = getFolds(obj,obj.transD,obj.G);
           EvalData = cell(obj.Kfold,1);
           Info = cell(obj.Kfold,1);
           parfor_progress(obj.Kfold);
           parfor k=1:obj.Kfold
               TestInd = find(F==k);
               TrInd = setdiff(1:obj.nrS,TestInd); %#ok<*PFBNS>
               out = trainClassifier(obj,TrInd,TestInd);
               EvalData{k} = out.Eval;
               Info{k} = out.Info;
               parfor_progress;
           end
           parfor_progress(0);
           obj.TrainingOut = TwoClassDemiModel.extractEvaluations(EvalData);
           obj.TrainingOut.Info = TwoClassDemiModel.extractInfo(Info);
           obj.Status = 2;
           %bagging(obj);
        end
    end
    methods % implementing classifiers
        function out = trainClassifier(obj,TrInd,TestInd)
            out.Classifiers = cell(obj.Kfold,obj.BalanceRuns);
            out.Info.Removed = zeros(1,obj.BalanceRuns);
            for b=1:1:obj.BalanceRuns
                switch obj.Balance
                    case 'none' % pass on the same (original) data every time
                        imbalanced(obj,TrInd);
                    case 'ADASYN'
                        ADASYN(obj,TrInd);
                end
                switch obj.DataCleaning
                    case 'none'
                        out.Info.Removed(b) = 0;
                    case 'Tomec'
                        out.Info.Removed(b) = Tomec(obj);
                    case 'ENN'
                        out.Info.Removed(b) = ENN(obj);
                end
                switch obj.Classifier
                    case {'PLSR RIP ROC' 'PLSR RIP NB'}
                        runout = trainPLSR_RIP(obj);
                    case 'BRIM RIP ROC'
                        runout = trainBRIM_RIP(obj);
                    case 'KNN'
                        runout = trainKNN(obj);
                    case 'Decision Tree'
                        runout = trainDecisionTree(obj);
                    case 'LSVM'
                        runout = trainLSVM(obj);
                    case 'KSVM'
                        runout = trainKSVM(obj);
                end
                out.Classifiers(:,b) = runout(:);
            end
            if obj.Status == 1% in cross validation mode
               [out.Classes,out.Condidence] = TwoClassDemiModel.classifyWithModel(out.Classifiers,obj.transD(TestInd,:),obj.Classifier,obj.DRef);
               out.Eval = getEvaluation(obj,out,obj.G(TestInd));
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
                    class.TestRip = TwoClassDemiModel.getRIP(obj.DRef',obj.trainD(TestInd,:)',class.M);
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
        function out = trainBRIM_RIP(obj)
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
                    % BRIM LOOP
                    W = ones(1,size(A,1));
                    RIP = A(:,end);
                    %Kappa = 2.5;
                    for i=1:1:3
                        index = find(W==1);
                        A(:,end) = RIP(:);
                        [~,~,~,~,class.M] = plsregress(A(index,:),B(index,:),size(A,2));
                        RIP = TwoClassDemiModel.getRIP(obj.DRef',B',class.M);
                        W = ones(1,size(A,1));
                        for j=1:1:3
                            index = find(W==1);
                            % normalizing between -1 and 1
                            avgRIP = mean(RIP(index));RIP = RIP-avgRIP;
                            mRIP = min(RIP(index));MRIP = max(RIP(index));RIP = RIP/(MRIP-mRIP);
                            stdRIP = std(RIP(index));
                            % taking care of outliers
                            W = ones(1,size(A,1));
                            W(abs(RIP)>2.5*stdRIP) = 0;
                        end
                    end
                    class.TestRip = TwoClassDemiModel.getRIP(obj.DRef',obj.trainD(TestInd,:)',class.M);
                    switch obj.Classifier
                        case 'BRIM RIP ROC'
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
                        case 'BRIM RIP NB'
                    end
                    out{k2,1} = class;
                end
        end
        function out = trainKNN(obj)
            % devide into folds to extract training parameters
                F = getFolds(obj,obj.trainD,obj.trainG);
                out = cell(obj.Kfold,1);
                % for every fold train a classifier
                for k2=1:1:obj.Kfold
                    TestInd = find(F==k2);
                    TrInd = setdiff(1:size(obj.trainD,1),TestInd);
                    classG = zeros(1,10);
                    for KK=1:10% finding optimal number of nearest neighbors
                        mdl = ClassificationKNN.fit(obj.trainD(TrInd,:),obj.trainG(TrInd),'NumNeighbors',KK);
                        pred.Classes = predict(mdl,obj.trainD(TestInd,:));
                        testout = getEvaluation(obj,pred,obj.trainG(TestInd));
                        classG(KK) = testout.G;
                    end
                    [class.Gmean,class.kk] = max(classG);
                    class.mdl = ClassificationKNN.fit(obj.trainD(TrInd,:),obj.trainG(TrInd),'NumNeighbors',class.kk);
                    out{k2,1} = class;
                end
        end
        function out = trainDecisionTree(obj)
            % devide into folds to extract training parameters
                F = getFolds(obj,obj.trainD,obj.trainG);
                out = cell(obj.Kfold,1);
                % for every fold train a classifier
                for k2=1:1:obj.Kfold
                    TestInd = find(F==k2);
                    TrInd = setdiff(1:size(obj.trainD,1),TestInd);
                    leafs = logspace(1,2,10);
                    leafs = round((leafs/100)*size(obj.trainD,2));
                    nr = numel(leafs);
                    classG = zeros(1,nr);
                    trees = cell(1,nr);
                    for l=1:1:nr
                        trees{l} = ClassificationTree.fit(obj.trainD(TrInd,:),obj.trainG(TrInd),'minleaf',leafs(l));
                        [~,~,~,bestlevel] = loss(trees{l},obj.trainD(TestInd,:),obj.trainG(TestInd),'SUBTREES','all');
                        trees{l} = prune(trees{l},'Level',bestlevel);    
                        pred.Classes = predict(trees{l},obj.trainD(TestInd,:));
                        testout = getEvaluation(obj,pred,obj.trainG(TestInd));
                        classG(l) = testout.G;
                    end
                    [class.Gmean,class.l] = max(classG);
                    class.Tree = trees{class.l};
                    out{k2,1} = class;
                    %view(class.Tree,'mode','graph');
                end
            
        end
        function out = trainLSVM(obj)
                F = getFolds(obj,obj.trainD,obj.trainG);
                out = cell(obj.Kfold,1);
                % for every fold train a classifier
                for k2=1:1:obj.Kfold
                    TestInd = find(F==k2);
                    TrInd = setdiff(1:size(obj.trainD,1),TestInd);
                    svmmodel = svmtrain(obj.trainD(TrInd,:),obj.trainG(TrInd),'Kernel_Function','linear','autoscale',false);
                    pred.Classes = svmclassify(svmmodel,obj.trainD(TestInd,:));
                    testout = getEvaluation(obj,pred,obj.trainG(TestInd));
                    class.svm = svmmodel;
                    class.Gmean = testout.G;
                    out{k2,1} = class;
                end
        end
        function out = trainKSVM(obj)
                F = getFolds(obj,obj.trainD,obj.trainG);
                out = cell(obj.Kfold,1);
                % for every fold train a classifier
                for k2=1:1:obj.Kfold
                    TestInd = find(F==k2);
                    TrInd = setdiff(1:size(obj.trainD,1),TestInd);
                    nrscales = 10;
                    scales = 1:nrscales;
                    classG = zeros(1,nrscales);
                    models = cell(1,nrscales);
                    for i=1:1:nrscales
                        models{i} = svmtrain(obj.trainD(TrInd,:),obj.trainG(TrInd),'Kernel_Function','rbf','autoscale',false,'rbf_sigma',scales(i));
                        pred.Classes = svmclassify(models{i},obj.trainD(TestInd,:));
                        testout = getEvaluation(obj,pred,obj.trainG(TestInd));
                        classG(i) = testout.G;
                    end
                    [class.Gmean,class.sigma] = max(classG);
                    class.svm = models{class.sigma};
                    out{k2,1} = class;
                end
        end
        function out = getEvaluation(obj,res,G)
            % by definition, minority is the positive group, and the
            % majority the negative group
            posest = res.Classes(G==obj.MinorLab);
            negest = res.Classes(G==obj.MajorLab);
            TP = length(find(posest==obj.MinorLab));
            FN = length(posest)-TP;
            TN = length(find(negest==obj.MajorLab));
            FP = length(negest)-TN;
            out.OA = (TP+TN)/length(G);
            out.Precision = TP/(TP+FP);
            out.Recall = TP/(TP+FN);
            out.G = sqrt((TP/(TP+FN))*(TN/(TN+FP)));
        end
        function out = classify(obj,D)
            % First perform necessary data transformations!
           switch obj.DataTransformation
               case 'none'
               case 'seuclidean'
                   D = ((obj.DW)*D')';
               case 'LMNN'
                   D = (obj.L*D')';
           end     
           out = TwoClassDemiModel.classifyWithModel(obj.BaggedClassifiers,D,obj.Classifier,obj.DRef);
        end
    end
    methods % implemented Balancing methods
        function imbalanced(obj,TrInd)
           obj.trainD = obj.transD(TrInd,:);
           obj.trainG = obj.G(TrInd,:);
           if ~obj.nrC==0, obj.trainC = obj.C(TrInd,:);end 
        end
        function ADASYN(obj,TrInd)
            D = obj.transD(TrInd,:);
            G = obj.G(TrInd,:);
            K = obj.Kadasyn;% number of closest examples to be found within ADASYN
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
            distances = pdist2(D(indm,:),D);
            [~,index] = sort(distances','ascend'); %#ok<UDIM>
            index = index(2:K+1,:);
            r = ismember(index,indM);
            r = sum(r)/K;
            r = r/sum(r);
            g = round(r*nrSyn);
            totg = sum(g);
            Dadd = zeros(totg,size(D,2));
            Gadd = obj.MinorLab*ones(totg,1);
            Dm = D(indm,:);
            if ~obj.nrC == 0, 
               Cm = C(indm,:);
               Cadd = zeros(totg,size(C,2));
            end
            counter = 0;
            distances = squareform(pdist(Dm));
            [~,index] = sort(distances,'ascend');
            index = index(2:K+1,:);
            for i=1:1:nrm
                if g(i)==0, continue; end
                for l=1:g(i)
                    counter = counter+1;
                    nb = randi(K,1,1);
                    alpha = rand(1);
                    Dadd(counter,:) = Dm(i,:) + (Dm(index(nb,i),:)-Dm(i,:))*alpha;
                    if ~obj.nrC==0, Cadd(counter,:) = Cm(i,:) + (Cm(index(nb,i),:)-Cm(i,:))*alpha;end
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
        function info = Tomec(obj)
            % remove all minor and major class members whose distance is
            % closer than any other co-class member
            D = obj.trainD;
            G = obj.trainG;
            if ~obj.nrC==0, C = obj.trainC; end
            indm = find(G==obj.MinorLab);
            indM = find(G==obj.MajorLab);
            Dm = D(indm,:);
            DM = D(indM,:);
            distances = pdist2(Dm,D);
            [sortdistances,index] = sort(distances','ascend');
            index = index(2,:);
            sortdistances = sortdistances(2,:);
            Minst = find(ismember(index,indM));
            sortdistances = sortdistances(Minst);
            DMsub = D(index(Minst),:);
            Mdistances = pdist2(DMsub,DM);
            [sortMdistances,~] = sort(Mdistances','ascend'); %#ok<*UDIM>
            sortMdistances = sortMdistances(2,:);
            Dmsub = Dm(Minst,:);
            mdistances = pdist2(Dmsub,Dm);
            [sortmdistances,~] = sort(mdistances','ascend'); %#ok<*UDIM>
            sortmdistances = sortmdistances(2,:);
            minst1 = find(sortdistances<sortMdistances);
            minst2 = find(sortdistances<sortmdistances);
            minst = intersect(minst1,minst2);
            % minor class instances to remove
            mremove = indm(Minst(minst));
            % major class instances to remove
            Mremove = index(Minst(minst));
            info = ((length(mremove)*2)/size(D,1))*100;
            keep = setdiff(1:size(D,1),mremove);
            keep = setdiff(keep,Mremove);
            obj.trainD = D(keep,:);
            obj.trainG = G(keep);
            if ~obj.nrC==0, obj.trainC = C(keep,:); end
        end
        function info = ENN(obj)
            % remove all minor and major class members whos is different to
            % 2 out of 3 closest instances
            D = obj.trainD;
            G = obj.trainG;
            if ~obj.nrC==0, C = obj.trainC; end
            indm = find(G==obj.MinorLab);
            distances = squareform(pdist(D));
            [~,index] = sort(distances','ascend');
            index = index(1:4,:);
            members = ismember(index,indm)';
            diff = abs(repmat(members(:,1),1,3)-members(:,2:end));
            diff = sum(diff,2);
            remove = find(diff>1);
            info = ((length(remove))/size(D,1))*100;
            keep = setdiff(1:size(D,1),remove);
            obj.trainD = D(keep,:);
            obj.trainG = G(keep);
            if ~obj.nrC==0, obj.trainC = C(keep,:); end
        end
    end
    methods (Static = true)
        function [Classes,Confidence] = classifyWithModel(model,D,method,DRef)
            model = model(:);
            nrM = numel(model);
            nrT = size(D,1);
            out = nan*zeros(nrT,nrM);
            %parfor i=1:nrM
            for i=1:nrM
                class = model{i};
                switch method
                    case {'PLSR RIP ROC' 'BRIM RIP ROC'}% simple thresholding on rip values
                        rip = TwoClassDemiModel.getRIP(DRef',D',class.M);
                        est = class.lower*ones(1,nrT);
                        est(rip>class.T) = class.higher;
                    case 'KNN'
                        est = predict(class.mdl,D);
                    case 'Decision Tree'
                        est = predict(class.Tree,D);
                    case {'LSVM' 'KSVM'}
                        est = svmclassify(class.svm,D);
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
        function out = extractInfo(info)
            info = info(:);
            nr = length(info);
            out.AvgRemoved = zeros(1,nr);
            out.StdRemoved = zeros(1,nr);
            for i=1:1:length(info)
                out.AvgRemoved(i) = mean(info{i}.Removed);
                out.StdRemoved(i) = std(info{i}.Removed);
            end
        end
        function out = getRIP(Ref,in,M)
            [~,n] = size(in); % determine input size
            out = nan*zeros(1,n);% allocate memory
            % distance between input faces and used reference
            Dir = repmat(Ref,1,n);
            Dir = (in-Dir);
            dist = sqrt(sum(((Dir).^2)));
            % Direction between input face and used reference
            coeff2 = M(end,:)';
            for i=1:1:n
                coeff1 = Dir(:,i)/norm(Dir(:,i));
                % Angle between direction and model direction
                T = coeff1'*coeff2;
                N = sqrt((coeff1'*coeff1)*(coeff2'*coeff2));
                angle = T/N;
                % parallell distance (RIP)
                out(i) = angle*dist(i);
            end   
        end
    end
end