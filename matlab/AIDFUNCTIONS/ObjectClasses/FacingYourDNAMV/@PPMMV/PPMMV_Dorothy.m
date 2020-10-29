classdef PPMMV <superClassLight
    % PROPERTIES
    properties % GENERAL INTERFACING
        nLev = 6;
    end
    properties (Dependent = true)
        nLC;
        nS;
        AvgX;
        StdX;
        PercNan;
    end
    properties (Hidden = true)
        Levels;
        Clusters;
    end
    properties % ASSOCIATION TESTING
        R2;
        pR2;
        A;
        pA;
        pAR;
        FA;
        pFA;
        EER;
        ppR2;           % created by Dorothy
        ppA;            % created by Dorothy
        ppAR;           % created by Dorothy
        
    end
    properties (Hidden = true)
        TpA = single(0.05);
        TpAR = single(0.05);
        TpFA = single(0.05);
    end
    properties (Hidden = true, Dependent = true)
        XMOD; % X used in model (excluding nan values in XREG)
        XMODInd; % Index of non nan XREG values
        MatrixR2;
        MatrixpR2;
        MatrixA;
        MatrixpA;
        MatrixpAR;
        MatrixFA;
        MatrixpFA;
        Hits;
        MatrixHits;
        nHits;
        LHits;
        wpA;
        wpAR;
        wpAF;
        wwpA;
        wwpAR;
        wwpAF;
        HMatrixpR2;              % created by Dorothy
        HMatrixpA;               % created by Dorothy
        HMatrixpAR;              % created by Dorothy
        HMatrixppR2;             % created by Dorothy
        HMatrixppA;              % created by Dorothy
        HMatrixppAR;             % created by Dorothy
    end
    properties % FEATURE LEARNING
        FeatureType = 'ANGLETEST'; % ANGLETEST, PLSR or BRIM
    end
    properties (Hidden = true)
        HTest = true;
        HTestp = single(0.0001);
        RIPNorm = true;
        VarInfo; % INFO needed for RIP normalization
        FeatM; % FEATURE REGRESSION PATH
        FeatH; % FEATURE REGRESSION HYPOTHESIS RESULT
        SamplePerc = single(1); % Sampling percentage for each regression run (ROBUST REGRESSION)
        SampleWR = true; % Sampling with or without replacement (ROBUST REGRESSION)
        RegSampling = true;% Bootstrap/None (ROBUST REGRESSION)
        RegBalance = true; % balance the data yes or no (ROBUST REGRESSION)
        RegBalanceTH = single(0.3); % desired balance factor (ROBUST REGRESSION)
        RegBalanceMethod = 'UpSample'; % UpSample/DownSample/ADASYN (SNPBRIM)
        BRIMFold = uint8(10); % number of  Folds (BRIM)
        MaxIterations = uint8(5); % Maximum number of BRIM iterations (BRIM)
        StopCorr = single(0.98); % Stopping correlation between subsequent iterations (BRIM)
        ShapeDistance = 'seuclidean';
    end
    properties (Hidden = true, Dependent = true)
        DIM;
        MatrixDIM;
        RegRuns; % Number of Regression runs (ROBUST REGRESSION)
    end
    properties % PHENOMIC PREDICTOR
        Predict = 'SOFT';
    end
    properties (Hidden = true)
        FSelection = 'COR'; % COR, ROC, MRMTR, SPECCMI, CMIM, CIFE, SVMRFE
        FSInd;
        FSInfo;
        maxFS = 200;
        Quantiles = 5;
        nrTrees = uint8(100);
    end
    properties(Hidden = true, Dependent = true)
        nFS;
        nrSplits;
    end
    properties (Abstract = true)
        X;
        XREG; % X USED FOR REGRESSION/ SUBCLASS DEPENDENT
        Predictor;
        treeTempl;
        predictortype;
    end
    % OBJECT METHODS
    methods % CONSTRUCTOR
        function obj = PPMMV(varargin)
            obj = obj@superClassLight(varargin{:});
            storeLevels(obj);
            storeClusters(obj);
        end
    end
    methods % GENERAL GETTING
        function out = get.nS(obj)
            out = length(obj.X);
        end
        function out = get.XMODInd(obj)
            out = find(~isnan(obj.XREG));
        end
        function out = get.XMOD(obj)
            if isempty(obj.X), out = []; return; end
            out = double(obj.XREG(obj.XMODInd));
        end
        function out = get.PercNan(obj)
            if isempty(obj.X), out = []; return; end
            out = round((sum(isnan(obj.X))/obj.nS)*100);
        end
        function out = get.DIM(obj)
            if isempty(obj.FeatM), out = []; return; end
            out = nan*zeros(1,obj.nLC);
            for i=1:1:obj.nLC
                out(i) =  length(obj.FeatM{i});
            end
        end
        function out = get.nLC(obj)
            out = length(obj.Levels);
        end
        function out = get.MatrixR2(obj)
            out = List2Matrix(obj,obj.R2);
        end
        function out = get.MatrixpR2(obj)
            out = List2Matrix(obj,obj.pR2);
        end
        function out = get.MatrixA(obj)
            out = List2Matrix(obj,obj.A);
        end
        function out = get.MatrixFA(obj)
            out = List2Matrix(obj,obj.FA);
        end
        function out = get.MatrixpFA(obj)
            out = List2Matrix(obj,obj.pFA);
        end
        function out = get.MatrixpA(obj)
            out = List2Matrix(obj,obj.pA);
        end
        function out = get.MatrixpAR(obj)
            out = List2Matrix(obj,obj.pAR);
        end
        function out = get.MatrixDIM(obj)
            out = List2Matrix(obj,obj.DIM);
        end
        function out = get.Hits(obj)
            if isempty(obj.A), out = []; return;end
            out = (obj.pA<=obj.TpA).*(obj.pAR<=obj.TpAR);
        end
        function out = get.MatrixHits(obj)
            out = List2Matrix(obj,obj.Hits);
        end
        function out = get.nHits(obj)
            if isempty(obj.A), out = []; return; end
            out = sum(obj.Hits);
        end
        function out = get.LHits(obj)
            if isempty(obj.A), out = []; return; end
            out = nansum(obj.MHits,1);
        end
        function out = get.nFS(obj)
            out = length(obj.FSInd);
        end
        function out = get.nrSplits(obj)
            if obj.nFS==0, out = 0; return; end
            if obj.nFS < 50;
                out = 2;
            elseif obj.nFS < 200;
                out = 3;
            else
                out = 5;
            end
        end
        function out = get.AvgX(obj)
            out = nanmean(obj.X);
        end
        function out = get.StdX(obj)
            out = nanstd(obj.X);
        end
        function out = get.RegRuns(obj)
            switch obj.FeatureType
                case 'PLSR'
                    out = 1000;
                case 'BRIM'
                    out = 50;
                otherwise
                    out = 1;
            end
        end
        function out = get.RegBalance(obj)
            if strcmp(obj.predictortype,'regressor'), out = false; return; end
            out = obj.RegBalance;
        end
        function out = get.wpA(obj)
            if isempty(obj.A), out = 0; return; end
            out = nanmean(-log(obj.pA));
        end
        function out = get.wpAR(obj)
            if isempty(obj.A), out = 0; return; end
            out = nanmean(-log(obj.pAR));
        end
        function out = get.wpAF(obj)
            if isempty(obj.A), out = 0; return; end
            out = obj.wpA+obj.wpAR;
        end
        function out = get.wwpA(obj)
            if isempty(obj.A), out = 0; return; end
            %out = nanmean(-log(obj.pA));
            val = -log(obj.pA);
            index = find(~isnan(val));
            val = val(index);
            w = 1./double(obj.Levels(index));
            out = sum(w.*val)/sum(w);
            
        end
        function out = get.wwpAR(obj)
            if isempty(obj.A), out = 0; return; end
            val = -log(obj.pAR);
            index = find(~isnan(val));
            val = val(index);
            w = 1./double(obj.Levels(index));
            out = sum(w.*val)/sum(w);
        end
        function out = get.wwpAF(obj)
            if isempty(obj.A), out = 0; return; end
            out = obj.wwpA+obj.wwpAR;
        end
        function out = get.HMatrixpR2(obj)       % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.pR2);
        end
        function out = get.HMatrixpA(obj)        % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.pA);
        end
        function out = get.HMatrixpAR(obj)       % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.pAR);
        end
        function out = get.HMatrixppR2(obj)      % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.ppR2);
        end
        function out = get.HMatrixppA(obj)       % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.ppA);
        end
        function out = get.HMatrixppAR(obj)      % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.ppAR);
        end
    end
    methods % GENERAL SETTING
        function obj = set.nLev(obj,in)
            obj.nLev = single(in);
            storeLevels(obj);
            storeClusters(obj);
        end
        function obj = set.FSInd(obj,in)
            obj.FSInd = PPMMV.convertUInt(in);
        end
        function obj = set.maxFS(obj,in)
            obj.maxFS = PPMMV.convertUInt(in);
        end
        function obj = set.nrTrees(obj,in)
            obj.nrTrees = PPMMV.convertUInt(in);
        end
    end
    methods % ASSOCIATION TESTING
        function out = runR2Test(obj,DepVar,COV,t)
            R2 = nan*zeros(1,obj.nLC);
            pR2 = nan*zeros(1,obj.nLC);
            parfor i=1:obj.nLC
                [R2(i),pR2(i)] =  PPMMV.testPartialPLSR(obj.XREG,getDepVar(obj,DepVar,i),getCOV(obj,COV,i),t); %#ok<*PFBNS>
            end
            out.R2 = R2;out.pR2 = pR2;
            obj.R2 = R2;obj.pR2 = pR2;
        end
        function out = runAngleTest(obj,DepVar,COV,maxM,AngleTable,t)
            out = cell(1,obj.nLC);
            parfor i=1:obj.nLC
                [out{i}.A,out{i}.M,out{i}.STAT] = PPMMV.AngleTest([getCOV(obj,COV,i) obj.XREG],getDepVar(obj,DepVar,i),AngleTable,t,maxM,0);
            end
            obj.A = nan*zeros(1,obj.nLC);
            obj.pA = nan*zeros(1,obj.nLC);
            obj.pAR = nan*zeros(1,obj.nLC);
            obj.FA = nan*zeros(1,obj.nLC);
            obj.pFA = nan*zeros(1,obj.nLC);
            for i=1:1:obj.nLC
                obj.A(i) = out{i}.STAT.AvgA;
                obj.pA(i) = out{i}.STAT.pA;
                obj.pAR(i) = out{i}.STAT.pAR;
                obj.FA(i) = out{i}.STAT.FA;
                obj.pFA(i) = out{i}.STAT.pFA;
            end
        end
        function out = runPropagatingR2Test(obj,DepVar,COV,t)        % created by Dorothy
            ppR2 = nan*zeros(1,obj.nLC);
            ppR2Perms = nan*zeros(t,obj.nLC);
            nS = size(DepVar,1);
            forcedPerms = getForcedPermutations(t,nS);
            
            % level with no children
            l = obj.nLev;
            for cl = 1:1:2^(l-1),
                clind = getClind(obj,l,cl);
                [ppR2(clind), ppR2Perms(:,clind)] = testPartialPLSR_Permutated(obj.XREG,getDepVar(obj,DepVar,clind),getCOV(obj,COV,clind),t,forcedPerms);
            end
            
            % levels with children
            for l = obj.nLev-1:-1:1,
                for cl = 1:1:2^(l-1),
                    % get the p values for the parent cluster
                    clind = getClind(obj,l,cl);
                    [ppR2_temp, ppR2Perms_temp] = testPartialPLSR_Permutated(obj.XREG,getDepVar(obj,DepVar,clind),getCOV(obj,COV,clind),t,forcedPerms);
                    
                    % get the p values for the two children clusters
                    children = getChildren(obj,l,cl);
                    clind1 = getClind(obj,l+1,children(1));
                    clind2 = getClind(obj,l+1,children(2));
                    
                    % nonparametric combination of the p values of the parent and the two children
                    [ppR2(clind), ppR2Perms(:,clind)] = nonParametricCombination ([ppR2_temp,ppR2(clind1),ppR2(clind2)], [ppR2Perms_temp,ppR2Perms(:,clind1),ppR2Perms(:,clind2)]); 
                end
            end
            out.ppR2 = ppR2;
            obj.ppR2 = ppR2;
            
        end
        function out = runPropagatingAngleTest(obj,DepVar,COV,maxM,AngleTable,t,t2)    % created by Dorothy (long computation time)
            ppA = nan*zeros(1,obj.nLC);
            ppAPerms = nan*zeros(t2,obj.nLC);
            ppAR = nan*zeros(1,obj.nLC);
            ppARPerms = nan*zeros(t2,obj.nLC);
            nS = size(DepVar,1);
            forcedPerms = getForcedPermutations(t2,nS);
            
            % level with no children
            l = obj.nLev;
            for cl = 1:1:2^(l-1),
                clind = getClind(obj,l,cl);
                [ppA(clind), ppAPerms(:,clind), ppAR(clind), ppARPerms(:,clind)] = AngleTest_Permutated([getCOV(obj,COV,clind) obj.XREG],getDepVar(obj,DepVar,clind),AngleTable,t,maxM,0,t2,forcedPerms);
            end
            
            % levels with children
            for l = obj.nLev-1:-1:1,
                for cl = 1:1:2^(l-1),
                    % get the p values for the parent cluster
                    clind = getClind(obj,l,cl);
                    [ppA_temp, ppAPerms_temp, ppAR_temp, ppARPerms_temp] = AngleTest_Permutated([getCOV(obj,COV,clind) obj.XREG],getDepVar(obj,DepVar,clind),AngleTable,t,maxM,0,t2,forcedPerms);
            
                    % get the p values for the two children clusters
                    children = getChildren(obj,l,cl);
                    clind1 = getClind(obj,l+1,children(1));
                    clind2 = getClind(obj,l+1,children(2));
                    
                    % nonparametric combination of the p values of the parent and the two children
                    [ppA(clind), ppAPerms(:,clind)] = nonParametricCombination ([ppA_temp,ppA(clind1),ppA(clind2)], [ppAPerms_temp,ppAPerms(:,clind1),ppAPerms(:,clind2)]); 
                    [ppAR(clind), ppARPerms(:,clind)] = nonParametricCombination ([ppAR_temp,ppAR(clind1),ppAR(clind2)], [ppARPerms_temp,ppARPerms(:,clind1),ppARPerms(:,clind2)]); 
                end
            end
            out.ppA = ppA;
            obj.ppA = ppA;
            out.ppAR = ppAR;
            obj.ppAR = ppAR;
        end
    end
    methods % FEATURE LEARNING
        function featureLearning(obj,DepVar,varargin)
            %DepVar = TrDepVar;
            info = cell(1,obj.nLC);
            M = cell(1,obj.nLC);
            H = cell(1,obj.nLC);
            switch obj.FeatureType
                case 'ANGLETEST'
                    ATest = varargin{1};
                    COV = [];
                case {'PLSR' 'BRIM'}
                    COV = varargin{1};
                    ATest = cell(1,obj.nLC);
            end
            parfor i=1:obj.nLC
                forDepVar = getDepVar(obj,DepVar,i);
                % VARIABLE INFORMATION
                tmpinfo = getVarInfo(obj,forDepVar);
                info{i} = tmpinfo;
                % M CONSTRUCTION
                switch obj.FeatureType
                    case 'ANGLETEST'
                        [M{i},H{i}] = learnFromAngleTest(obj,ATest{i});
                    case 'PLSR'
                        [M{i},H{i}] = learnFromPLSR(obj,forDepVar,getCOV(obj,COV,i));
                    case 'BRIM'
                        [M{i},H{i}] = learnFromBRIM(obj,forDepVar,getCOV(obj,COV,i),tmpinfo);
                end
            end
            obj.VarInfo = info;
            obj.FeatM = M;
            obj.FeatH = H;
            postFeatureLearning(obj);
        end
        function [M,H] = learnFromAngleTest(obj,ATest)
            % M CONSTRUCTION
            fullM = [squeeze(ATest.M(:,1,:)) squeeze(ATest.M(:,2,:))];
            M = nanmedian(fullM,2)';
            H = PPMMV.wilcoxonM(fullM',obj.HTestp);
        end
        function [M,H] = learnFromPLSR(obj,DepVar,COV)
            % INITIALIZATION
            GG = obj.XMOD;
            switch obj.predictortype
                case 'classifier'
                    GM = GG;
                case 'regressor'
                    GM = [];
            end
            DepVar = DepVar(obj.XMODInd,:);
            if isempty(COV)
                IndVar = GG;
            else
                IndVar = [COV(obj.XMODInd,:) GG];
            end
            % M CONSTRUCTION
            [M,H] = robustPLSRegression(obj,IndVar,DepVar,GM,obj.RegRuns);
        end
        function [M,H] = learnFromBRIM(obj,DepVar,COV,info)
            % DepVar = getDepVar(obj,TrDepVar,1);info = getVarInfo(obj,DepVar);
            % To bootstrap or not to bootstrap that is the question
            Bootstrap = true;if obj.MaxIterations == 0, Bootstrap = false;end
            % INITIALIZATION
            GG = obj.XMOD;
            switch obj.predictortype
                case 'classifier'
                    GM = GG;
                case 'regressor'
                    GM = ones(size(GG));
            end
            DepVar = DepVar(obj.XMODInd,:);
            [n,~] = size(DepVar);
            if isempty(COV)
                IndVar = GG;
            else
                IndVar = [COV(obj.XMODInd,:) GG];
            end
            % BRIM
            BootProgress = zeros(1,obj.MaxIterations);
            ContBoot = true;counter = 0;
            while ContBoot && Bootstrap>0
                % Keeping track of iterations
                counter = counter + 1;
                % seperate data into inner folds
                Finner = PPMMV.DOB_SCV_DM(obj.BRIMFold,DepVar,GM,obj.ShapeDistance);
                TMPIndVar = IndVar(:,end);% Allocate Memory, only last collumn
                for fi=1:obj.BRIMFold % Fi...
                    FiTestInd = find(Finner==fi);
                    FiTrInd = setdiff(1:n,FiTestInd);
                    [M,H] = robustPLSRegression(obj,IndVar(FiTrInd,:),DepVar(FiTrInd,:),GM(FiTrInd),obj.RegRuns);% Compute BootSampled Regression
                    if obj.HTest, M = H.*M;end% Only keep significant regression coefficients
                    rip = PPMMV.updateRIP(DepVar(FiTestInd,:),M)';% Get RIP scores
                    if obj.RIPNorm, rip = PPMMV.normalizeRIP(rip,M,info);end% rescale RIP scores
                    TMPIndVar(FiTestInd) = rip;% Store updated RIP scores, always in the last collumn
                end
                % monitoring progress of Booting
                tmpC = corrcoef(IndVar(:,end),TMPIndVar);
                BootProgress(counter) = tmpC(1,2);
                if BootProgress(counter)>=obj.StopCorr, ContBoot = false; end
                if counter >= obj.MaxIterations, ContBoot = false; end
                IndVar(:,end) = TMPIndVar;% Update Training IndVar for next round
            end
            % Extracting final M from outer training data
            [M,H] = robustPLSRegression(obj,IndVar,DepVar,GM,obj.RegRuns*20);
        end
        function [M,H] = robustPLSRegression(obj,indvar,depvar,gg,runs)
            % initialize
            sampling = obj.RegSampling;balance = obj.RegBalance;% only read once
            if runs==1, sampling = false; end% in a single run no point in sampling
            if balance % determine group memberships and balances, within test such that not executed when not necessary
                balanceTH = obj.RegBalanceTH;balanceMethod = obj.RegBalanceMethod;
                labgg = unique(gg);nrgg = length(labgg);% group labels
                indgg = cell(1,nrgg);nrwgg = zeros(1,nrgg);
                for i=1:1:nrgg
                    indgg{i} = find(gg==labgg(i));
                    nrwgg(i) = length(indgg{i});
                end
                [nrM,indM] = max(nrwgg);% Majority class
                indm = setdiff(1:nrgg,indM);nrm = nrwgg(indm);% minority class(es)
                nrmgg = length(indm);dcrit = zeros(1,nrmgg);nrSyn = zeros(1,nrmgg);nrKeep = nrM;% Initialze minority classes
                for i=1:1:nrmgg
                    dcrit(i) = nrm(i)/nrM;
                    if dcrit(i)>=balanceTH; continue; end
                    nrSyn(i) = (nrM-nrm(i))*balanceTH;% UpSampling
                    nrKeep = min(nrKeep,round(nrm(i)/balanceTH));% DownSampling
                end
                balInfo = cell(1,nrgg);
                for i=1:1:nrgg % gathering information to be used during runs
                    Info.Ind = indgg{i};
                    switch balanceMethod
                        case 'DownSample'
                            Info.nrKeep = nrKeep;
                            if i==1, upsampling = false; end
                        case {'UpSample' 'ADASYN'}
                            if i==1, upsampling = true; end
                            if i==indM
                                Info.nrSyn = 0;
                            else
                                K = 10;% number of closest neighbors
                                subind = find(indm==i);
                                Info.nrSyn = round(nrSyn(subind));
                                if Info.nrSyn>0
                                    tmpindm = indgg{indm(subind)};
                                    % within group K closest neighbors
                                    distances = squareform(pdist(depvar(tmpindm,:),obj.ShapeDistance));
                                    [~,index] = sort(distances,'ascend');
                                    mK = min(length(tmpindm)-1,K);% it is possible that you do not have enough samples
                                    index = index(2:mK+1,:);
                                    Info.SampleIndB = index;
                                    Info.K = mK;
                                    % establishing sampling distribution
                                    if strcmp(balanceMethod,'UpSample')
                                        % all within sample are equally likely to be chosen
                                        Info.SampleInd = 1:length(tmpindm);
                                    else
                                        % samples closer to other group members get higher weigths (ADASYN)
                                        tmpindM = indgg{indM};
                                        for j=1:1:nrmgg
                                            if j==subind; continue; end
                                            tmpindM = [tmpindM; indgg{indm(j)}];
                                        end
                                        % between group K closest neighbors
                                        distances = pdist2(depvar(tmpindm,:),depvar);
                                        [~,index] = sort(distances','ascend'); %#ok<UDIM>
                                        index = index(2:K+1,:); % probably no K problem
                                        r = ismember(index,tmpindM);r = sum(r)/K;r = r/sum(r);
                                        g = round(r*nrSyn(subind));totg = sum(g);
                                        if totg==0, g = ones(1,length(tmpindm)); totg = length(tmpindm); end % perfectly seperatable groups
                                        index = zeros(1,totg);
                                        counter = 1;
                                        for j=1:1:length(g)
                                            if g(j)==0, continue; end
                                            for k=1:1:g(j)
                                                index(counter) = j;
                                                counter = counter +1;
                                            end
                                        end
                                        Info.SampleInd = index;
                                    end
                                end
                            end
                    end
                    balInfo{i} = Info;
                end
            end
            nrD = size(depvar,2);Mfor = zeros(runs,nrD);
            % execute runs
            for s=1:1:runs
                if balance % Balance the given input
                    if upsampling % Balancing by Increasing Minority Groups
                        SampleIndVar = indvar;SampleDepVar = depvar;
                        for i=1:1:nrgg
                            if balInfo{i}.nrSyn==0, continue; end
                            [DepVarNew,IndVarNew] = PPMMV.regUpSample(depvar(balInfo{i}.Ind,:),indvar(balInfo{i}.Ind,:),balInfo{i});
                            SampleDepVar = [SampleDepVar; DepVarNew]; %#ok<*AGROW>
                            SampleIndVar = [SampleIndVar; IndVarNew];
                        end
                    else % Balancing by Decreasing Majority Groups
                        SampleIndVar = [];SampleDepVar = [];
                        for i=1:1:nrgg
                            [DepVarNew,IndVarNew] = PPMMV.regDownSample(depvar(balInfo{i}.Ind,:),indvar(balInfo{i}.Ind,:),balInfo{i}.nrKeep);
                            SampleDepVar = [SampleDepVar; DepVarNew];
                            SampleIndVar = [SampleIndVar; IndVarNew];
                        end
                    end
                else % No Balancing
                    SampleIndVar = indvar;SampleDepVar = depvar;
                end
                nrS = size(SampleIndVar,1);
                if sampling % Bootstrap Sampling of given input
                    SampleInd = randsample(nrS,round(obj.SamplePerc*nrS),obj.SampleWR);
                else
                    SampleInd = 1:nrS;
                end
                Mfor(s,:) = PPMMV.getRegression(SampleIndVar(SampleInd,:),SampleDepVar(SampleInd,:));
            end
            % finalizing output
            M = median(Mfor,1); % extract the median of all regression coefficients
            if ~obj.HTest, H=ones(1,nrD);return; end
            if runs==1, H=ones(1,nrD); return; end
            if nargout<2, return; end
            H = PPMMV.wilcoxonM(Mfor,obj.HTestp);% wilcoxon test for median
        end
        function out = getRIP(obj,DepVar,clind)
            if nargin<3 % FULL RIP CALCULATION
                out = nan*zeros(size(DepVar{1}.DepVar{1},1),obj.nLC);
                for i=1:obj.nLC
                    out(:,i) = getRIP(obj,DepVar,i);
                end
                return;
            end
            if ~isscalar(clind)% SELECTED RIP CALCULATION
                out = nan*zeros(size(DepVar{1}.DepVar{1},1),length(clind));
                for i=1:1:length(clind)
                    out(:,i) = getRIP(obj,DepVar,clind(i));
                end
                return;
            end
            % SINGLE RIP CALCULATION
            out = PPMMV.updateRIP(getDepVar(obj,DepVar,clind),getFeatM(obj,clind));
            if obj.RIPNorm, out = PPMMV.normalizeRIP(out,getFeatM(obj,clind),obj.VarInfo{clind});end
            out = out';
        end
        function out = getFeatM(obj,clind)
            out = obj.FeatM{clind};
            if obj.HTest, out = out.*obj.FeatH{clind};end
        end
    end
    methods % PHENOMIC PREDICTOR
        function out = featureSelection(obj,X,XC,FS)
            rng(1);
            switch obj.FSelection
                case 'COR'
                    [out,obj.FSInfo] = PPMMV.fsCOR(FS,XC,obj.maxFS);
                case 'ROC'
                    [out,obj.FSInfo] = PPMMV.fsROC(FS,X,obj.maxFS);
                case 'MRMTR'
                    [~,~,H] = PPMMV.prepMIfs(FS,X,obj.Quantiles);
                    out = PPMMV.fsMRMTR(H,min(obj.maxFS,size(FS,2)));
                    obj.FSInfo = [];
                case 'SPECCMI'
                    [~,~,H] = PPMMV.prepMIfs(FS,X,obj.Quantiles);
                    [out,obj.FSInfo] = PPMMV.fsSPECCMI(H,min(obj.maxFS,size(FS,2)));
                case 'CMIM'
                    [~,~,H] = PPMMV.prepMIfs(FS,X,obj.Quantiles);
                    out = PPMMV.fsCMIM(H,min(obj.maxFS,size(FS,2)));
                    obj.FSInfo = [];
                case 'CIFE'
                    [~,~,H] = PPMMV.prepMIfs(FS,X,obj.Quantiles);
                    out = PPMMV.fsCIFE(H,min(obj.maxFS,size(FS,2)));
                    obj.FSInfo = [];
                case 'MIXED'
                    nFS = size(FS,2);
                    ranks = zeros(3,nFS);
                    info = zeros(3,nFS);
                    [ranks(1,:),info(1,:)] = PPMMV.fsCOR(FS,XC,nFS);
                    [ranks(2,:),info(2,:)] = PPMMV.fsROC(FS,X,nFS);
                    [~,~,H] = PPMMV.prepMIfs(FS,X,obj.Quantiles);
                    [ranks(3,:),info(3,:)] = PPMMV.fsSPECCMI(H,nFS);
                    r = 1:nFS;
                    new = repmat(r,3,1);
                    for i=1:1:3
                        new(i,ranks(i,:)) = r;
                    end
                    test = mean(new,1);
                    [~,out] = sort(test,'ascend');
                    out = out(1:obj.maxFS);
                    obj.FSInfo = info;
                case 'SVMRFE'
                case 'NONE'
                    out = 1:size(FS,2);
                    obj.FSInfo = [];
                otherwise
                    error('Wrong feature selection procedure')
            end
            obj.FSInd = out;
        end
        function out = selectFeatures(obj,FS)
            if isempty(obj.FSInd), out = []; return; end
            out = FS(:,obj.FSInd);
        end
    end
    methods % BIOMETRICS
    end
    methods % INTERFACING
        function storeLevels(obj)
            out = [];
            for i=1:1:obj.nLev
                out = [out i*ones(1,2^(i-1))];  %#ok<AGROW>
            end
            obj.Levels = PPMMV.convertUInt(out);
        end
        function storeClusters(obj)
            out = [];
            for i=1:1:obj.nLev
                out = [out 1:2^(i-1)]; %#ok<AGROW>
            end
            obj.Clusters = PPMMV.convertUInt(out);
        end
        function out = getDepVar(obj,in,clind)
            if ~iscell(in), out = in;return;end
            if ~isscalar(clind)
                out = [];
                for i=1:1:length(clind)
                    out = [out, getDepVar(obj,in,clind(i))];
                end
                return;
            end
            if ~clind==0, out = in{obj.Levels(clind)}.DepVar{obj.Clusters(clind)};return;end
            out = [];
            for l=1:1:length(in)
                for i=1:1:length(in{l}.DepVar)
                    out = [out, in{l}.DepVar{i}]; %#ok<AGROW>
                end
            end
        end
        function out = getFS(obj,FS)
            if ~iscell(FS), out = FS; return; end
            out = getDepVar(obj,FS,find(obj.Hits));
        end
        function out = getCOV(obj,in,clind)
            if ~iscell(in), out = in;return;end
            if ~clind==0, out = in{obj.Levels(clind)}.COV{obj.Clusters(clind)};return;end
            out = [];
            for l=1:1:length(in)
                for i=1:1:length(in{l}.COV)
                    out = [out, in{l}.COV{i}]; %#ok<AGROW>
                end
            end
        end
        function out = getSymSpace(obj,in,clind)
            if ~iscell(in), out = in;return;end
            if ~obj.Level==0,out = in{obj.Levels(clind)}.SymSpace{obj.Clusters(clind)};return;end
            out = [];% cannot retrieve a single SymSpace for level 0;
        end
        function out = List2Matrix(obj,in)
            if isempty(in), out = []; return; end
            out = nan*zeros(obj.nLev,max(obj.Clusters));
            %out(obj.Levels,obj.Clusters) = in;
            for i=1:obj.nLC
                out(obj.Levels(i),obj.Clusters(i)) = in(i);
            end
            out = out';
        end
        function out = getTestX(obj,COV,GB,GT,RS,nrT)
            switch obj.ID
                case {1 2 3 4 5 6 7 8 9 10}
                    out = COV(:,obj.ID);
                case 100
                    out = GB(:,obj.GBID);
                case 1000
                    ind = find(strcmp(obj.RS,RS));
                    if isempty(ind),out = nan*zeros(nrT,1);return;end
                    out = GT(:,ind);
                otherwise
                    out = nan*zeros(nrT,1);
            end
        end
        function out = List2HierachicalMatrix(obj,in,significance)    % created by Dorothy
            if nargin < 3 || isempty(significance)==1, significance = true; end
            if islogical(significance)==0, out = []; return; end
            nCl = 2^(obj.nLev-1);
            out = zeros(obj.nLev,nCl);
            for cl = 1:nCl,
                parentEvolution = getParentEvolution(obj,obj.nLev,cl);
                for l = 1:obj.nLev,
                    clind = getClind(obj,l,parentEvolution(l));
                    value = in(clind);
                    if significance == false;        % fill the matrix with the value
                        out(l,cl) = value;
                    else                             % fill the matrix with 1 of the value is significant
                        p_crit = 0.05;
                        if value <= p_crit;
                            out(l,cl) = 1;
                        end
                    end
                end
            end
            % Plot the figure with imagesc(out);
        end
        
    end
    methods (Abstract = true)
        out = getVarInfo(obj,DepVar);
        postFeatureLearning(obj);
        trainPredictor(obj,X,REGF,PCF,ADDF);
        out = testPredictor(obj,X,REGF,PCF,ADDF);
        [H,S] = predict(obj,REGF,PCF,ADDF);
    end
    % STATIC METHODS
    methods (Static = true) % ASSOCIATION ANALYSIS
        function [R2,pR2] = testPartialPLSR(X,Y,C,t)
            index = intersect(PPMMV.notNAN(C),PPMMV.notNAN(X));
            if ~isempty(C)
                E = PPMMV.getResiduals(C(index,:),Y(index,:));
                X = PPMMV.getResiduals(C(index,:),X(index,:));
            else
                E = Y(index,:);
                X = X(index,:);
            end
            [~,~,~,~,~,var] = plsregress(X,E,1);R2 = var(2);
            if t==0, pR2 = nan; return; end
            R2Count = false(t,1);
            % incremental permutation
            for p=1:1:t
                pind = randperm(length(index));
                [~,~,~,~,~,forvar] = plsregress(X,E(pind,:),1);
                R2Count(p) = forvar(2)>=R2;
                if mod(p,100)==0
                    pR2 = (sum(R2Count(1:p))+1)/(p+1);
                    acc = 10/p;
                    if pR2>acc, return;end
                end
            end
            pR2 = (sum(R2Count)+1)/(t+1);
        end
        function [A,M,STAT] = AngleTest(IndVar,DepVar,AngleTable,t,maxM,Aval)
            if nargin<6,Aval=0;end
            if nargin<5,maxM=t;end
            A = nan*zeros(t,1);
            [nS,ShapeDim] = size(DepVar);
            M = nan*zeros(ShapeDim,2,maxM);
            el = 0;
            for i=1:t
                el = el+1;
                F = crossvalind('Kfold',nS,2);
                forM = nan*zeros(ShapeDim,2);
                for k=1:1:2
                    ind = find(F==k);
                    forM(:,k) = PPMMV.getRegression(IndVar(ind,:),DepVar(ind,:));
                end
                A(i) = angle(forM(:,1),forM(:,2));
                M(:,:,el) = forM;
                if mod(i,maxM)==0, el = 0; end
                if mod(i,50)==0
                    tmppA = (length(find(A(1:i)<=0))+1)/(i+1);
                    acc = 10/i;
                    if tmppA>acc, break;end
                end
            end
            A = A(1:i);
            index = find(~isnan(A));
            A = A(index);
            i = length(A);
            STAT.AvgA = nanmean(A);
            STAT.pA = zeros(length(Aval),1);
            for j=1:1:length(Aval)
                STAT.pA(j) = (length(find(A<=Aval(j)))+1)/(i+1);
            end
            STAT.pAR = PPMMV.lookUppRA(STAT.AvgA,ShapeDim,AngleTable);
            [STAT.FA,STAT.pFA] = PPMMV.lookUpFA(A,ShapeDim,AngleTable);
            STAT.MedA = median(A);
            STAT.StdA = std(A);
            STAT.MadA = mad(A);
            sortA = sort(A);
            STAT.UpperA = sortA(round(0.975*i));
            STAT.LowerA = sortA(round(0.025*i+1));
            if i<maxM,M = M(:,:,index);end
        end
        function M = getRegression(A,B)
            [A,B] = eliminateNAN(A,B);
            [~,~,~,~,M] = plsregress(A,B,size(A,2));
            M = M(end,:);
        end
        function out = getResiduals(X,Y)
            [~,~,~,~,~,~,~,stats] = plsregress(X,Y,min(size(X,2),size(Y,2)));
            out = stats.Yresiduals;
        end
        function out = lookUppRA(A,Dim,Table)
            out = interp1(Table.Angles,Table.pAngles(Dim,:),A);
        end
        function [FA,pFA] = lookUpFA(A,Dim,Table)
            nA = length(A);
            G = [ones(nA,1);zeros(size(Table.Distr,2),1)];
            [~,T] = anova1([A;Table.Distr(Dim,:)'],G,'off');
            FA = T{2,5};
            pFA  = T{2,6};
        end
        function [ppR2, ppR2Perms] = testPartialPLSR_Permutated(X,Y,C,t,forcedPerms)        % Created by Dorothy
            index = intersect(PPMMV.notNAN(C),PPMMV.notNAN(X));
            if ~isempty(C)
                E = PPMMV.getResiduals(C(index,:),Y(index,:));
                X = PPMMV.getResiduals(C(index,:),X(index,:));
            else
                E = Y(index,:);
                X = X(index,:);
            end
            
            % calculate R2 for the true data
            [~,~,~,~,~,var] = plsregress(X,E,1);R2True = var(2);
            if t==0, ppR2 = nan; return; end
            
            % calculate R2 for the permutated data
            R2Perm = zeros(t,1);
            for p=1:1:t
                if nargin < 5 || isempty(forcedPerms)==1;
                    perm = randperm(length(index));
                else
                    perm = forcedPerms(p,:);
                end
                [~,~,~,~,~,forvar] = plsregress(X,E(perm,:),1);
                R2Perm(p) = forvar(p);
            end
            
            % calculate the p-value for the true data and the permutated data
            ppR2 = sum([R2Perm;R2True] >= repmat(R2True,t+1,1),1)/(t+1);
            ppR2Perms = sum(repmat([R2Perm;R2True]',t,1) >= repmat(R2Perm,1,t+1),2)./(t+1);
        end
        function [pA, pAPerm, pAR, pARPerm] = AngleTest_Permutated(IndVar,DepVar,AngleTable,t,maxM,Aval,t2,forcedPerms)    % Created by Dorothy (long computation time)
            % calculate the p-value for the true data
            [~,~,STAT] = AngleTest(IndVar,DepVar,AngleTable,t,maxM,Aval);
            pA = STAT.pA;
            pAR = STAT.pAR;
            
            % calculate the p-value for the permutated data
            pAPerm = zeros(t2,1);
            pARPerm = zeros(t2,1);
            for p = 1:1:t2; 
                if nargin < 8 || isempty(forcedPerms)==1;
                    perm = randperm(size(IndVar,1));
                else
                    perm = forcedPerms(p,:);
                end
                [~,~,STAT] = AngleTest(IndVar(perm,:),DepVar,AngleTable,t,maxM,Aval);
                pAPerm(p) = STAT.pA;
                pARPerm(p) = STAT.pAR;
            end
        end
    end
    methods (Static = true) % FEATURE LEARNING
        function H = wilcoxonM(M,pT)
            nrD = size(M,2);
            P = zeros(1,nrD);
            for j=1:1:nrD
                P(j) = signrank(M(:,j));% wilcoxon test for median
            end
            H = logical(P<=pT);
        end
        function [out] = normalizeRIP(rip,M,var)
            Mrip = dot(var.MDepVar',M'/norm(M'));
            Prip = dot(var.PDepVar',M'/norm(M'));
            out = ((rip-Mrip)/(Prip-Mrip))*var.Range+var.Mel;
        end
        function [out] = updateRIP(in,M)
            out = dot(in',repmat(M'/norm(M'),1,size(in,1)));
        end
        function [GRed,CRed] = regDownSample(GO,CO,nrKeep)
            nr = size(GO,1);
            if nr>nrKeep % randomly reduce
                s = randsample(nr,nrKeep,false);
                GRed = GO(s,:);CRed = CO(s,:);
            else % do not reduce
                GRed = GO;CRed = CO;
            end
        end
        function [GAdd,CAdd] = regUpSample(GO,CO,info)
            s = randsample(info.SampleInd,info.nrSyn,true);% sampling starting faces
            rK = randsample(info.K,info.nrSyn,true);% sampling neighbor ID 1 - K
            nb = info.SampleIndB(sub2ind(size(info.SampleIndB),rK(:),s(:)));% sampling neighbors
            alpha = rand(1,info.nrSyn);% sampling interpolation factor
            GAdd = GO(s,:) + (GO(nb,:)-GO(s,:)).*repmat(alpha(:),1,size(GO,2));% creating new depvar
            CAdd = CO(s,:) + (CO(nb,:)-CO(s,:)).*repmat(alpha(:),1,size(CO,2));% creating new indvar
        end
        function F = DOB_SCV_DM(K,DepVar,G,type)
            D = squareform(pdist(DepVar,type));
            n = size(D,1);
            if nargin < 3, G = ones(n,1); end% simply distribute similar data features accross
            if isempty(G), G = ones(n,1); end% one big group
            % dealing with unknown group memberships
            indnan = find(isnan(G)); %#ok<*EFIND>
            if isempty(indnan)
                F = nan*zeros(1,n);
                Fid = (2:K);
                lC = unique(G);
                nC = length(lC);
                Gfor = G;
                for i=1:nC
                    indC = find(Gfor==lC(i));
                    nrC = length(indC);
                    while nrC>0
                        e = randsample(1:nrC,1);
                        F(indC(e)) = 1;% assign to first fold
                        Gfor(indC(e)) = nan;% reducing the group members to deplete them and stop the while
                        if nrC==1,break;end
                        re = setdiff(1:nrC,e);
                        distances = D(indC(e),indC(re));
                        %distances = sqrt(sum((repmat(D(indC(e),:),length(re),1)-D(indC(re),:)).^2,2));
                        [~,index] = sort(distances,'ascend');
                        if length(index)<K-1
                            F(indC(re(index))) = Fid(1:length(index));
                            Gfor(indC(re(index))) = nan;
                        else
                            % assign to folds
                            F(indC(re(index(1:K-1)))) = Fid;
                            Gfor(indC(re(index(1:K-1)))) = nan;
                        end
                        indC = find(Gfor==lC(i));
                        nrC = length(indC);
                    end
                end
            else
                indrest = setdiff((1:n),indnan);
                Fnan = DOB_SCV_DM(K,D(indnan,indnan));% distribute nan values according to similarity in D only
                Frest = DOB_SCV_DM(K,D(indrest,indrest),G(indrest));% distribute as originally planned
                F(indnan) = Fnan;
                F(indrest) = Frest;
            end
        end
    end
    methods (Static = true) % FEATURE SELECTION
        function [out,COR] = fsCOR(FS,X,maxFeature)
            if isempty(FS), out = []; COR = [];return; end
            nFS = size(FS,2);
            COR = nan*zeros(1,nFS);
            parfor f=1:nFS
                tmp = corrcoef(X,FS(:,f)); %#ok<*PFBNS>
                COR(f) = tmp(1,2);
            end
            COR = abs(COR);
            [~,sortind] = sort(COR,'descend');
            if nFS>=maxFeature
                out = sortind(1:maxFeature);
            else
                out = sortind;
            end
            COR = COR(out);
        end
        function [out,EER] = fsROC(FS,X,maxFeature)
            if isempty(FS), out = []; EER = []; return; end
            nFS = size(FS,2);
            EER = nan*zeros(1,nFS);
            Tind = find(X==1);
            Find = find(X==-1);
            parfor f=1:nFS
                tmp = PPMMV.getEER(FS(Tind,f),FS(Find,f));
                if tmp>0.5, tmp = PPMMOD.getEER(-1*FS(Tind,f),-1*FS(Find,f));end
                EER(f) = tmp;
            end
            [~,sortind] = sort(EER,'ascend');
            if nFS>=maxFeature
                out = sortind(1:maxFeature);
            else
                out = sortind;
            end
            EER = EER(out);
            EER = 0.5-EER;% the higher the better for optimization reasons later
        end
        function out = fsMRMTR(H,maxFeature)
            %[~,~,H] = PPMMV.prepMIfs(FS,X,qs);
            %if maxFeature>size(FS,2), maxFeature = size(FS,2); end
            %select first feature as the one with max MI with C
            max_MI=0;firstFeature=1;
            dim = size(H,1);
            for i=1:dim
                CMI=H(i,i);
                if CMI>max_MI
                    max_MI=CMI;
                    firstFeature=i;
                end
            end
            for i=1:dim %create the JMI matrix
                for j=1:dim
                    if j==i, continue;end;
                    H(i,j)=H(i,i)+H(i,j);  %I(Xi;C) + I(Xj;C|Xi)= I(Xi,Xj;C)
                end
            end
            best_fs=zeros(1,maxFeature);
            best_fs(1)=firstFeature;
            selected=zeros(1,dim);
            selected(best_fs(1))=1;
            for j=2:maxFeature
                max_inc=-inf;
                bestFeature=0;
                for i=1:dim
                    if selected(i), continue;end;
                    totalJMI=sum(H(i,best_fs(1:j-1)));
                    if totalJMI>max_inc
                        max_inc=totalJMI;
                        bestFeature=i;
                    end
                end
                best_fs(j)=bestFeature;
                selected(bestFeature)=1;
            end
            out = best_fs;
        end
        function [out,weights] = fsSPECCMI(H,maxFeature)
            %[~,~,H] = PPMMV.prepMIfs(FS,X,qs);
            %if maxFeature>size(FS,2), maxFeature = size(FS,2); end
            dim = size(H,1);
            H=(H+H')/2;
            [V,~] = eigs(H,1);
            x=V(:,1);
            x=abs(x/norm(x));
            y=zeros(dim,2);
            y(:,1)=-x;
            for i=1:dim, y(i,2)=i;end;
            y=sortrows(y);
            y=[y(:,2) -y(:,1)];
            SPECCMI_Fs=y(:,1);
            weights=y(:,2);
            weights = weights(1:maxFeature);
            %weights = weights-min(weights);
            %weights = weights/max(weights);
            out = SPECCMI_Fs(1:maxFeature);
        end
        function out = fsCMIM(H,maxFeature)
            %[~,~,H] = PPMMV.prepMIfs(FS,X,qs);
            %if maxFeature>size(FS,2), maxFeature = size(FS,2); end
            dim = size(H,1);
            H=H'; %so that H(i,j)=I(X_i;C|Xj)??
            %select first feature as the one with max MI with C
            max_MI=0;firstFeature=1;
            for i=1:dim
                CMI=H(i,i);
                if CMI>max_MI
                    max_MI=CMI;
                    firstFeature=i;
                end
            end
            best_fs=zeros(1,maxFeature);
            best_fs(1)=firstFeature;
            %fprintf('Adding first feature %d\n',firstFeature);
            selected=zeros(1,dim);
            selected(best_fs(1))=1;
            for j=2:maxFeature
                max_red=-inf; %max of min conditional relevancy
                bestFeature=0;
                for i=1:dim
                    if selected(i) continue;end;
                    [mini]=min(H(i,best_fs(1:j-1)));
                    if max_red<mini
                        max_red=mini;
                        bestFeature=i;
                    end
                end
                
                best_fs(j)=bestFeature;
                selected(bestFeature)=1;
                %     fprintf('Adding %d\n',bestFeature);
            end
            out = best_fs;
        end
        function out = fsCIFE(H,maxFeature)
            %[~,~,H] = PPMMV.prepMIfs(FS,X,qs);
            %if maxFeature>size(FS,2), maxFeature = size(FS,2); end
            dim = size(H,1);
            H=H'; %so that H(i,j)=I(X_i;C|Xj)??
            %select first feature as the one with max MI with C
            max_MI=0;firstFeature=1;
            for i=1:dim
                CMI=H(i,i);
                if CMI>max_MI
                    max_MI=CMI;
                    firstFeature=i;
                end
            end
            best_fs=zeros(1,maxFeature);
            best_fs(1)=firstFeature;
            %fprintf('Adding first feature %d\n',firstFeature);
            selected=zeros(1,dim);
            selected(best_fs(1))=1;
            for j=2:maxFeature
                bestobj=-inf; %minimum conditional redundancy
                bestFeature=0;
                for i=1:dim
                    if selected(i) continue;end;
                    obj=(2-j)*H(i,i)+ sum(H(i,best_fs(1:j-1))); %CIFE criterion, as in PR paper
                    if bestobj<obj
                        bestobj=obj;
                        bestFeature=i;
                    end
                end
                best_fs(j)=bestFeature;
                selected(bestFeature)=1;
                %     fprintf('Adding %d\n',bestFeature);
            end
            out = best_fs;
        end
        function [a,C,H] = prepMIfs(FS,X,qs)
            index = find(~isnan(X));
            C = X(index);
            C(C==-1) = 2;
            data = normalize_max(FS(index,:));
            a=myQuantileDiscretize(data,qs);
            H=computeCMImatrix_4([a C]);
        end
        function out = standarizeFS(fs)
            n = size(fs,1);
            avgs = mean(fs,1);
            out = fs-repmat(avgs,n,1);
            stds = nanstd(out,1,1);
            out = fs./repmat(stds,n,1);
        end
    end
    methods (Static = true) % PHENOMIC PREDICTOR
        function out = matchthresholding(matches,T)
            out = matches>=T;
        end
        function out = class2score(classes)
            out = 1-classes;
        end
        function out = normVAL(val)
            out = val;
            for i=1:1:size(val,2)
                m = min(out(:,i));
                out(:,i) = out(:,i)-m;
                M = max(out(:,i));
                out(:,i) = out(:,i)/M;
            end
        end
    end
    methods (Static = true) % BIOMETRICS
        function out = pRAWeight(A,Dim,Table)
            min = -log(BASEMOD.lookUppRA(0,Dim,Table));
            out = -log(BASEMOD.lookUppRA(A,Dim,Table));
            out = max(out-min,0);
        end
        function out = pAWeight(pA,minpA)
            min = -log(minpA);
            out = -log(pA);
            out = max(out-min,0);
        end
        function [EER,G,AUC,x,y,TH,Y] = getEER(tmatches,fmatches)
            g = [ones(1,length(tmatches)), -1*ones(1,length(fmatches))];
            [sorted,order] = sort([tmatches;fmatches],'ascend');
            g = g(order);
            true_neg = g == -1;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
            true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
            y = tpf(1:end-1)+dy./2;
            x = fpf(1:end-1)+dx./2;
            [Y,indy] = max(y-x);
            yn = 1-y;d = abs(x-yn);
            [~,ind] = min(d);
            EER = ((x(ind)+yn(ind))/2);
            AUC = sum(dx.*y);
            FN = tpf(ind+1)*na;TN = fpf(ind+1)*nn;
            TP = na-FN;FP = nn-TN;
            G = 1-sqrt((TP/(TP+FN))*(TN/(TN+FP)));
            T(1,1) = TP;
            T(1,2) = FP;
            T(2,1) = FN;
            T(2,2) = TN;
            %disp(num2str(T));
            %disp(num2str(1-sorted(ind)));
            TH(1) = 1-sorted(ind);
            TH(2) = 1-sorted(indy);
            %disp(['EER TH: ' num2str(TH(1))]);
            %disp(['Y TH: ' num2str(TH(2))]);
        end
        function [out,R] = getRANK(tmatches,fmatches)
            R = sum(repmat(tmatches,1,length(fmatches))>repmat(fmatches',length(tmatches),1),2);
            R = R./(length(fmatches)+1);
            RM = mean(R);
            R1 = sum(R<=0.01)/length(R);
            R10 = sum(R<=0.10)/length(R);
            R20 = sum(R<=0.20)/length(R);
            out = [R1 R10 R20 RM];
        end
        function [R1,R10,R20,CR] = getRANKSTAT(tmatches,fmatches)
            %tmatches = TMatches;
            %fmatches = FMatches;
            [nT,nF] = size(fmatches);
            R = (sum(fmatches<repmat(tmatches,1,nF),2)+1);%/(nF+1);
            R = (R/(nF+1)).*100;
            CR = zeros(1,100);
            for cr=1:1:100
                CR(1,cr) = (sum(R<=cr)/nT).*100;
            end
            R1 = CR(1);
            R10 = CR(10);
            R20 = CR(20);
        end
    end
    methods (Static = true) % INTERFACING
        function out = notNAN(in)
            index = (1:size(in,1));
            [i,~] = find(isnan(in));
            out = setdiff(index,unique(i));
        end
        function out = getL1CL1DepVar(in)
            if ~iscell(in), out = in;return;end
            out = in{1}.DepVar{1};
        end
        function [PREC,REC,G,TNF,XH,YH,ACC] = getHardClass(Tclass,Fclass)
            Pc = length(Tclass);Nc = length(Fclass);
            TP = length(find(Tclass==1));
            TN = length(find(Fclass==0));TNF = TN/Nc;
            FN = Pc-TP;FP = Nc-TN;
            PREC = TP/(TP+FP);
            REC = TP/(TP+FN);
            G = sqrt((TP/(TP+FN))*(TN/(TN+FP)));
            YH = TP/Pc;
            XH = 1-TN/Nc;
            ACC = (TP+TN)/(Pc+Nc);
        end
        function out = convertUInt(in)
            if isempty(in), out = []; return; end
            M = max(in);
            if M<=intmax('uint8'),out = uint8(in);return;end
            if M<=intmax('uint16'), out = uint16(in);return;end
            if M<=intmax('uint32'), out = uint32(in);return;end
            out = uint64(in);
        end
        function out = getForcedPermutations(t,nS)           % Created by Dorothy
            out = zeros(t,nS);
            for p = 1:1:t,
                out(p,:) = randperm(nS);
            end
        end
        function out = getParent(level,cluster)              % Created by Dorothy
            if level == 1,  out = []; return; end
            if cluster > 2^(level-1), out = []; return; end
            out = ceil(cluster/2);
        end
        function out = getChildren(obj,level,cluster)        % Created by Dorothy
            if level >= obj.nLev; out = []; return; end
            if cluster > 2^(level-1), out = []; return; end
            out = [(2*cluster)-1 2*cluster];
        end
        function out = getParentEvolution(obj,level,cluster) % Created by Dorothy
            if level > obj.nLev; out = []; return; end
            if cluster > 2^(level-1), out = []; return; end
            out = zeros(1,level);
            out(level) = cluster;
            for l = level-1:-1:1,
                out(l) = getParent(l+1,out(l+1));
            end
        end
        function out = getClind(obj,level,cluster)           % Created by Dorothy
            if level > obj.nLev; out = []; return; end
            if cluster > 2^(level-1), out = []; return; end
            out = sum(2.^([1:1:level-1]-1));
            out = out + cluster;
        end
        function [p, pPerms] = nonParametricCombination (p2combine, pPerms2combine)     % Created by Dorothy
            t = size(p2combine,1);
            % Calculate the combined statistics: "fisher"
            stat = -1* sum(log(p2combine),2);
            statPerm = -1* sum(log(pPerms2combine),2);
            % Calculate the p_values of the true data and the permutated data
            p = sum([statPerm;stat] >= repmat(stat,t+1,1),1)/(t+1);
            if nargout >= 2,
                pPerms = sum(repmat([statPerm;stat]',t,1) >= repmat(statPerm,1,t+1),2)./(t+1);
            end
        end
    end
end