classdef PPMMV2 <superClassLight
% PROPERTIES    
    properties % GENERAL INTERFACING
       nLev = 6;
       PosClassBal = [];
       nTraining = 0;
    end
    properties (Dependent = true)
       nLC;
       nS;
       AvgX;
       StdX;
       PercNan;
       PosClassBalTot;
    end
    properties (Hidden = true)
       Levels;
       Clusters;
    end
    properties % ASSOCIATION TESTING
       T;
       pT;
       AvgT;
       pFastT;
       R2;
       pR2;
       ppR2;
       CCA;
       pCCA;
       npCCA;
       nppCCA;
       bCCA;
       bPLSR;
       A;
       pA;
       ppA;
       pAR;
       ppAR;
       pAP;
       ppAP;
       FA;
       pFA;
       EER;
    end
    properties (Dependent = true)
       MaxLogpT;
       LogpFastT
    end
    properties (Hidden = true)
       TpA = single(0.05);
       TpAR = single(0.05);
       TpFA = single(0.05);
       pCrit = [];
       HitType = 'pAR';
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
       HMatrixA;              % created by Dorothy
       HMatrixR2;              % created by Dorothy
       HMatrixpR2;              % created by Dorothy
       HMatrixpA;               % created by Dorothy
       HMatrixpAR;              % created by Dorothy
       HMatrixppR2;             % created by Dorothy
       HMatrixppA;              % created by Dorothy
       HMatrixppAR;             % created by Dorothy
       HMatrixppAP;             % created by Dorothy 2
       HMatrixpAP;             % created by Dorothy 2
       HMatrixT;
       HMatrixFastpT;
       HMatrixMaxpT;
       HMatrixCCA;
       HMatrixpCCA;
       %HMatrixnpCCA;
       %HMatrixnppCCA;
       nHitspA;
       nHitspAR;
       nHitsppA;
       nHitsppAR;
       nHitsppAP;
       nHitspAP;
       nHitspR2;
       nHitsppR2;
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
       function obj = PPMMV2(varargin)
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
           %if isempty(obj.A), out = []; return;end
           %out = (obj.pA<=obj.TpA).*(obj.pAR<=obj.TpAR);
           if strcmp(obj.HitType,'none'), out = ones(obj.nLC,1); return; end
           eval(['tmp = obj.' obj.HitType ';']);
           out = tmp<=obj.pCrit;
        end
       function out = get.MatrixHits(obj)
            out = List2Matrix(obj,obj.Hits);
        end
       function out = get.nHits(obj)
            %if isempty(obj.A), out = []; return; end
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
       function out = get.HMatrixR2(obj)
            bk = obj.pCrit;
            obj.pCrit = [];
            out = List2HierachicalMatrix(obj,obj.R2);
            obj.pCrit = bk;
       end
       function out = get.HMatrixA(obj)
            bk = obj.pCrit;
            obj.pCrit = [];
            out = List2HierachicalMatrix(obj,obj.A);
            obj.pCrit = bk;
       end
       function out = get.HMatrixpR2(obj)       % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.pR2);
       end
       function out = get.HMatrixpA(obj)        % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.pA);
       end
       function out = get.HMatrixpAP(obj)        % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.pAP);
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
       function out = get.HMatrixT(obj)      % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.AvgT);
       end
       function out = get.HMatrixFastpT(obj)      % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.LogpFastT);
       end
       function out = get.HMatrixMaxpT(obj)      % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.MaxLogpT);
       end
       function out = get.HMatrixCCA(obj)      % created by Dorothy
            out = List2HierachicalMatrix(obj,obj.CCA);
       end
       function out = get.HMatrixpCCA(obj)      % created by Dorothy
            out = List2HierachicalMatrix(obj,-log10(obj.pCCA));
       end
       function out = get.nHitspA(obj)
                out = sum(obj.pA<=obj.pCrit);
       end
       function out = get.nHitspAR(obj)
                out = sum(obj.pAR<=obj.pCrit);
       end
       function out = get.nHitsppAR(obj)
                out = sum(obj.ppAR<=obj.pCrit);
       end
       function out = get.nHitsppAP(obj)
                out = sum(obj.ppAP<=obj.pCrit);
       end
       function out = get.nHitspAP(obj)
                out = sum(obj.pAP<=obj.pCrit);
       end
       function out = get.nHitsppA(obj)
                out = sum(obj.ppA<=obj.pCrit);
       end
       function out = get.nHitspR2(obj)
                out = sum(obj.pR2<=obj.pCrit);
       end
       function out = get.nHitsppR2(obj)
                out = sum(obj.ppR2<=obj.pCrit);
       end
       function out = get.PosClassBal(obj)
           if ~isempty(obj.PosClassBal), out = obj.PosClassBal; return; end
           if isempty(obj.X), out = 0; return; end
           tmp = obj.XREG;
           tmp(isnan(tmp)) = [];
           Tind = tmp==obj.PosClass;
           out = sum(Tind)/length(tmp);
       end
       function out = get.PosClassBalTot(obj)
           if obj.nTraining==0, out = 0; return; end
           perc = obj.nTraining/obj.nS;
           out = obj.PosClassBal.*perc;
       end
       function out = get.MaxLogpT(obj)
          out = -log10(obj.pT);
          out = max(out,[],2);
       end
       function out = get.LogpFastT(obj)
                out = -log10(obj.pFastT); 
       end
    end
    methods % GENERAL SETTING
       function obj = set.nLev(obj,in)
           obj.nLev = single(in);
           storeLevels(obj);
           storeClusters(obj);
       end
       function obj = set.FSInd(obj,in)
            obj.FSInd = PPMMV2.convertUInt(in);
       end
       function obj = set.maxFS(obj,in)
          obj.maxFS = PPMMV2.convertUInt(in); 
       end
       function obj = set.nrTrees(obj,in)
          obj.nrTrees = PPMMV2.convertUInt(in); 
       end
    end
    methods % ASSOCIATION TESTING
       function out = runR2Test(obj,DepVar,COV,t)
           R2 = nan*zeros(1,obj.nLC);
           pR2 = nan*zeros(1,obj.nLC);
           parfor i=1:obj.nLC
              [R2(i),pR2(i)] =  PPMMV2.testPartialPLSR(obj.XREG,getDepVar(obj,DepVar,i),getCOV(obj,COV,i),t); %#ok<*PFBNS>
           end
           out.R2 = R2;out.pR2 = pR2;
           obj.R2 = R2;obj.pR2 = pR2;
       end
       function out = runParReg(obj,DepVar,COV)
           tmp = PPMMV2.getL1CL1DepVar(DepVar);
           nVar = size(tmp,2);
           T = nan*zeros(obj.nLC,nVar);
           pT = nan*zeros(obj.nLC,nVar);
           AvgT = nan*zeros(1,obj.nLC);
           pFastT = nan*zeros(1,obj.nLC);
           parfor i=1:obj.nLC
              [T(i,:),pT(i,:),pFastT(i),AvgT(i)] =  PPMMV2.parReg(obj.XREG,getDepVar(obj,DepVar,i),getCOV(obj,COV,i)); %#ok<*PFBNS>
           end
           out.T = T;out.pT = pT;out.AvgT = AvgT;out.pFastT = pFastT;
           obj.T = T;obj.pT = pT;obj.AvgT = AvgT;obj.pFastT = pFastT;
       end
       function out = runParPLSReg(obj,DepVar,COV)
           tmp = PPMMV2.getL1CL1DepVar(DepVar);
           nVar = size(tmp,2);
           T = nan*zeros(obj.nLC,nVar);
           pT = nan*zeros(obj.nLC,nVar);
           AvgT = nan*zeros(1,obj.nLC);
           pFastT = nan*zeros(1,obj.nLC);
           B = cell(1,obj.nLC);
           parfor i=1:obj.nLC
              [T(i,:),pT(i,:),pFastT(i),AvgT(i),B{i}] =  PPMMV2.parPLSReg(obj.XREG,getDepVar(obj,DepVar,i),getCOV(obj,COV,i)); %#ok<*PFBNS>
           end
           out.T = T;out.pT = pT;out.AvgT = AvgT;out.pFastT = pFastT;
           obj.T = T;obj.pT = pT;obj.AvgT = AvgT;obj.pFastT = pFastT;obj.bPLSR = B;
       end
       function out = runParCCA(obj,DepVar,COV)
           F = nan*zeros(1,obj.nLC);
           pF = nan*zeros(1,obj.nLC);
           B = cell(1,obj.nLC);
           parfor i=1:obj.nLC
              [F(i),pF(i),B{i}] =  PPMMV2.parCCA(obj.XREG,getDepVar(obj,DepVar,i),getCOV(obj,COV,i)); %#ok<*PFBNS>
           end
           out.F = F;out.pF = pF;out.B = B;
           obj.CCA = F;obj.pCCA = pF;obj.bCCA = B;
       end
       function out = runNonParCCA(obj,DepVar,COV,t)
           F = nan*zeros(1,obj.nLC,'single');
           pF = nan*zeros(1,obj.nLC,'single');
           ppF = nan*zeros(1,obj.nLC,'single');
           ppFPerms = nan*zeros(t,obj.nLC,'single');
           L1CL1DepVar =  PPMMV2.getL1CL1DepVar(DepVar);
           nS = size(L1CL1DepVar,1);
           forcedPerms = PPMMV2.getForcedPermutations(t,nS);        
           % level with no children
           l = obj.nLev;
           clind = zeros(1,2^(l-1));
           tmpF = zeros(1,2^(l-1));
           %tmppF = zeros(1,2^(l-1));
           tmpppF = zeros(1,2^(l-1),'single');
           tmpppFPerms = zeros(t,2^(l-1),'single');
           parfor cl = 1:1:2^(l-1),
               tmpclind = getClind(obj,l,cl);
               [tmpppF(1,cl), tmpppFPerms(:,cl),tmpF(cl)] = PPMMV2.nonparCCA(obj.XREG,getDepVar(obj,DepVar,tmpclind),getCOV(obj,COV,tmpclind),t,forcedPerms);
               clind(cl) = tmpclind;
           end
           F(1,clind) = tmpF;
           pF(1,clind) = tmpppF;
           ppF(1,clind) = tmpppF;
           ppFPerms(:,clind) = tmpppFPerms; 
           % levels with children
           for l = obj.nLev-1:-1:1,
               %l=5
               clind = zeros(1,2^(l-1),'single');
               tmpF = zeros(1,2^(l-1),'single');
               tmppF = zeros(1,2^(l-1),'single');
                tmpppF = zeros(1,2^(l-1),'single');
                tmpppFPerms = zeros(t,2^(l-1),'single');
                parfor cl = 1:1:2^(l-1),
                    %cl = 1
                    % get the p values for the parent cluster
                    tmpclind = getClind(obj,l,cl);
                    [ppF_temp, ppFPerms_temp,tmpF(cl)] = PPMMV2.nonparCCA(obj.XREG,getDepVar(obj,DepVar,tmpclind),getCOV(obj,COV,tmpclind),t,forcedPerms);
                    clind(cl) = tmpclind;
                    % get the p values for the two children clusters
                    children = getChildren(obj,l,cl);
                    clind1 = getClind(obj,l+1,children(1));
                    clind2 = getClind(obj,l+1,children(2));
                    tmppF(1,cl) = ppF_temp;
                    % nonparametric combination of the p values of the parent and the two children
                    [tmpppF(1,cl), tmpppFPerms(:,cl)] = PPMMV2.nonParametricCombination ([ppF_temp,ppF(clind1),ppF(clind2)], [ppFPerms_temp,ppFPerms(:,clind1),ppFPerms(:,clind2)]); 
                end
                F(1,clind) = tmpF;
                pF(1,clind) = tmppF;
                ppF(1,clind) = tmpppF;
                ppFPerms(:,clind) = tmpppFPerms;
            end
            out.F = F;out.pF = pF;out.ppF = ppF;
            obj.CCA = F;obj.npCCA = pF;obj.nppCCA = ppF;
       end
       function out = runNonParCCAv2(obj,DepVar,COV,t)
           L1CL1DepVar =  PPMMV2.getL1CL1DepVar(DepVar);
           nS = size(L1CL1DepVar,1);
           if nS < intmax('uint8')
               str = 'uint8';
           elseif nS < intmax('uint16')
               str = 'uint16';
           elseif nS < intmax('uint32')
               str = 'uint32';
           else
               str = 'uint64';
           end
           %forcedPerms = PPMMV2.getForcedPermutations(t,nS);
           %str = class(forcedPerms);
           F = nan*zeros(1,obj.nLC,'single');
           pF = nan*zeros(1,obj.nLC,str);
           ppF = nan*zeros(1,obj.nLC,str);
           %ppFPerms = nan*zeros(t,obj.nLC,str);  
           % level with no children
           l = obj.nLev;
           %clind = zeros(1,2^(l-1),'uint16');
           %tmppF = zeros(1,2^(l-1));
           %tmpppF = zeros(1,2^(l-1),str);
           ChildrenFPerms = zeros(t,2^(l-1),str);
           for cl = 1:1:2^(l-1)
               disp(['Level: ' num2str(obj.nLev) ' Cluster: ' num2str(cl)]);
               clind = uint16(getClind(obj,l,cl));
               %disp(num2str(tmpclind));
               if cl<3, tic; end
               [pF(clind), ChildrenFPerms(:,cl),F(clind)] = PPMMV2.nonparCCAv2(obj.XREG,getDepVar(obj,DepVar,clind),getCOV(obj,COV,clind),t,str);
               if cl<3, toc; end
               ppF(clind) = pF(clind);
           end
           % levels with children
           for l = obj.nLev-1:-1:1,
               %l=5
               % l=4;
                %clind = zeros(1,2^(l-1),'uint16');
                %tmpF = zeros(1,2^(l-1),'single');
                %tmppF = zeros(1,2^(l-1),str);
                %tmpppF = zeros(1,2^(l-1),str);
                ParentFPerms = zeros(t,2^(l-1),str);
                for cl = 1:1:2^(l-1)
                    disp(['Level: ' num2str(l) ' Cluster: ' num2str(cl)]);
                    %cl = 1
                    % get the p values for the parent cluster
                    clind = uint16(getClind(obj,l,cl));
                    [pF(clind), ppFPerms_temp,F(clind)] = PPMMV2.nonparCCAv2(obj.XREG,getDepVar(obj,DepVar,clind),getCOV(obj,COV,clind),t,str);
                    %clind(cl) = tmpclind;
                    % get the p values for the two children clusters
                    children = getChildren(obj,l,cl);
                    clind1 = uint16(getClind(obj,l+1,children(1)));
                    clind2 = uint16(getClind(obj,l+1,children(2)));
                    % nonparametric combination of the p values of the parent and the two children
                    p2combine = [pF(clind),ppF(clind1),ppF(clind2)];
                    pPerms2combine = zeros(t,3,str);
                    pPerms2combine(:,1) = ppFPerms_temp;clear ppFPerms_temp;
                    pPerms2combine(:,2) = ChildrenFPerms(:,children(1));
                    pPerms2combine(:,3) = ChildrenFPerms(:,children(2));
                    [ppF(clind), ParentFPerms(:,cl)] = PPMMV2.nonParametricCombinationv2(p2combine,pPerms2combine,str); 
                end
                ChildrenFPerms = ParentFPerms;clear ParentFPerms;
            end
            out.F = F;out.pF = double(pF)/(t+1);out.ppF = double(ppF)/(t+1);
            obj.CCA = F;obj.npCCA = double(pF)/(t+1);obj.nppCCA = double(ppF)/(t+1);
       end
       function out = runAngleTest(obj,DepVar,COV,maxM,AngleTable,t)
           out = cell(1,obj.nLC);
           parfor i=1:obj.nLC
              [out{i}.A,out{i}.M,out{i}.STAT] = PPMMV2.AngleTest([getCOV(obj,COV,i) obj.XREG],getDepVar(obj,DepVar,i),AngleTable,t,maxM,0); 
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
            R2 = nan*zeros(1,obj.nLC);
            pR2 = nan*zeros(1,obj.nLC);
            ppR2 = nan*zeros(1,obj.nLC);
            ppR2Perms = nan*zeros(t,obj.nLC);
            L1CL1DepVar =  PPMMV2.getL1CL1DepVar(DepVar);
            nS = size(L1CL1DepVar,1);
            forcedPerms = PPMMV2.getForcedPermutations(t,nS);        
            % level with no children
            l = obj.nLev;
            clind = zeros(1,2^(l-1));
            tmpR2 = zeros(1,2^(l-1));
            %tmppR2 = zeros(1,2^(l-1));
            tmpppR2 = zeros(1,2^(l-1));
            tmpppR2Perms = zeros(t,2^(l-1));
            parfor cl = 1:1:2^(l-1),
                tmpclind = getClind(obj,l,cl);
                [tmpppR2(1,cl), tmpppR2Perms(:,cl),tmpR2(cl)] = PPMMV2.testPartialPLSR_Permutated(obj.XREG,getDepVar(obj,DepVar,tmpclind),getCOV(obj,COV,tmpclind),t,forcedPerms);
                clind(cl) = tmpclind;
            end
            R2(1,clind) = tmpR2;
            pR2(1,clind) = tmpppR2;
            ppR2(1,clind) = tmpppR2;
            ppR2Perms(:,clind) = tmpppR2Perms; 
            % levels with children
            for l = obj.nLev-1:-1:1,
                %l=5
                clind = zeros(1,2^(l-1));
                tmpR2 = zeros(1,2^(l-1));
                tmppR2 = zeros(1,2^(l-1));
                tmpppR2 = zeros(1,2^(l-1));
                tmpppR2Perms = zeros(t,2^(l-1));
                parfor cl = 1:1:2^(l-1),
                    %cl = 1
                    % get the p values for the parent cluster
                    tmpclind = getClind(obj,l,cl);
                    [ppR2_temp, ppR2Perms_temp,tmpR2(cl)] = PPMMV2.testPartialPLSR_Permutated(obj.XREG,getDepVar(obj,DepVar,tmpclind),getCOV(obj,COV,tmpclind),t,forcedPerms);
                    clind(cl) = tmpclind;
                    % get the p values for the two children clusters
                    children = getChildren(obj,l,cl);
                    clind1 = getClind(obj,l+1,children(1));
                    clind2 = getClind(obj,l+1,children(2));
                    tmppR2(1,cl) = ppR2_temp;
                    % nonparametric combination of the p values of the parent and the two children
                    [tmpppR2(1,cl), tmpppR2Perms(:,cl)] = PPMMV2.nonParametricCombination ([ppR2_temp,ppR2(clind1),ppR2(clind2)], [ppR2Perms_temp,ppR2Perms(:,clind1),ppR2Perms(:,clind2)]); 
                end
                R2(1,clind) = tmpR2;
                pR2(1,clind) = tmppR2;
                ppR2(1,clind) = tmpppR2;
                ppR2Perms(:,clind) = tmpppR2Perms;
            end
            out.R2 = R2;out.pR2 = pR2;out.ppR2 = ppR2;
            obj.R2 = R2;obj.pR2 = pR2;obj.ppR2 = ppR2;
        end
       function out = runPropagatingAngleTest(obj,DepVar,COV,maxM,AngleTable,t,t2)    % created by Dorothy (long computation time)
            A = nan*zeros(1,obj.nLC);
            pA = nan*zeros(1,obj.nLC);
            ppA = nan*zeros(1,obj.nLC);
            ppAPerms = nan*zeros(t2,obj.nLC);
            pAR = nan*zeros(1,obj.nLC);
            ppAR = nan*zeros(1,obj.nLC);
            ppARPerms = nan*zeros(t2,obj.nLC);
            pAP = nan*zeros(1,obj.nLC);
            ppAP = nan*zeros(1,obj.nLC);
            ppAPPerms = nan*zeros(t2,obj.nLC);
            L1CL1DepVar =  PPMMV2.getL1CL1DepVar(DepVar);
            nS = size(L1CL1DepVar,1);
            forcedPerms = PPMMV2.getForcedPermutations(t2,nS);
            % level with no children
            l = obj.nLev;
            clind = zeros(1,2^(l-1));
            %tmppA = zeros(1,2^(l-1));
            tmpppA = zeros(1,2^(l-1));
            tmpppAPerms = zeros(t2,2^(l-1));
            %tmppAR = zeros(1,2^(l-1));
            tmpppAR = zeros(1,2^(l-1));
            tmpppARPerms = zeros(t2,2^(l-1));
            %tmppAP = zeros(1,2^(l-1));
            tmpppAP = zeros(1,2^(l-1));
            tmpppAPPerms = zeros(t2,2^(l-1));
            tmpA = zeros(1,2^(l-1));
            parfor cl = 1:1:2^(l-1),
                tmpclind = getClind(obj,l,cl);
                [tmpppA(cl), tmpppAPerms(:,cl), tmpppAR(cl), tmpppARPerms(:,cl), tmpppAP(cl), tmpppAPPerms(:,cl),tmpA(cl)] = PPMMV2.AngleTest_Permutated([getCOV(obj,COV,tmpclind) obj.XREG],getDepVar(obj,DepVar,tmpclind),AngleTable,t,maxM,0,t2,forcedPerms);
                clind(cl) = tmpclind;
            end
            A(clind) = tmpA;
            pA(1,clind) = tmpppA;
            ppA(1,clind) = tmpppA;
            ppAPerms(:,clind) = tmpppAPerms;
            pAR(1,clind) = tmpppAR;
            ppAR(1,clind) = tmpppAR;
            ppARPerms(:,clind) = tmpppARPerms;
            pAP(1,clind) = tmpppAP;
            ppAP(1,clind) = tmpppAP;
            ppAPPerms(:,clind) = tmpppAPPerms;
            % levels with children
            for l = obj.nLev-1:-1:1,
                clind = zeros(1,2^(l-1));
                tmpA = zeros(1,2^(l-1));
                tmppA = zeros(1,2^(l-1));
                tmpppA = zeros(1,2^(l-1));
                tmpppAPerms = zeros(t2,2^(l-1));
                tmppAR = zeros(1,2^(l-1));
                tmpppAR = zeros(1,2^(l-1));
                tmpppARPerms = zeros(t2,2^(l-1));
                tmppAP = zeros(1,2^(l-1));
                tmpppAP = zeros(1,2^(l-1));
                tmpppAPPerms = zeros(t2,2^(l-1));
                parfor cl = 1:1:2^(l-1),
                    % get the p values for the parent cluster
                    tmpclind = getClind(obj,l,cl);
                    [ppA_temp, ppAPerms_temp, ppAR_temp, ppARPerms_temp, ppAP_temp, ppAPPerms_temp, tmpA(cl)] = PPMMV2.AngleTest_Permutated([getCOV(obj,COV,tmpclind) obj.XREG],getDepVar(obj,DepVar,tmpclind),AngleTable,t,maxM,0,t2,forcedPerms);
                    clind(cl) = tmpclind;
                    % get the p values for the two children clusters
                    children = getChildren(obj,l,cl);
                    clind1 = getClind(obj,l+1,children(1));
                    clind2 = getClind(obj,l+1,children(2));  
                    tmppA(cl) = ppA_temp;
                    tmppAR(cl) = ppAR_temp;
                    tmppAP(cl) = ppAP_temp;
                    % nonparametric combination of the p values of the parent and the two children
                    [tmpppA(cl), tmpppAPerms(:,cl)] = PPMMV2.nonParametricCombination ([ppA_temp,ppA(clind1),ppA(clind2)], [ppAPerms_temp,ppAPerms(:,clind1),ppAPerms(:,clind2)]); 
                    [tmpppAR(cl), tmpppARPerms(:,cl)] = PPMMV2.nonParametricCombination ([ppAR_temp,ppAR(clind1),ppAR(clind2)], [ppARPerms_temp,ppARPerms(:,clind1),ppARPerms(:,clind2)]); 
                    [tmpppAP(cl), tmpppAPPerms(:,cl)] = PPMMV2.nonParametricCombination ([ppAP_temp,ppAP(clind1),ppAP(clind2)], [ppAPPerms_temp,ppAPPerms(:,clind1),ppAPPerms(:,clind2)]); 
                end
                A(clind) = tmpA;
                pA(1,clind) = tmppA;
                ppA(1,clind) = tmpppA; %#ok<*PROPLC>
                ppAPerms(:,clind) = tmpppAPerms;
                pAR(1,clind) = tmppAR;
                ppAR(1,clind) = tmpppAR;
                ppARPerms(:,clind) = tmpppARPerms;
                pAP(1,clind) = tmppAP;
                ppAP(1,clind) = tmpppAP;
                ppAPPerms(:,clind) = tmpppAPPerms;
            end
            out.A = A;out.pA = pA;out.ppA = ppA;out.pAR = pAR;out.ppAR = ppAR;out.pAP = pAP;out.ppAP = ppAP;
            obj.A = A;obj.pA = pA;obj.ppA = ppA;obj.pAR = pAR;obj.ppAR = ppAR;obj.pAP = pAP;obj.ppAP = ppAP;
       end
       function out = runPropagatingAnglePermutatedTest(obj,DepVar,COV,t,t2)   % created by Dorothy 2
            ppAP = nan*zeros(1,obj.nLC);
            ppAPPerms = nan*zeros(t2,obj.nLC);
            L1CL1DepVar =  PPMMV2.getL1CL1DepVar(DepVar);
            nS = size(L1CL1DepVar,1);
            forcedPerms = PPMMV2.getForcedPermutations(t2,nS);
            
            % level with no children
            l = obj.nLev;
            clind = zeros(1,2^(l-1));
            tmpppAP = zeros(1,2^(l-1));
            tmpppAPPerms = zeros(t2,2^(l-1));
            parfor cl = 1:1:2^(l-1),
                tmpclind = getClind(obj,l,cl);
                [tmpppAP(1,cl), tmpppAPPerms(:,cl)] = PPMMV2.AnglePermutatedTest_Permutated([getCOV(obj,COV,tmpclind) obj.XREG],getDepVar(obj,DepVar,tmpclind),t,t2,forcedPerms);
                clind(cl) = tmpclind;
            end
            ppAP(1,clind) = tmpppAP;
            ppAPPerms(:,clind) = tmpppAPPerms; 
            
            % levels with children
            for l = obj.nLev-1:-1:1,
                clind = zeros(1,2^(l-1));
                tmpppAP = zeros(1,2^(l-1));
                tmpppAPPerms = zeros(t2,2^(l-1));
                parfor cl = 1:1:2^(l-1),
                    % get the p values for the parent cluster
                    tmpclind = getClind(obj,l,cl);
                    [ppAP_temp, ppAPPerms_temp] = PPMMV2.AnglePermutatedTest_Permutated([getCOV(obj,COV,tmpclind) obj.XREG],getDepVar(obj,DepVar,tmpclind),t,t2,forcedPerms);
                    clind(cl) = tmpclind;
                    % get the p values for the two children clusters
                    children = getChildren(obj,l,cl);
                    clind1 = getClind(obj,l+1,children(1));
                    clind2 = getClind(obj,l+1,children(2));
                    
                    % nonparametric combination of the p values of the parent and the two children
                    [tmpppAP(1,cl), tmpppAPPerms(:,cl)] = PPMMV2.nonParametricCombination ([ppAP_temp,ppAP(clind1),ppAP(clind2)], [ppAPPerms_temp,ppAPPerms(:,clind1),ppAPPerms(:,clind2)]); 
                end
                ppAP(1,clind) = tmpppAP;
                ppAPPerms(:,clind) = tmpppAPPerms;
            end
            out.ppAP = ppAP;
            obj.ppAP = ppAP;
       end 
       function out = assembleNonParCCA(obj,path,precF)
               % PREPARATION
                nS = length(PPMMV2.notNAN(obj.XREG));
                nLC = obj.nLC;
                files = dir(path);
                files = files(3:end);
                nrJobs = length(files);
                Jobs = nan*zeros(1,nrJobs);
                for i=1:1:nrJobs
                    ind = strfind(files(i).name,'_');
                    Jobs(i) = str2double(files(i).name(1:ind(1)-1));
                end
                in = load([path '/' files(Jobs==1).name]);
                in = in.in;
                F = in.F;out.F = F;
                nB = size(in.permF,1);
                nperm = nrJobs*nB;
                if nperm < intmax('uint8')
                   precP = 'uint8';
                elseif nperm < intmax('uint16')
                   precP = 'uint16';
                elseif nperm < intmax('uint32')
                   precP = 'uint32';
                else
                   precP = 'uint64';
                end
                % LOADING DATA
                disp('LOADING DATA');
                permF = zeros(nB,nrJobs,nLC,precF);
                permF(:,1,:) = in.permF; 
                [tmppath,tmpID] = setupParForProgress(nrJobs-1);
                good = ones(1,nrJobs);
                parfor pr=2:1:nrJobs
                   %pr = 2;
                   forname = files(Jobs==pr).name;
                   succ = str2double(forname(end-4));
                   if ~(succ==1), good(pr) = 0; continue; end
                   in = load([path '/' files(Jobs==pr).name]);
                   permF(:,pr,:) = in.in.permF;
                   parfor_progress;
                end
                if sum(good)<nrJobs, out = []; return; end
                permF = reshape(permF,nB*nrJobs,nLC);
                permF = [permF;F];
                closeParForProgress(tmppath,tmpID);
                disp('DONE');
                % P COMPUTATIONS
                disp('COMPUTING P VALUES');
                pF = zeros(nperm+1,nLC,precP);
                [tmppath,tmpID] = setupParForProgress(nLC);
                parfor clind = 1:nLC
                   %clind = 1;
                   pF(:,clind) = PPMMV2.getpFromStat(permF(:,clind),precP);
                   parfor_progress;
                end
                closeParForProgress(tmppath,tmpID);
                clear permF;
                % PROPAGATING P VALUES
                disp('PROPAGATING P VALUES');
                disp(['Level: ' num2str(obj.nLev) ' DONE']);
                out.pF = pF(end,:);
                for l=obj.nLev-1:-1:1
                    %l=5
                    disp(['PROCESSING Level: ' num2str(l)]);
                    [tmppath,tmpID] = setupParForProgress(2^(l-1));
                    for cl = 1:1:2^(l-1)
                        %cl=1
                        clind = uint16(getClind(obj,l,cl));
                        children = getChildren(obj,l,cl);
                        clind1 = uint16(getClind(obj,l+1,children(1)));
                        clind2 = uint16(getClind(obj,l+1,children(2)));
                        pF(:,clind) = PPMMV2.nonParametricCombinationv3(pF(:,[clind clind1 clind2]),precP);
                        parfor_progress;
                    end
                    closeParForProgress(tmppath,tmpID);
                end
                % FINAL OUTPUT
                out.ppF = pF(end,:);
                out.F = PPMMV2.inverseFPrec(out.F,precF);
                out.pF = double(out.pF)/(nperm+1);
                out.ppF = double(out.ppF)/(nperm+1);
                out.nperm = nperm;
                obj.npCCA = out.pF;
                obj.nppCCA = out.ppF;
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
               H = PPMMV2.wilcoxonM(fullM',obj.HTestp);
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
                      Finner = PPMMV2.DOB_SCV_DM(obj.BRIMFold,DepVar,GM,obj.ShapeDistance);
                      TMPIndVar = IndVar(:,end);% Allocate Memory, only last collumn
                      for fi=1:obj.BRIMFold % Fi...
                          FiTestInd = find(Finner==fi);
                          FiTrInd = setdiff(1:n,FiTestInd);
                          [M,H] = robustPLSRegression(obj,IndVar(FiTrInd,:),DepVar(FiTrInd,:),GM(FiTrInd),obj.RegRuns);% Compute BootSampled Regression
                          if obj.HTest, M = H.*M;end% Only keep significant regression coefficients
                          rip = PPMMV2.updateRIP(DepVar(FiTestInd,:),M)';% Get RIP scores
                          if obj.RIPNorm, rip = PPMMV2.normalizeRIP(rip,M,info);end% rescale RIP scores
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
                            [DepVarNew,IndVarNew] = PPMMV2.regUpSample(depvar(balInfo{i}.Ind,:),indvar(balInfo{i}.Ind,:),balInfo{i});
                            SampleDepVar = [SampleDepVar; DepVarNew]; %#ok<*AGROW>
                            SampleIndVar = [SampleIndVar; IndVarNew];
                        end
                    else % Balancing by Decreasing Majority Groups
                        SampleIndVar = [];SampleDepVar = [];
                        for i=1:1:nrgg
                            [DepVarNew,IndVarNew] = PPMMV2.regDownSample(depvar(balInfo{i}.Ind,:),indvar(balInfo{i}.Ind,:),balInfo{i}.nrKeep);
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
                 Mfor(s,:) = PPMMV2.getRegression(SampleIndVar(SampleInd,:),SampleDepVar(SampleInd,:));
             end
           % finalizing output
             M = median(Mfor,1); % extract the median of all regression coefficients
             if ~obj.HTest, H=ones(1,nrD);return; end
             if runs==1, H=ones(1,nrD); return; end
             if nargout<2, return; end
             H = PPMMV2.wilcoxonM(Mfor,obj.HTestp);% wilcoxon test for median
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
           out = PPMMV2.updateRIP(getDepVar(obj,DepVar,clind),getFeatM(obj,clind));
           if obj.RIPNorm, out = PPMMV2.normalizeRIP(out,getFeatM(obj,clind),obj.VarInfo{clind});end
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
                   [out,obj.FSInfo] = PPMMV2.fsCOR(FS,XC,obj.maxFS);
               case 'ROC'
                   [out,obj.FSInfo] = PPMMV2.fsROC(FS,X,obj.maxFS);
               case 'MRMTR'
                   [~,~,H] = PPMMV2.prepMIfs(FS,X,obj.Quantiles);
                   out = PPMMV2.fsMRMTR(H,min(obj.maxFS,size(FS,2)));
                   obj.FSInfo = [];
               case 'SPECCMI'
                   [~,~,H] = PPMMV2.prepMIfs(FS,X,obj.Quantiles);
                   [out,obj.FSInfo] = PPMMV2.fsSPECCMI(H,min(obj.maxFS,size(FS,2)));
                   obj.FSInfo = obj.FSInfo';
               case 'CMIM'
                   [~,~,H] = PPMMV2.prepMIfs(FS,X,obj.Quantiles);
                   out = PPMMV2.fsCMIM(H,min(obj.maxFS,size(FS,2)));
                   obj.FSInfo = [];
               case 'CIFE'
                   [~,~,H] = PPMMV2.prepMIfs(FS,X,obj.Quantiles);
                   out = PPMMV2.fsCIFE(H,min(obj.maxFS,size(FS,2)));
                   obj.FSInfo = [];
               case 'MIXED'
                   nFS = size(FS,2);
                   order = zeros(3,nFS);
                   info = zeros(3,nFS);
                   [order(1,:),info(1,:)] = PPMMV2.fsCOR(FS,XC,nFS);
                   [order(2,:),info(2,:)] = PPMMV2.fsROC(FS,X,nFS);
                   [~,~,H] = PPMMV2.prepMIfs(FS,X,obj.Quantiles);
                   [order(3,:),info(3,:)] = PPMMV2.fsSPECCMI(H,nFS);
                   r = 1:nFS;
                   ranks = repmat(r,3,1);
                   new = info;
                   for i=1:1:3
                       ranks(i,order(i,:)) = r;
                       new(i,order(i,:)) = info(i,:);
                   end
                   info = new;
                   test = mean(ranks,1);
                   [~,out] = sort(test,'ascend');
                   out = out(1:obj.maxFS);
%                    info = info(:,out);
%                    for i=1:1:3
%                        info(i,:) = info(i,:)-min(info(i,:));
%                        info(i,:) = info(i,:)/max(info(i,:));
%                    end
%                    info = sum(info,1);
                   
                     %obj.FSInfo = info(:,out);
%                      obj.FSInfo = info;
                     obj.FSInfo = 1./(1:length(out));
                   
               case 'SVMRFE'
               case 'NONE'
                   out = 1:size(FS,2);
                   obj.FSInfo = [];
               otherwise
                   error('Wrong feature selection procedure')
           end
           obj.FSInd = out;
           %obj.FSInfo = 1./(1:length(out));
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
            obj.Levels = PPMMV2.convertUInt(out);
        end
       function storeClusters(obj)
           out = [];
           for i=1:1:obj.nLev
               out = [out 1:2^(i-1)]; %#ok<AGROW>
           end
           obj.Clusters = PPMMV2.convertUInt(out);
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
       function [out,Beta] = redDepVar(obj,in,COV)
           out = in;
           Beta = in;
           nObs = size(COV,1);
           for c=1:1:obj.nLC
               %c=1;
               Y = in{obj.Levels(c)}.DepVar{obj.Clusters(c)};
               [~,~,~,~,M] = plsregress(COV,Y,size(COV,2));
               Y_est = [ones(nObs,1) COV]*M;
               out{obj.Levels(c)}.DepVar{obj.Clusters(c)} = Y-Y_est;
               Beta{obj.Levels(c)}.DepVar{obj.Clusters(c)} = M;
           end
       end
       function out = redDepVarBeta(obj,in,COV)
           out = in;
           nObs = size(COV,1);
           for c=1:1:obj.nLC
               Y = in{obj.Levels(c)}.DepVar{obj.Clusters(c)};
               Y_est = [ones(nObs,1) COV]*obj.Beta{obj.Levels(c)}.DepVar{obj.Clusters(c)};
               out{obj.Levels(c)}.DepVar{obj.Clusters(c)} = Y-Y_est;
           end
       end
       function out = getFS(obj,FS)
          if ~iscell(FS), out = FS; return; end
          out = getDepVar(obj,FS,find(obj.FSHits));
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
           out = nan*zeros(nrT,1);
            switch obj.ID
                case {1 2 3 4 5 6 7 8 9 10}
                    if isempty(COV), return; end
                    out = COV(:,obj.ID);
                case 100
                    if isempty(GB), return; end
                    out = GB(:,obj.GBID);
                case 1000
                    ind = find(strcmp(obj.RS,RS));
                    if isempty(ind),return;end
                    out = GT(:,ind);
                otherwise
                    return;      
            end 
       end
       function out = List2HierachicalMatrix(obj,in)    % created by Dorothy
            nCl = 2^(obj.nLev-1);
            out = zeros(obj.nLev,nCl);
            for cl = 1:nCl,
                parentEvolution = getParentEvolution(obj,obj.nLev,cl);
                for l = 1:obj.nLev,
                    clind = getClind(obj,l,parentEvolution(l));
                    value = in(clind);
                    if isempty(obj.pCrit);        % fill the matrix with the value
                        out(l,cl) = value;
                    else                             % fill the matrix with 1 of the value is significant
                        if value <= obj.pCrit;
                            out(l,cl) = 1;
                        end
                    end
                end
            end
            % Plot the figure with imagesc(out);
       end
       function out = getParent(obj,level,cluster)              %#ok<INUSL> % Created by Dorothy
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
                out(l) = getParent(obj,l+1,out(l+1));
            end
        end
       function out = getClind(obj,level,cluster)           % Created by Dorothy
            if level > obj.nLev; out = []; return; end
            if cluster > 2^(level-1), out = []; return; end
            out = sum(2.^([1:1:level-1]-1));
            out = out + cluster;
       end
       function out = get.HMatrixppAP(obj)      % created by Dorothy 2
            out = List2HierachicalMatrix(obj,obj.ppAP);
       end
       function imageHits(obj)
           figure;subplot(2,2,1);imagesc(obj.HMatrixpR2);
           title(['pR2 ' num2str(obj.nHitspR2)]);set(gca,'clim',[0 1]);
           subplot(2,2,2);imagesc(obj.HMatrixppR2);
           title(['ppR2 ' num2str(obj.nHitsppR2)]);set(gca,'clim',[0 1]);
           subplot(2,2,3);imagesc(obj.HMatrixppAR);
           title(['ppAR ' num2str(obj.nHitsppAR)]);set(gca,'clim',[0 1]);
           subplot(2,2,4);imagesc(obj.HMatrixppAP);
           title(['ppAP ' num2str(obj.nHitsppAP)]);set(gca,'clim',[0 1]);
       end
       function imageR2Hits(obj)
           figure;subplot(1,3,1);imagesc(obj.HMatrixR2);axis square;%colormap('jet');
           title(['R2 ' num2str(max(obj.R2))]);set(gca,'clim',[0 max(obj.R2)]);
           subplot(1,3,2);imagesc(obj.HMatrixpR2);axis square;
           title(['pR2 ' num2str(obj.nHitspR2)]);set(gca,'clim',[0 1]);
           subplot(1,3,3);imagesc(obj.HMatrixppR2);axis square;
           title(['ppR2 ' num2str(obj.nHitsppR2)]);set(gca,'clim',[0 1]); 
       end
       function imageAHits(obj)
           figure;subplot(3,3,1);imagesc(obj.HMatrixA);axis square;
           title(['A ' num2str(max(obj.A))]);set(gca,'clim',[0 max(obj.A)]);
           %obj.pCrit = [];
           subplot(3,3,2);imagesc(obj.HMatrixpA);axis square;
           %obj.pCrit = 0.05;
           title(['pA ' num2str(obj.nHitspA)]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(3,3,3);imagesc(obj.HMatrixppA);axis square;
           %obj.pCrit = 0.05;
           title(['ppA ' num2str(obj.nHitsppA)]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(3,3,5);imagesc(obj.HMatrixpAR);axis square;
           %obj.pCrit = 0.05;
           title(['pAR ' num2str(obj.nHitspAR)]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(3,3,6);imagesc(obj.HMatrixppAR);axis square;
           %obj.pCrit = 0.05;
           title(['ppAR ' num2str(obj.nHitsppAR)]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(3,3,8);imagesc(obj.HMatrixpAP);axis square;
           %obj.pCrit = 0.05;
           title(['pAP ' num2str(obj.nHitspAP)]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(3,3,9);imagesc(obj.HMatrixppAP);axis square;
           %obj.pCrit = 0.05;
           title(['ppAP ' num2str(obj.nHitsppAP)]);set(gca,'clim',[0 1]);
           
       end
       function imageAHitsRed(obj)
           figure;subplot(1,3,1);imagesc(obj.HMatrixA);axis square;
           title(['A ' num2str(max(obj.A))]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(1,3,2);imagesc(obj.HMatrixpA);axis square;
           %obj.pCrit = 0.05;
           title(['pA ' num2str(obj.nHitspA)]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(1,3,3);imagesc(obj.HMatrixpAR);axis square;
           %obj.pCrit = 0.05;
           title(['pAR ' num2str(obj.nHitspAR)]);set(gca,'clim',[0 1]);
       end
       function imageR2AHits(obj)
           figure;subplot(2,3,1);imagesc(obj.HMatrixR2);axis square;%colormap('jet');
           title(['R2 ' num2str(max(obj.R2))]);set(gca,'clim',[0 max(obj.R2)]);
           subplot(2,3,2);imagesc(obj.HMatrixpR2);axis square;
           title(['pR2 ' num2str(obj.nHitspR2)]);set(gca,'clim',[0 1]);
           subplot(2,3,3);imagesc(obj.HMatrixppR2);axis square;
           title(['ppR2 ' num2str(obj.nHitsppR2)]);set(gca,'clim',[0 1]);
           subplot(2,3,4);imagesc(obj.HMatrixA);axis square;
           title(['A ' num2str(max(obj.A))]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(2,3,5);imagesc(obj.HMatrixpA);axis square;
           %obj.pCrit = 0.05;
           title(['pA ' num2str(obj.nHitspA)]);set(gca,'clim',[0 1]);
           %obj.pCrit = [];
           subplot(2,3,6);imagesc(obj.HMatrixpAR);axis square;
           %obj.pCrit = 0.05;
           title(['pAR ' num2str(obj.nHitspAR)]);set(gca,'clim',[0 1]); 
       end
       function imageTR2A(obj)
           bk = obj.pCrit;
           obj.pCrit = [];
           figure;subplot(3,3,1);imagesc(obj.HMatrixR2);axis square;%colormap('jet');
           title(['R2 ' num2str(max(obj.R2))]);set(gca,'clim',[0 max(obj.R2)]);
           subplot(3,3,2);imagesc(-log10(obj.HMatrixpR2));axis square;
           title(['pR2 ' num2str(max(-log10(obj.pR2)))]);set(gca,'clim',[0 max(-log10(obj.pR2))]);
           subplot(3,3,3);imagesc(-log10(obj.HMatrixppR2));axis square;
           title(['ppR2 ' num2str(max(-log10(obj.pR2)))]);set(gca,'clim',[0 max(-log10(obj.pR2))]);
           subplot(3,3,4);imagesc(obj.HMatrixA);axis square;
           title(['A ' num2str(max(obj.A))]);set(gca,'clim',[0 1]);
           subplot(3,3,5);imagesc(-log10(obj.HMatrixpA));axis square;
           title(['pA ' num2str(max(-log10(obj.pA)))]);set(gca,'clim',[0 max(-log10(obj.pA))]);
           subplot(3,3,6);imagesc(-log10(obj.HMatrixpAR));axis square;
           title(['pAR ' num2str(max(-log10(obj.pAR)))]);set(gca,'clim',[0 max(-log10(obj.pAR))]);
           subplot(3,3,7);imagesc(obj.HMatrixT);axis square;
           title(['T ' num2str(max(obj.AvgT))]);set(gca,'clim',[0 max(obj.AvgT)]);
           subplot(3,3,8);imagesc(obj.HMatrixMaxpT);axis square;
           title(['MpT ' num2str(max(obj.MaxLogpT))]);set(gca,'clim',[0 max(obj.MaxLogpT)]);
           subplot(3,3,9);imagesc(obj.HMatrixFastpT);axis square;
           title(['FpT ' num2str(max(obj.LogpFastT))]);set(gca,'clim',[0 max(obj.LogpFastT)]);
           obj.pCrit = bk;
       end
       function f = imageTR2ACCA(obj)
           bk = obj.pCrit;
           obj.pCrit = [];
           f = figure;subplot(4,3,1);imagesc(obj.HMatrixR2);axis square;%colormap('jet');
           title(['R2 ' num2str(max(obj.R2))]);set(gca,'clim',[0 max(obj.R2)]);
           subplot(4,3,2);imagesc(-log10(obj.HMatrixpR2));axis square;
           title(['pR2 ' num2str(max(-log10(obj.pR2)))]);set(gca,'clim',[0 max(-log10(obj.pR2))]);
           subplot(4,3,3);imagesc(-log10(obj.HMatrixppR2));axis square;
           title(['ppR2 ' num2str(max(-log10(obj.pR2)))]);set(gca,'clim',[0 max(-log10(obj.pR2))]);
           subplot(4,3,4);imagesc(obj.HMatrixA);axis square;
           title(['A ' num2str(max(obj.A))]);set(gca,'clim',[0 0.7]);
           subplot(4,3,5);imagesc(-log10(obj.HMatrixpA));axis square;
           title(['pA ' num2str(max(-log10(obj.pA)))]);set(gca,'clim',[0 max(-log10(obj.pA))]);
           subplot(4,3,6);imagesc(-log10(obj.HMatrixpAR));axis square;
           title(['pAR ' num2str(max(-log10(obj.pAR)))]);set(gca,'clim',[0 max(-log10(obj.pAR))]);
           subplot(4,3,7);imagesc(obj.HMatrixT);axis square;
           title(['T ' num2str(max(obj.AvgT))]);set(gca,'clim',[0 max(obj.AvgT)]);
           subplot(4,3,8);imagesc(obj.HMatrixMaxpT);axis square;
           title(['MpT ' num2str(max(obj.MaxLogpT))]);set(gca,'clim',[0 max(obj.MaxLogpT)]);
           subplot(4,3,9);imagesc(obj.HMatrixFastpT);axis square;
           title(['FpT ' num2str(max(obj.LogpFastT))]);set(gca,'clim',[0 max(obj.LogpFastT)]);
           subplot(4,3,10);imagesc(obj.HMatrixCCA);axis square;
           title(['CCA ' num2str(max(obj.CCA))]);set(gca,'clim',[0 max(obj.CCA)]);
           subplot(4,3,11);imagesc(obj.HMatrixpCCA);axis square;
           title(['pCCA ' num2str(max(-log10(obj.pCCA)))]);set(gca,'clim',[0 max(-log10(obj.pCCA))]);
           obj.pCrit = bk;
       end
       function f = imageTR2ACCAv2(obj)
           bk = obj.pCrit;
           obj.pCrit = [];
           f = figure;subplot(2,4,1);imagesc(obj.HMatrixR2);axis square;%colormap('jet');
           title(['R2 ' num2str(max(obj.R2))]);set(gca,'clim',[0 max(obj.R2)]);
           xlabel([obj.RS],'FontSize',12,'FontWeight','Bold','Color','b');
           subplot(2,4,5);imagesc(-log10(obj.HMatrixppR2));axis square;
           title(['pR2 ' num2str(max(-log10(obj.pR2)))]);set(gca,'clim',[0 max(-log10(obj.pR2))]);
           xlabel(['p11: ' num2str(obj.ppR2(1))],'FontWeight','Bold');
           subplot(2,4,2);imagesc(obj.HMatrixA);axis square;
           title(['A ' num2str(max(obj.A))]);set(gca,'clim',[0 0.7]);
           xlabel(['TI: ' num2str(obj.TestInd)],'FontSize',12,'FontWeight','Bold','Color','b');
           subplot(2,4,6);imagesc(-log10(obj.HMatrixpAR));axis square;
           title(['pA ' num2str(max(-log10(obj.pAR)))]);set(gca,'clim',[0 max(-log10(obj.pAR))]);
           xlabel(['N: ' num2str(sum(obj.pAR<=0.05))],'FontWeight','Bold');
           subplot(2,4,3);imagesc(obj.HMatrixT);axis square;
           title(['T ' num2str(max(obj.AvgT))]);set(gca,'clim',[0 max(obj.AvgT)]);
           xlabel(['mALF: ' num2str(obj.mAlF)],'FontSize',12,'FontWeight','Bold','Color','b');
           subplot(2,4,7);imagesc(obj.HMatrixFastpT);axis square;
           title(['pT ' num2str(max(obj.LogpFastT))]);set(gca,'clim',[0 max(obj.LogpFastT)]);
           xlabel(['N: ' num2str(sum(obj.LogpFastT>=-log10(10^-3)))],'FontWeight','Bold');
           subplot(2,4,4);imagesc(obj.HMatrixCCA);axis square;
           title(['CCA ' num2str(max(obj.CCA))]);set(gca,'clim',[0 max(obj.CCA)]);
           xlabel(['Bal: ' num2str(obj.PosClassBal)],'FontSize',12,'FontWeight','Bold','Color','b');
           subplot(2,4,8);imagesc(obj.HMatrixpCCA);axis square;
           title(['pCCA ' num2str(max(-log10(obj.pCCA)))]);set(gca,'clim',[0 max(-log10(obj.pCCA))]);
           xlabel(['N: ' num2str(sum(obj.pCCA<=10^-3))],'FontWeight','Bold');
           obj.pCrit = bk;
       end
       function f = imageCCA(obj,t)
           bk = obj.pCrit;
           obj.pCrit = [];
           if nargin<2, t= 1000;end
           f = figure;
           setup = [2 2];
           in = List2HierachicalMatrix(obj,obj.CCA);
           subplot(setup(1),setup(2),1);imagesc(in);axis square;
           title(['F ' num2str(max(in(:)))]);set(gca,'clim',[-log10(0.05) max(in(:))]);colorbar;
           if strcmp(obj.Type,'SNPMV2')
              xlabel([obj.VarName ' T: ' num2str(obj.TestInd) ' Bal: ' num2str(obj.PosClassBal)],'FontSize',12,'FontWeight','Bold','Color','b');
           else
              xlabel(obj.VarName,'FontSize',12,'FontWeight','Bold','Color','b');
           end
           
           if~isempty(obj.nppCCA)
               in = List2HierachicalMatrix(obj,-log10(obj.nppCCA));
               subplot(setup(1),setup(2),2);imagesc(in);axis square;
               %title(['nppF ' num2str(obj.nppCCA(1))]);set(gca,'clim',[-log10(0.05) -log10(1/t)]);colorbar;
               title(['nppF ' num2str(obj.nppCCA(1))]);set(gca,'clim',[-log10(0.05) -log10(5*10^-8)]);colorbar;
               xlabel(['N: ' num2str(sum(obj.nppCCA<=0.05))],'FontWeight','Bold');
           end
           
           in = List2HierachicalMatrix(obj,-log10(obj.pCCA));
           subplot(setup(1),setup(2),3);imagesc(in);axis square;
           title(['pF ' num2str(obj.pCCA(1))]);set(gca,'clim',[-log10(0.05) -log10(5*10^-8)]);colorbar;
           xlabel(['N: ' num2str(sum(obj.pCCA<=0.05)) ' Max: ' num2str(max(in(:)))],'FontWeight','Bold');
           
           if ~isempty(obj.npCCA)
               in = List2HierachicalMatrix(obj,-log10(obj.npCCA));
               subplot(setup(1),setup(2),4);imagesc(in);axis square;
               %title(['npF ' num2str(obj.npCCA(1))]);set(gca,'clim',[-log10(0.05) -log10(1/t)]);colorbar;
               title(['npF ' num2str(obj.npCCA(1))]);set(gca,'clim',[-log10(0.05) -log10(5*10^-8)]);colorbar;
               xlabel(['N: ' num2str(sum(obj.npCCA<=0.05))],'FontWeight','Bold');
           end
           obj.pCrit = bk;
       end
       function f = imageA(obj,pCrit)
           if nargin<2 pCrit = 0.05; end
           f = figure;
           setup = [1 2];
           in = List2HierachicalMatrix(obj,obj.A);
           subplot(setup(1),setup(2),1);imagesc(in);axis square;
           title(['A ' num2str(max(in(:)))]);set(gca,'clim',[0 max(in(:))]);colorbar;
           xlabel([obj.RS ' T: ' num2str(obj.TestInd) ' Bal: ' num2str(obj.PosClassBal)],'FontSize',12,'FontWeight','Bold','Color','b');
           
           in = List2HierachicalMatrix(obj,-log10(obj.pAR));
           subplot(setup(1),setup(2),2);imagesc(in);axis square;
           title(['pA ' num2str(max(in(:)))]);set(gca,'clim',[0 -log10(pCrit)]);colorbar;
           xlabel(['N: ' num2str(sum(obj.pAR<=pCrit))],'FontWeight','Bold');
           
       end
       function imageT(obj)
           bk = obj.pCrit;
           obj.pCrit = [];
           figure;subplot(1,3,1);imagesc(obj.HMatrixT);axis square;
           title(['T ' num2str(max(obj.AvgT))]);set(gca,'clim',[0 max(obj.AvgT)]);
           subplot(1,3,2);imagesc(obj.HMatrixMaxpT);axis square;
           title(['MpT ' num2str(max(obj.MaxLogpT))]);set(gca,'clim',[0 max(obj.MaxLogpT)]);
           subplot(1,3,3);imagesc(obj.HMatrixFastpT);axis square;
           title(['FpT ' num2str(max(obj.LogpFastT))]);set(gca,'clim',[0 max(obj.LogpFastT)]);
           obj.pCrit = bk;
       end
       function illustrateHits(obj,RefScan,label)
           index = find(obj.Hits);
           for i=1:1:obj.nHits
               scan = clone(RefScan);
               scan.Material = 'Dull';
               L = obj.Levels(index(i));
               CL = obj.Clusters(index(i));
               val = zeros(1,scan.nrV);
               val(label(L,:)==CL) = 1;
               scan.Value = val;
               scan.ColorMode = 'Indexed';
               v = viewer(scan);
               v.SceneLightVisible = true;
               v.SceneLightLinked = true;
               colormap(v.RenderAxes,'summer');
               set(v.RenderAxes,'clim',[0 1]);
               str = ['L' num2str(obj.Levels(index(i))) ' C' num2str(obj.Clusters(index(i)))];
               v.Tag = str;
           end
       end
       function f = illustrateModValues(obj,RefScan,label,type,th,range)
          if nargin<5, th = []; end
          if nargin<6, range = []; end
          if ischar(type),
             val = eval(['obj.' type]);
          else
             val = type;
          end
          %val = -log(val);
          val = List2Matrix(obj,val);
          if ~isempty(th), val = val<=th; end
          if isempty(range), range = [0 max(val(:))]; end
          scanvalues = nan*zeros(obj.nLev,RefScan.nrV);
          for i=1:1:obj.nLev
              forval = val(:,i);
              lab = label(i,:);
              labc = unique(lab);
              nlabc = length(labc);
              for c=1:1:nlabc
                  %c=1;
                  scanvalues(i,lab==c) = forval(c);
              end
          end
          f = figure;hold on;
          subplot(2,3,1);hold on;camlight HEADLIGHT;axis square;
          set(gca,'clim',range);axis off;colorbar;
          scan = clone(RefScan);
          scan.Value = scanvalues(1,:);
          scan.ColorMode = 'Indexed';
          scan.Axes = gca;
          scan.Visible = true;
          
          subplot(2,3,2);hold on;camlight HEADLIGHT;axis square;
          set(gca,'clim',range);axis off;colorbar;
          scan = clone(RefScan);
          scan.Value = scanvalues(2,:);
          scan.ColorMode = 'Indexed';
          scan.Axes = gca;
          scan.Visible = true;
          
          subplot(2,3,3);hold on;camlight HEADLIGHT;axis square;
          set(gca,'clim',range);axis off;colorbar;
          scan = clone(RefScan);
          scan.Value = scanvalues(3,:);
          scan.ColorMode = 'Indexed';
          scan.Axes = gca;
          scan.Visible = true;
          
          subplot(2,3,4);hold on;camlight HEADLIGHT;axis square;
          set(gca,'clim',range);axis off;colorbar;
          scan = clone(RefScan);
          scan.Value = scanvalues(4,:);
          scan.ColorMode = 'Indexed';
          scan.Axes = gca;
          scan.Visible = true;
          
          subplot(2,3,5);hold on;camlight HEADLIGHT;axis square;
          set(gca,'clim',range);axis off;colorbar;
          scan = clone(RefScan);
          scan.Value = scanvalues(5,:);
          scan.ColorMode = 'Indexed';
          scan.Axes = gca;
          scan.Visible = true;
          
          subplot(2,3,6);hold on;camlight HEADLIGHT;axis square;
          set(gca,'clim',range);axis off;colorbar;
          scan = clone(RefScan);
          scan.Value = scanvalues(6,:);
          scan.ColorMode = 'Indexed';
          scan.Axes = gca;
          scan.Visible = true;

       end
       function illustrateAngleTest(obj,Atest,AngleTable)
             index = find(obj.Hits);
             for i=1:1:obj.nHits
                 L = obj.Levels(index(i));
                 CL = obj.Clusters(index(i));
                 f = figure;hist(Atest{index(i)}.A,100);
                 set(gca,'xlim',[-1 1]);
                 set(gca,'ylim',[0 50]);
                 str = ['L' num2str(L) ' C' num2str(CL) ' Dim: ' num2str(size(Atest{index(i)}.M,1))];
                 title(str);  
                 xlabel(['Avg: ' num2str(Atest{index(i)}.STAT.AvgA)...
                         ' pAR:' num2str(Atest{index(i)}.STAT.pAR)...
                         ' Std:' num2str(Atest{index(i)}.STAT.StdA)...
                         ' U:' num2str(Atest{index(i)}.STAT.UpperA)...
                         ' L:' num2str(Atest{index(i)}.STAT.LowerA)...
                         ' pL:' num2str(PPMMV2.lookUppRA(Atest{index(i)}.STAT.LowerA,size(Atest{index(i)}.M,1),AngleTable))]);
             end

       end
       function out = getPosClassBal(obj,X)
           X(isnan(X)) = [];
           Tind = X==obj.PosClass;
           out = sum(Tind)/length(X);
       end
       function out = compareM(obj,Robj,AngleTable,vis,pCrit)
           A = nan*zeros(1,obj.nLC);
           pA = nan*zeros(1,obj.nLC);
           for i=1:obj.nLC
               %i=1
               A(i) = angle(obj.FeatM{i}',Robj.FeatM{i}');
               pA(i) = PPMMV2.lookUppRA(A(i),length(obj.FeatM{i}),AngleTable);
           end
           out.A = A;
           out.pA = pA;
           if vis
               out.f = figure;
               setup = [1 2];
               in = List2HierachicalMatrix(obj,A);
               subplot(setup(1),setup(2),1);imagesc(in);axis square;
               title(['A ' num2str(max(in(:)))]);set(gca,'clim',[0 max(in(:))]);colorbar;
               xlabel([obj.RS ' T: ' num2str(obj.TestInd) ' Bal: ' num2str(obj.PosClassBal)],'FontSize',12,'FontWeight','Bold','Color','b');

               in = List2HierachicalMatrix(obj,-log10(pA))>=-log10(pCrit);
               subplot(setup(1),setup(2),2);imagesc(in);axis square;colorbar
               title(['N: ' num2str(sum(pA<=pCrit))]);
               set(gca,'clim',[0 max(in(:))]);
               xlabel(num2str(pCrit),'FontWeight','Bold');
           end
       end
       function out = compareCCA(obj,Robj,AngleTable,vis,pCrit)
           A = nan*zeros(1,obj.nLC);
           pA = nan*zeros(1,obj.nLC);
           for i=1:obj.nLC
               %i=1
               A(i) = angle(obj.bCCA{i},Robj.bCCA{i});
               pA(i) = PPMMV2.lookUppRA(A(i),length(obj.bCCA{i}),AngleTable);
           end
           out.A = A;
           out.pA = pA;
           if vis
               out.f = figure;
               setup = [1 2];
               in = List2HierachicalMatrix(obj,A);
               subplot(setup(1),setup(2),1);imagesc(in);axis square;
               title(['A ' num2str(max(in(:)))]);set(gca,'clim',[0 max(in(:))]);colorbar;
               xlabel([obj.RS ' T: ' num2str(obj.TestInd) ' Bal: ' num2str(obj.PosClassBal)],'FontSize',12,'FontWeight','Bold','Color','b');

               in = List2HierachicalMatrix(obj,-log10(pA))>=-log10(pCrit);
               subplot(setup(1),setup(2),2);imagesc(in);axis square;colorbar
               title(['N: ' num2str(sum(pA<=pCrit))]);
               set(gca,'clim',[0 max(in(:))]);
               xlabel(num2str(pCrit),'FontWeight','Bold');
           end
       end
       function out = comparePLSR(obj,Robj,AngleTable,vis,pCrit)
           bk = obj.pCrit;
           obj.pCrit = [];
           A = nan*zeros(1,obj.nLC);
           pA = nan*zeros(1,obj.nLC);
           for i=1:obj.nLC
               %i=1
               A(i) = angle(obj.bPLSR{i}',Robj.bPLSR{i}');
               pA(i) = PPMMV2.lookUppRA(A(i),length(obj.bPLSR{i}),AngleTable);
           end
           out.A = A;
           out.pA = pA;
           if vis
               out.f = figure;
               setup = [1 2];
               in = List2HierachicalMatrix(obj,A);
               subplot(setup(1),setup(2),1);imagesc(in);axis square;
               title(['A ' num2str(max(in(:)))]);set(gca,'clim',[0 max(in(:))]);colorbar;
               xlabel([obj.RS ' T: ' num2str(obj.TestInd) ' Bal: ' num2str(obj.PosClassBal)],'FontSize',12,'FontWeight','Bold','Color','b');

               in = List2HierachicalMatrix(obj,-log10(pA))>=-log10(pCrit);
               subplot(setup(1),setup(2),2);imagesc(in);axis square;colorbar
               title(['N: ' num2str(sum(pA<=pCrit))]);
               set(gca,'clim',[0 max(in(:))]);
               xlabel(num2str(pCrit),'FontWeight','Bold');
           end
           obj.pCrit = bk;
       end
       function out = imageDisRepl(obj,robj,pDis,pRepl,AngleTable)
           
           
           
           
           
           
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
           index = intersect(PPMMV2.notNAN(C),PPMMV2.notNAN(X));
           if ~isempty(C)
            E = PPMMV2.getResiduals(C(index,:),Y(index,:));
            X = PPMMV2.getResiduals(C(index,:),X(index,:));
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
       function [T,pT,pFastT,avgT] = parReg(X,Y,C)
           index = intersect(PPMMV2.notNAN(C),PPMMV2.notNAN(X));
           Y = Y(index,:);
           CX = [C(index,:) X(index)];
           nVar = size(Y,2);
           T = zeros(1,nVar);
           pT = zeros(1,nVar);
           for i=1:1:nVar
               STATS = regstats(Y(:,i),CX,'linear','tstat');
               T(i) = STATS.tstat.t(end);
               pT(i) = STATS.tstat.pval(end);
           end
           pFastT = pfast(pT);
           W = -log10(pT);
           avgT = sum(abs(T).*W)/sum(W);
       end
       function [T,pT,pFastT,avgT,B] = parPLSReg(X,Y,C)
           index = intersect(PPMMV2.notNAN(C),PPMMV2.notNAN(X));
           if ~isempty(C)
            Y = PPMMV2.getResiduals(C(index,:),Y(index,:));
            X = PPMMV2.getResiduals(C(index,:),X(index,:));
           else
            Y = Y(index,:);
            X = X(index,:);
           end
           [nObs,~] = size(Y);
           [~,~,~,~,betha] = plsregress(X,Y,1);
           Y_est = [ones(nObs,1) X]*betha;
           B = betha(end,:);
           %T = betha(2,:)./sqrt(sum((Y_est-Y).^2,1)/(nObs-2))/sqrt(sum((X-mean(X)).^2));
           T = betha(2,:)./(sqrt(sum((Y_est-Y).^2,1)/(nObs-2))/sqrt(sum((X-mean(X)).^2)));
           pT = 2*(1-tcdf(abs(T),nObs-2));
           pT(pT==0) = 10^-32;
           pFastT = pfast(pT);
           W = -log10(pT);
           avgT = sum(abs(T).*W)/sum(W);
       end
       function [F,pF,B] = parCCA(X,Y,C)
           index = intersect(PPMMV2.notNAN(C),PPMMV2.notNAN(X));
           if ~isempty(C)
            Y = PPMMV2.getResiduals(C(index,:),Y(index,:));
            X = PPMMV2.getResiduals(C(index,:),X(index,:));
           else
            Y = Y(index,:);
            X = X(index,:);
           end
           %Y = Y(index,:);
           %CX = [C(index,:) X(index)];
           %CX = X(index);
           [~,B,~,~,~,STATS] = canoncorr(X,Y);
           F = STATS.F(end);
           pF = STATS.pF(end);
       end
       function [ppF, ppFPerms,FTrue] = nonparCCA(X,Y,C,t,forcedPerms)        % Created by Dorothy
            if nargin < 5, forcedPerms = []; end
            index = intersect(PPMMV2.notNAN(C),PPMMV2.notNAN(X));
            if ~isempty(C)
                Y = PPMMV2.getResiduals(C(index,:),Y(index,:));
                X = PPMMV2.getResiduals(C(index,:),X(index,:));
            else
                Y = Y(index,:);
                X = X(index,:);
            end
            % calculate F for the true data
            [~,~,~,~,~,STATS] = canoncorr(X,Y);
            FTrue = single(STATS.F(end));
            %parpF = STATS.pF(end);
            if t==0, ppF = nan; return; end
            % calculate F for the permutated data
            FPerm = zeros(t,1);
            %parpFPerms = zeros(t,1);
            parfor p=1:1:t
                if isempty(forcedPerms)
                    perm = randperm(length(index));
                else
                    perm = forcedPerms(p,:);
                    perm(perm>length(index)) = [];
                end
                [~,~,~,~,~,STATS] = canoncorr(X,Y(perm,:));
                FPerm(p) = single(STATS.F(end));
                %parpFPerms(p) = STATS.pF(end);
            end
            % calculate the p-value for the true data and the permutated data
            
            ppF = single(sum([FPerm;FTrue] >= repmat(FTrue,t+1,1),1)/(t+1));
            ppFPerms = single(sum(repmat([FPerm;FTrue]',t,1) >= repmat(FPerm,1,t+1),2)./(t+1));  
            
            %ppF = sum([FPerm;FTrue] >= repmat(FTrue,t+1,1),1)/(t+1);
            %if ppF==(1/t+1), ppF = parpF; end
            %ppFPerms = sum(repmat([FPerm;FTrue]',t,1) >= repmat(FPerm,1,t+1),2)./(t+1);
            %index = find(ppFPerms==(1/t+1));
            %if isempty(index), return; end
            %ppFPerms(index) = parpFPerms(index);
       end
       function [ppF, ppFPerms,FTrue] = nonparCCAv2(X,Y,C,t,str)        % Created by Dorothy
            index = intersect(PPMMV2.notNAN(C),PPMMV2.notNAN(X));
            if ~isempty(C)
                Y = PPMMV2.getResiduals(C(index,:),Y(index,:));
                X = PPMMV2.getResiduals(C(index,:),X(index,:));
            else
                Y = Y(index,:);
                X = X(index,:);
            end
            % calculate F for the true data
            [~,~,~,~,~,STATS] = canoncorr(X,Y);
            FTrue = single(STATS.F(end));
            %parpF = STATS.pF(end);
            if t==0, ppF = nan; return; end
            % calculate F for the permutated data
            FPerm = zeros(t,1,'single');
            %parpFPerms = zeros(t,1);
            parfor p=1:1:t
                rng(p);% Fix permutations
                perm = randperm(length(index));
                [~,~,~,~,~,STATS] = canoncorr(X,Y(perm,:));
                FPerm(p) = single(STATS.F(end));
                %parpFPerms(p) = STATS.pF(end);
            end
            % calculate the p-value for the true data and the permutated data
            switch str
                case 'uint8'
                    ppF = uint8(sum(FPerm>=FTrue)+1);
                    ppFPerms = zeros(t,1,str);
                    parfor p=1:t
                        forFTrue = FPerm(p);
                        ppFPerms(p) = uint8(sum(FPerm>=forFTrue)+(FTrue>=forFTrue));
                    end 
                case 'uint16'
                    ppF = uint16(sum(FPerm>=FTrue)+1);
                    ppFPerms = zeros(t,1,str);
                    parfor p=1:t
                        forFTrue = FPerm(p);
                        ppFPerms(p) = uint16(sum(FPerm>=forFTrue)+(FTrue>=forFTrue));
                    end
                case 'uint32'
                    ppF = uint32(sum(FPerm>=FTrue)+1);
                    ppFPerms = zeros(t,1,str);
                    parfor p=1:t
                        forFTrue = FPerm(p);
                        ppFPerms(p) = uint32(sum(FPerm>=forFTrue)+(FTrue>=forFTrue));
                    end
                case 'uint64'
                    ppF = uint64(sum(FPerm>=FTrue)+1);
                    ppFPerms = zeros(t,1,str);
                    parfor p=1:t
                        forFTrue = FPerm(p);
                        ppFPerms(p) = uint64(sum(FPerm>=forFTrue)+(FTrue>=forFTrue));
                    end
            end
       end
       function out = singleCCA(X,Y,precision)
                [~,~,~,~,~,STATS] = canoncorr(X,Y);
                switch precision
                    case 'uint8'
                        out = uint8(STATS.F(end)*10);
                    case 'uint16'
                        out = uint16(STATS.F(end)*100);
                    otherwise
                        out = STATS.F(end);
                end    
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
                       forM(:,k) = PPMMV2.getRegression(IndVar(ind,:),DepVar(ind,:));
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
               STAT.pAR = PPMMV2.lookUppRA(STAT.AvgA,ShapeDim,AngleTable);
               [STAT.FA,STAT.pFA] = PPMMV2.lookUpFA(A,ShapeDim,AngleTable);
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
                %[~,~,~,~,~,~,~,stats] = plsregress(X,Y,min(size(X,2),size(Y,2)));
                [~,~,~,~,M] = plsregress(X,Y,min(size(X,2),size(Y,2)));
                Y_est = [ones(size(X,1),1) X]*M;
                out = Y-Y_est;
                %out = stats.Yresiduals;
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
       function [ppR2, ppR2Perms,R2True] = testPartialPLSR_Permutated(X,Y,C,t,forcedPerms)        % Created by Dorothy
            if nargin < 5, forcedPerms = []; end
            index = intersect(PPMMV2.notNAN(C),PPMMV2.notNAN(X));
            %diff = setdiff(1:length(X),index);
            %E = nan*zeros(size(Y));
            if ~isempty(C)
                E = PPMMV2.getResiduals(C(index,:),Y(index,:));
                X = PPMMV2.getResiduals(C(index,:),X(index,:));
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
                if isempty(forcedPerms)
                    perm = randperm(length(index));
                else
                    perm = forcedPerms(p,:);
                    perm(perm>length(index)) = [];
                    %perm = intersect(forcedPerms(p,:),index,'stable');
                    
                end
                [~,~,~,~,~,forvar] = plsregress(X,E(perm,:),1);
                R2Perm(p) = forvar(2);
            end
            % calculate the p-value for the true data and the permutated data
            ppR2 = sum([R2Perm;R2True] >= repmat(R2True,t+1,1),1)/(t+1);
            ppR2Perms = sum(repmat([R2Perm;R2True]',t,1) >= repmat(R2Perm,1,t+1),2)./(t+1);
       end
       function [pA, pAPerm, pAR, pARPerm,pAP,pAPPerm,ATrue] = AngleTest_Permutated(IndVar,DepVar,AngleTable,t,maxM,Aval,t2,forcedPerms)    % Created by Dorothy (long computation time)
            % calculate the p-value for the true data
            [A,~,STAT] = PPMMV2.AngleTest(IndVar,DepVar,AngleTable,t,maxM,Aval);
            ATrue = mean(A);
            pA = STAT.pA;
            pAR = STAT.pAR;       
            % calculate the p-value for the permutated data
            pAPerm = zeros(t2,1);
            pARPerm = zeros(t2,1);
            APerm = zeros(t2,1);
            for p = 1:1:t2; 
                if nargin < 8 || isempty(forcedPerms)==1;
                    perm = randperm(size(IndVar,1));
                else
                    perm = forcedPerms(p,:);
                end
                [Afor,~,STAT] = PPMMV2.AngleTest(IndVar(perm,:),DepVar,AngleTable,t,maxM,Aval);
                pAPerm(p) = STAT.pA;
                pARPerm(p) = STAT.pAR;
                APerm(p) = mean(Afor);
            end
            pAP = sum([APerm;ATrue] >= repmat(ATrue,t2+1,1),1)/(t2+1);
            pAPPerm = sum(repmat([APerm;ATrue]',t2,1) >= repmat(APerm,1,t2+1),2)./(t2+1);
       end
       function [pAP, pAPPerm] = AnglePermutatedTest_Permutated(IndVar,DepVar,t,t2,forcedPerms)     % Created by Dorothy 2
           if nargin < 5, forcedPerms = []; end
           A = nan*zeros(t,1);                         % stores the angle of the different 2-group-divitions
           APerm = nan*zeros(t2,1);                    % stores the angle of the different permutations
           [nS,ShapeDim] = size(DepVar);    
           % calculate the 'stable' p-value for the true data
           for i=1:t
               F = crossvalind('Kfold',nS,2);
               forM = nan*zeros(ShapeDim,2);
               for k=1:1:2,
                   ind = find(F==k);
                   forM(:,k) = PPMMV2.getRegression(IndVar(ind,:),DepVar(ind,:));
               end
               A(i) = angle(forM(:,1),forM(:,2));
           end
           ATrue = mean(A);
           
           % calculate the 'unstable' p-values for the permutated data
           for p = 1:1:t2,
               if isempty(forcedPerms)
                   perm = randperm(size(IndVar,1));
               else
                   perm = forcedPerms(p,:);
               end
               IndVarPerm = IndVar(perm,:);
               F = crossvalind('Kfold',nS,2);
               forM = nan*zeros(ShapeDim,2);
               for k=1:1:2,
                   ind = find(F==k);
                   forM(:,k) = PPMMV2.getRegression(IndVarPerm(ind,:),DepVar(ind,:));
               end
               APerm(p) = angle(forM(:,1),forM(:,2));
           end
           
           % calculate the p-value for the true data and the permutated data
           pAP = sum([APerm;ATrue] >= repmat(ATrue,t2+1,1),1)/(t2+1);
           pAPPerm = sum(repmat([APerm;ATrue]',t2,1) >= repmat(APerm,1,t2+1),2)./(t2+1);
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
                 tmp = PPMMV2.getEER(FS(Tind,f),FS(Find,f));
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
            %[~,~,H] = PPMMV2.prepMIfs(FS,X,qs);
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
            %[~,~,H] = PPMMV2.prepMIfs(FS,X,qs);
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
            %[~,~,H] = PPMMV2.prepMIfs(FS,X,qs);
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
            %[~,~,H] = PPMMV2.prepMIfs(FS,X,qs);
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
       function [EER,G,AUC,x,y,TH,Y,pAUC] = getEER(tmatches,fmatches)
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
                se = PPMMV2.standardError(AUC,na,nn);
                pAUC = normpdf((AUC-0.5)/se,0,1);
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
       function se = standardError(a,na,nn)
           q1 = a/(2-a);
           a2 = a*a;
           q2 = 2*a2/(1+a);
           se = sqrt( (a-a2 + (na-1)*(q1-a2) + (nn-1)*(q2-a2))/(na*nn) );    
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
           M = max(in(:));
           if M<=intmax('uint8'),out = uint8(in);return;end
           if M<=intmax('uint16'), out = uint16(in);return;end
           if M<=intmax('uint32'), out = uint32(in);return;end
           out = uint64(in);
       end
       function out = getForcedPermutations(t,nS,runningindex)           % Created by Dorothy
            if nargin<3,runningindex=1:t;end
            if nS<=intmax('uint8')
               str = 'uint8';
            elseif nS<=intmax('uint16')
               str = 'uint16';
            elseif nS<=intmax('uint32')
               str = 'uint32';
            else
               str = 'uint64';
            end
            out = zeros(t,nS,str);
            for p = 1:1:t,
                rng(runningindex(p));
                out(p,:) = PPMMV2.convertUInt(randperm(nS));
            end
       end
       function [p, pPerms] = nonParametricCombination(p2combine, pPerms2combine)     % Created by Dorothy
            t = size(pPerms2combine,1);
            % Calculate the combined statistics: "fisher"
            stat = -1* sum(log(p2combine),2);
            statPerm = -1* sum(log(pPerms2combine),2);
            % Calculate the p_values of the true data and the permutated data 
            p = sum([statPerm;stat] >= repmat(stat,t+1,1),1)/(t+1);
            if nargout >= 2,
                pPerms = sum(repmat([statPerm;stat]',t,1) >= repmat(statPerm,1,t+1),2)./(t+1);
            end
       end
       function [p, pPerms] = nonParametricCombinationv2(p2combine, pPerms2combine,str)     % Created by Dorothy
            t = size(pPerms2combine,1);
            % Calculate the combined statistics: "fisher
            stat = single(-1*sum(log(single(p2combine)/(t+1)),2));
            statPerm = single(-1*sum(log(single(pPerms2combine)/(t+1)),2));
            %statPerm = zeros(t,1,'single');
            %parfor i=1:t;
            %    pPerm2c = single(pPerms2combine(i,:))/(t+1);
            %    statPerm(i) = single(-1*sum(log(pPerm2c),2));
            %end
            % Calculate the p_values of the true data and the permutated data
            switch str
                case 'uint8'
                    p = uint8(sum(statPerm>=stat)+1);
                    pPerms = zeros(t,1,str);
                    parfor i=1:t
                        forstat = statPerm(i);
                        pPerms(i) = uint8(sum(statPerm>=forstat)+(stat>=forstat));
                    end 
                case 'uint16'
                    p = uint16(sum(statPerm>=stat)+1);
                    pPerms = zeros(t,1,str);
                    parfor i=1:t
                        forstat = statPerm(i);
                        pPerms(i) = uint16(sum(statPerm>=forstat)+(stat>=forstat));
                    end
                case 'uint32'
                    p = uint32(sum(statPerm>=stat)+1);
                    pPerms = zeros(t,1,str);
                    parfor i=1:t
                        forstat = statPerm(i);
                        pPerms(i) = uint32(sum(statPerm>=forstat)+(stat>=forstat));
                    end
                case 'uint64'
                    p = uint64(sum(statPerm>=stat)+1);
                    pPerms = zeros(t,1,str);
                    parfor i=1:t
                        forstat = statPerm(i);
                        pPerms(i) = uint64(sum(statPerm>=forstat)+(stat>=forstat));
                    end
            end                 
%             switch str
%                 case 'uint8'
%                     p = uint8(sum([statPerm;stat] >= repmat(stat,t+1,1),1));
%                     if nargout >= 2,
%                         pPerms = uint8(sum(repmat([statPerm;stat]',t,1) >= repmat(statPerm,1,t+1),2));
%                     end 
%                 case 'uint16'
%                     p = uint16(sum([statPerm;stat] >= repmat(stat,t+1,1),1));
%                     if nargout >= 2,
%                         pPerms = uint16(sum(repmat([statPerm;stat]',t,1) >= repmat(statPerm,1,t+1),2));
%                     end 
%                 case 'uint32'
%                     p = uint32(sum([statPerm;stat] >= repmat(stat,t+1,1),1));
%                     if nargout >= 2,
%                         pPerms = uint32(sum(repmat([statPerm;stat]',t,1) >= repmat(statPerm,1,t+1),2));
%                     end 
%                 case 'uint64'
%                     p = uint64(sum([statPerm;stat] >= repmat(stat,t+1,1),1));
%                     if nargout >= 2,
%                         pPerms = uint64(sum(repmat([statPerm;stat]',t,1) >= repmat(statPerm,1,t+1),2));
%                     end 
%             end
       end
       function out = nonParametricCombinationv3(p2combine,precP)     % Created by Dorothy
            t = size(p2combine,1);
            % Calculate the combined statistics: "fisher
            stat = single(-1*sum(log(single(p2combine)/(t)),2));
            stat = stat*100;
            M = ceil(max(stat));
            if M < intmax('uint8')
                 stat = uint8(stat);
            elseif M < intmax('uint16')
                 stat = uint16(stat);
            elseif M < intmax('uint32')
                 stat = uint32(stat);
            else
                 stat = uint64(stat);
            end
            out = PPMMV2.getpFromStat(stat,precP); 
%             out = zeros(t,1,precP);
%             [tmppath,tmpID] = setupParForProgress(10000);
%             parfor pr=1:10000
%                out(pr) = sum(stat>=stat(pr));
%                parfor_progress;
%             end
%             closeParForProgress(tmppath,tmpID);
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
       function out = inverseFPrec(in,precision)
                switch precision
                    case 'uint8'
                        out = double(in)/10;
                    case 'uint16'
                        out = double(in)/100;
                    otherwise
                        out = in;
                end
       end
       function out = applyUintPrec(in,precision)
                switch precision
                    case 'uint8'
                        out = uint8(in);
                    case 'uint16'
                        out = uint16(in);
                    case 'uint32'
                       out = uint32(in);
                    case 'uint64'
                       out = uint64(in);
                    otherwise
                        out = in;
                end
       end
       function out = getpFromStat(in,precision)
                [U,~,IU] = unique(in);
                nU = length(U);
                tmp = zeros(nU,1);
                parfor u=1:nU
                    %u=1;
                    val = sum(in>=U(u));
                    tmp(u) = PPMMV2.applyUintPrec(val,precision);
                end
                out = tmp(IU);
       end
    end
end