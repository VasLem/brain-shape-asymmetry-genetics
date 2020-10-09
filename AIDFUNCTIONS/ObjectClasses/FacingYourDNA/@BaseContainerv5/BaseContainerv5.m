classdef BaseContainerv5 < superClass
   % A container class containing all the parameter values
   properties
      TrackID = 'X';
      ShapeSpace = [];
      RedShapeSpace = [];
      RedShapeSpaceType = 'RIP'; % RIP (using ripped values) or X (using original values)
      COV = [];
      COVNames = {};
      GB = [];% genetic background axes
      Method = 'PLSR'; % PLSR BRIM or FOLDED BRIM
      Match = 'QD'; % kind of matching score (QD = quadratic or likelyhood ratio, INLIER = bayesian inlier belief)
      Classify = 'SOFT';
      CV = [];
      EV = [];
      CondType = 'X';
   end
   properties (Dependent = true)
       GBNames;
       COVRIP;
       GBRIP;
       COVCond;
       GBCond;
   end
   properties (Dependent = true, Hidden = true)
      nrCOV;
      nrGB;
      nrCOVGB;
      DepVar;
      OrigDepVar;
      RedDepVar;
      OrigRedDepVar;
      AvgCOV;
      AvgGB;
      AvgCOVRIP;
      AvgGBRIP;
      n;
      ShapeDim;
      RedShapeDim;
      indCOV;
   end
   properties (Hidden = true)
       GBindex = [];
   end
   % STAGE 1: FINDING RELEVANT GB AXES, REGRESSION INDEPENDENT
   properties (Hidden = true)
       ST1R2 = [];
       ST1R2P = [];
       ST1R2PT = 0.05;
   end
   properties (Hidden = true, Dependent = true)
       ST1Code;
   end
   % STAGE 2: FINDING RELEVANT GB AXES, REGRESSION DEPENDENT
   properties (Hidden = true)
       ST2R2 = [];
       ST2R2P = [];
       ST2R2PT = 0.05;
   end
   properties (Hidden = true, Dependent = true)
       ST2Code;
   end
   % STAGE 3: FINDING RELEVANT GB AXES, NGBRIM
   properties (Hidden = true)
       ST3P = [];
       ST3Stat = [];
       ST3PT = 0.001;
   end
   properties (Hidden = true, Dependent = true)
       ST3Code;
       ST4Ref;
   end
   properties
       REF = [];
       COVMOD = [];
       GBMOD = [];
   end
   properties %(Hidden = true)% (NGBRIM)
      PartitionDM = [];% Distance Matrix to drive Fold partitioning
      PartitionGM = [];% Group membership to drive Fold partitioning
      RegRuns = 50; % Number of Regression runs (NGBRIM)
      SamplePercentage = 1;% Percentage of data to sample in each run
      OuterFold = 10; % number of Outer Folds (NGBRIM)
      InnerFold = 10; % number of Inner Folds (NGBRIM)
      MaxIterations = 5; % Maximum number of BRIM iterations (NGBRIM)
      Htest = true; % Perform Wilcoxon test om partial regression coefficients (NGBRIM)
      HtestP = 0.001;% significance level for Wilcoxon rank test of regression coefficients
      RIPNormalize = true; % normalize RIP values based on group averages (NGBRIM)
      StopCorr = 0.98; % Stopping correlation between subsequent iterations (NGBRIM)
      RIPStatPerm = 1000; % Number of permutations in testing significance of RIP values (NGBRIM)
      ParOuter = true;% execute parallel computation on outer fold or not (inner fold then) (NGBRIM)
   end
   methods % Constructor
        function obj = BaseContainerv5(varargin)
            obj = obj@superClass(varargin{:});         
        end
   end
   methods % GENERAL GETTING/SETTING
       function out = get.GBNames(obj)
           out = cell(1,obj.nrGB);
           for i=1:1:obj.nrGB
               out{i} = ['PC' num2str(obj.GBindex(i))];
           end
       end
       function out = get.nrCOV(obj)
                out = size(obj.COV,2); 
       end
       function out = get.indCOV(obj)
                if isempty(obj.COV), out = []; return; end
                index = (1:size(obj.COV,1));
                [i,~] = find(isnan(obj.COV));
                i = unique(i);out = setdiff(index,i);
       end
       function out = get.nrGB(obj)
                out = size(obj.GB,2); 
       end
       function out = get.AvgCOV(obj)
           if isempty(obj.COV), out = []; return; end
           out = nanmean(obj.COV);
       end
       function out = get.AvgGB(obj)
           if isempty(obj.GB), out = []; return; end
           out = nanmean(obj.GB);
       end
       function out = get.AvgCOVRIP(obj)
           if isempty(obj.COVRIP), out = []; return; end
           out = nanmean(obj.COVRIP);
       end
       function out = get.AvgGBRIP(obj)
           if isempty(obj.GBRIP), out = []; return; end
           out = nanmean(obj.GBRIP);
       end
       function out = get.ShapeSpace(obj)
           out = obj.ShapeSpace;
           if ~superClass.isH(out), out = []; end
       end
       function out = get.RedShapeSpace(obj)
           out = obj.RedShapeSpace;
           if ~superClass.isH(out), out = []; end
       end
       function out = get.n(obj)
          if isempty(obj.ShapeSpace), out = 0; return; end
          out = obj.ShapeSpace.n;
       end 
       function out = get.DepVar(obj)
            if isempty(obj.ShapeSpace), out = [];return; end
            out = obj.ShapeSpace.Tcoeff./repmat(obj.ShapeSpace.EigStd',obj.n,1);
       end
       function out = get.RedDepVar(obj)
            if isempty(obj.RedShapeSpace), out = [];return; end
            out = obj.RedShapeSpace.Tcoeff./repmat(obj.RedShapeSpace.EigStd',obj.n,1);
       end
       function out = get.OrigDepVar(obj)
            if isempty(obj.ShapeSpace), out = [];return; end
            out = obj.ShapeSpace.Tcoeff;
       end
       function out = get.OrigRedDepVar(obj)
            if isempty(obj.RedShapeSpace), out = [];return; end
            out = obj.RedShapeSpace.Tcoeff;
       end
       function out = get.ShapeDim(obj)
           if isempty(obj.ShapeSpace), out = 0; return; end
           out = obj.ShapeSpace.nrEV;
       end
       function out = get.RedShapeDim(obj)
           if isempty(obj.RedShapeSpace), out = 0; return; end
           out = obj.RedShapeSpace.nrEV;
       end
       function out = get.COVRIP(obj)
           if isempty(obj.COVMOD), out = []; return; end
           out = obj.COVMOD.RIP;
       end
       function out = get.GBRIP(obj)
           if isempty(obj.MOD), out = []; return; end
           out = obj.MOD.RIP(:,obj.nrCOV+1:end);
       end
       function out = get.COVCond(obj)
          switch obj.CondType
              case 'X'
                  out = obj.COV;
              case 'RIP'
                  out = obj.COVRIP;
          end
       end
       function out = get.GBCond(obj)
          switch obj.CondType
              case 'X'
                  out = obj.GB;
              case 'RIP'
                  out = obj.GBRIP;
          end
       end
       function out = get.nrCOVGB(obj)
           out = obj.nrCOV+obj.nrGB;
       end
       function out = get.ST4Ref(obj)
          out = obj.REF; 
       end
   end
   methods % STAGE 1 GETTING/SETTING
       function out = get.ST1Code(obj)
                if isempty(obj.ST1R2P), out = []; return; end
                out = obj.ST1R2P<=obj.ST1R2PT;
       end
   end
   methods % STAGE 2 GETTING/SETTING
       function out = get.ST2Code(obj)
                if isempty(obj.ST2R2P), out = []; return; end
                out = obj.ST2R2P<=obj.ST2R2PT;
       end
   end
   methods % STAGE 3 GETTING/SETTING
       function out = get.ST3Code(obj)
                if isempty(obj.ST3P), out = []; return; end
                out = obj.ST3P<=obj.ST3PT;
       end
   end
   methods % STAGE 4 GETTING/SETTING
   end
   methods % STAGE INTERFACE FUNCTIONS    
       function runStage1(obj,t)
           % traditional stage 1, using the partial regression statistic to
           % see if there is a significant effect or not
           E = BaseContainerv5.getResiduals(obj.COV(obj.indCOV,:),obj.ShapeSpace.Tcoeff(obj.indCOV,:));
           ind = obj.indCOV;COV = obj.COV(ind,:);gb = obj.GB(ind,:);nrGB = obj.nrGB; %#ok<*PROP>
           [path,ID] = setupParForProgress(nrGB);
           r = nan*zeros(1,nrGB);rp = nan*zeros(1,nrGB);n = length(ind);
           parfor i=1:1:nrGB
               T = gb(:,i);T = BaseContainerv5.getResiduals(COV,T);
               [~,~,~,~,~,pctvar] = plsregress(T,E,1);R2Part = pctvar(2);
               R2PartCount = false(t,1);
               for p=1:1:t
                   pind = randperm(n);
                   [~,~,~,~,~,pvar] = plsregress(T,E(pind,:),1);
                   R2PartCount(p) = (pvar(2)>=R2Part);
               end
               r(i) = R2Part;rp(i) = sum(R2PartCount)/t;
               parfor_progress;
           end
           closeParForProgress(path,ID);
           obj.ST1R2 = r;obj.ST1R2P = rp;
           %index = find(obj.ST1Code);
           %reduceGB(obj,index);
       end
       function runStage2(obj,t)
           ind = obj.indCOV;COV = obj.COV(ind,:);E = obj.ShapeSpace.Tcoeff(ind,:);
           gb = obj.GB(ind,:);nrGB = obj.nrGB;all = 1:nrGB; %#ok<*PROP>
           [path,ID] = setupParForProgress(nrGB);
           r = nan*zeros(1,nrGB);rp = nan*zeros(1,nrGB);n = length(ind);
           parfor i=1:1:nrGB
               tr = setdiff(all,i);
               T = BaseContainerv5.getResiduals([COV, gb(:,tr)],gb(:,i));
               ET = BaseContainerv5.getResiduals([COV, gb(:,tr)],E);
               [~,~,~,~,~,pctvar] = plsregress(T,ET,1);R2Part = pctvar(2);
               R2PartCount = false(t,1);
               for p=1:1:t
                   pind = randperm(n);
                   [~,~,~,~,~,pvar] = plsregress(T,ET(pind,:),1);
                   R2PartCount(p) = (pvar(2)>=R2Part);
               end
               r(i) = R2Part;rp(i) = sum(R2PartCount)/t;
               parfor_progress;
           end
           closeParForProgress(path,ID);
           obj.ST2R2 = r;obj.ST2R2P = rp;
           index = find(obj.ST2Code);
           if length(index)==nrGB; disp('All Survived'); return; end
           % reduce and do it again, untill all survive
           disp('Reducing and Running Again');
           reduceGB(obj,index);runStage2(obj,t);
       end
       function runStage3(obj,K,t)
           crossValidate(obj,K,t);
           obj.ST3P = obj.CV.pfast(obj.nrCOV+1:end);
           obj.ST3Stat = mean(obj.CV.Stat(obj.nrCOV+1:end,:),2)';
           index = find(obj.ST3Code);
           if length(index)==obj.nrGB; disp('All Survived');return;end
           % reduce and do it again, untill all survive
           disp('Reducing and Running Again');
           reduceGB(obj,index);runStage3(obj,K,t);
       end
       function runStage4(obj)
            build(obj); 
          % Setting up reference scan 
            obj.REF.Scan = clone(obj.ShapeSpace.Average);obj.REF.Coeff = (obj.ShapeSpace.AvgCoeff./obj.ShapeSpace.EigStd)';
            obj.REF.AvgCOV = obj.AvgCOV;obj.REF.AvgGB = obj.AvgGB;
            RIPref = BaseContainerv5.updateRIP(obj.REF.Coeff,obj.MOD.M)';
            if obj.RIPNormalize, RIPref = BaseContainerv5.normalizeRIP(RIPref,obj.MOD.M,obj.MOD.Var); end
            obj.REF.COVRIP = RIPref(1:obj.nrCOV);obj.REF.GBRIP = RIPref(obj.nrCOV+1:end);
           % Obtaining regressions to contruct morphs
            [obj.MOD.MX,obj.MOD.MIX] = getRegression([obj.COV, obj.GB],obj.DepVar,(1:size([obj.COV, obj.GB],2)));
            [obj.MOD.MRIP,obj.MOD.MIRIP] = getRegression(obj.MOD.RIP,obj.DepVar,(1:size(obj.MOD.RIP,2)));
           % Creating Reduced Space
            createReducedSpace(obj);
            obj.REF.RedScan = clone(obj.RedShapeSpace.Average);obj.REF.RedCoeff = (obj.RedShapeSpace.AvgCoeff./obj.RedShapeSpace.EigStd)';
       end
       function out = runStage1PPM(obj,K1,K2,t)
          % instead of the traditional stage 1, do something similar (axis per axis), in which the data is split and tested accross the folds in classification ability 
          ind = obj.indCOV;COV = obj.COV(ind,:);gb = obj.GB(ind,:);nrGB = obj.nrGB;
          [path,ID] = setupParForProgress(nrGB);       
          parfor i=1:1:nrGB;
              for k1=1:1:K1
                 %k1 = 1;i=1;

                 F = DOB_SCV_DM(K2,obj.PartitionDM,obj.PartitionGM); 




              end
          end
          closeParForProgress(path,ID); 
           
       end
   end
   methods % HIGH LEVEL INTERFACE FUNCTIONS
       function out = crossvalidateCOV(obj,K1,K2,t)
           FD = squareform(pdist(obj.DepVar,'euclidean'));% Facial Similarity Matrix
           allind = 1:obj.n;
           folds = cell(K1,K2);
           for k1=1:1:K1
               % k1=1;
               %F = DOB_SCV_DM(K2,FD,obj.COV(:,1));% Partition according to Sex (grouping) and facial similarity (Continuous)
               F = DOB_SCV_DM(K2,obj.n,[]);% Partition according to Sex (grouping) and facial similarity (Continuous)
               for k2=1:1:K2
                  % k2= 1; 
                  TestInd = find(F==k2);TrInd = setdiff(allind,TestInd);
                 % create sub object 
                  folds{k1,k2}.obj = reduceSamples(obj,TrInd);
                 % build the cov model 
                  buildCOV(folds{k1,k2}.obj);
                 % cross fold RIP statistics
                  [folds{k1,k2}.Stat,folds{k1,k2}.pStat] = COVRIPStatistics(folds{k1,k2}.obj,obj.DepVar(TestInd,:),obj.COV(TestInd,:),t);
                 % cross fold facial sex biometrics
                  [folds{k1,k2}.EER,folds{k1,k2}.R] = COVSubBiometrics(obj,obj.DepVar(TestInd,:),obj.COV(TestInd,:));
               end
           end
           EER = nan*zeros(K1,K2,obj.nrCOV);
           R1 = nan*zeros(K1,K2,obj.nrCOV);
           R10 = nan*zeros(K1,K2,obj.nrCOV);
           R20 = nan*zeros(K1,K2,obj.nrCOV);
           STAT = nan*zeros(K1,K2,obj.nrCOV);
           pSTAT = nan*zeros(K1,K2,obj.nrCOV);
           M = nan*zeros(K1,K2,K2-1,obj.nrCOV);
           allK2 = 1:K2;
           for k1=1:K1
               for k2=1:K2
                   for i=1:1:obj.nrCOV
                        EER(k1,k2,i) = folds{k1,k2}.EER(i);
                        R1(k1,k2,i) = folds{k1,k2}.R(i,1);
                        R10(k1,k2,i) = folds{k1,k2}.R(i,2);
                        R20(k1,k2,i) = folds{k1,k2}.R(i,3);
                        STAT(k1,k2,i) = folds{k1,k2}.Stat(i);
                        pSTAT(k1,k2,i) = folds{k1,k2}.pStat(i);
                   end
               end
               for i=1:1:obj.nrCOV
                   for k2=1:K2
                       remK2 = setdiff(allK2,k2);
                       for k3=1:K2-1
                           M(k1,k2,k3,i) = angle(folds{k1,k2}.obj.COVMOD.M(i,:)',folds{k1,remK2(k3)}.obj.COVMOD.M(i,:)');
                       end
                   end
               end
           end
           out.EER = EER;out.R1 = R1;out.R10 = R10;out.R20 = R20;out.STAT = STAT;out.pSTAT = pSTAT;out.M = M;
       end
       function buildCOV(obj)
                obj.COVMOD = getMOD(obj,obj.DepVar,obj.COV,1:obj.nrCOV,[],1);
       end
       function crossValidate(obj,K,t,ExtDepVar,ExtCOV,ExtGB)
                if isscalar(K)
                   F = DOB_SCV(K,obj.DepVar,obj.COV(:,1));
                else
                   F = K;K = length(unique(F));
                end
                Fold = cell(1,K);
                Res = cell(1,K);
                allind = 1:obj.n;
                parfor k=1:1:K
                       % k= 1;
                       TestInd = find(F==k);
                       TestCOV = obj.COV(TestInd,:);
                       TestGB = obj.GB(TestInd,:);
                       TestDepVar = obj.DepVar(TestInd,:);
                       TrInd = setdiff(allind,TestInd); %#ok<*PROP>
                       Fold{k} = clone(obj);
                       reduceSamples(Fold{k},TrInd);
                       Fold{k}.MOD = getMOD(Fold{k},Fold{k}.DepVar,[Fold{k}.COV, Fold{k}.GB],1:Fold{k}.nrCOV+Fold{k}.nrGB,[],1);
                       RIP = BaseContainerv5.updateRIP(TestDepVar,Fold{k}.MOD.M)';
                       if Fold{k}.RIPNormalize, RIP = BaseContainerv5.normalizeRIP(RIP,Fold{k}.MOD.M,Fold{k}.MOD.Var);end
                       [Res{k}.Stat,~,Res{k}.pStat] = performRIPStat([TestCOV, TestGB],RIP,Fold{k}.MOD.Var,t);
                       Res{k}.Acc = getBiometricAcc(Fold{k},TestDepVar,TestCOV,TestGB);
                end
                obj.CV.Stat = nan*zeros(obj.nrCOVGB,K);
                obj.CV.pStat = nan*zeros(obj.nrCOVGB,K);
                obj.CV.Acc = nan*zeros(obj.nrCOVGB,K);
                for k=1:1:K
                    obj.CV.Stat(:,k) = Res{k}.Stat';
                    obj.CV.pStat(:,k) = Res{k}.pStat';
                    obj.CV.Acc(:,k) = Res{k}.Acc';
                end
                obj.CV.pfast = nan*zeros(1,obj.nrCOVGB);
                for i=1:1:obj.nrCOVGB
                    p = obj.CV.pStat(i,:);
                    p(find(p==0)) = 0.1/t;
                    p = p(find(~isnan(obj.CV.Stat(i,:))));
                    if ~isempty(p), obj.CV.pfast(i) = BaseContainerv5.pfast(p); end
                end
                if nargin < 4, return; end
                Res = cell(1,K);
                parfor k=1:K
                    Res{k} = extTest(Fold{k},ExtDepVar,ExtCOV,ExtGB)
                    Res{k}.Acc = getBiometricAcc(Fold{k},ExtDepVar,ExtCOV,ExtGB);
                end
                obj.EV.Stat = nan*zeros(obj.nrCOVGB,K);
                obj.EV.pStat = nan*zeros(obj.nrCOVGB,K);
                obj.EV.Acc = nan*zeros(obj.nrCOVGB,K);
                for k=1:1:K
                    obj.EV.Stat(:,k) = Res{k}.Stat';
                    obj.EV.pStat(:,k) = Res{k}.pStat';
                    obj.EV.Acc(:,k) = Res{k}.Acc';
                end
                obj.EV.pfast = nan*zeros(1,obj.nrCOVGB);
                for i=1:1:obj.nrCOVGB
                    p = obj.EV.pStat(i,:);
                    p(find(p==0)) = 0.1/t;
                    p = p(find(~isnan(obj.EV.Stat(i,:))));
                    if ~isempty(p), obj.EV.pfast(i) = BaseContainerv5.pfast(p); end
                end
       end
       function build(obj)
                obj.MOD = getMOD(obj,obj.DepVar,[obj.COV, obj.GB],1:obj.nrCOVGB,[],1);
       end
       function out = getMOD(obj,DepVar,IndVar,BootInd,CondInd,PartInd)
           switch obj.Method
               case 'PLSR'
                   out = PLSR(obj,DepVar,IndVar,BootInd,CondInd);
               case 'BRIM'
                   out = BRIM(obj,DepVar,IndVar,BootInd,CondInd,PartInd);
               case 'FOLDEDBRIM'
                   out = FOLDEDBRIM(obj,DepVar,IndVar,BootInd,CondInd,PartInd);
           end
       end
       function out = FOLDEDBRIM(obj,DepVar,IndVar,BootInd,CondInd,PartInd)
                % Initialize
                  Bootstrap = true;if obj.MaxIterations == 0,Bootstrap = false;end
                  [n,~] = size(DepVar);Ind = (1:n);                    
                % Building Independent and Dependent Variables
                  nr2Condition = length(CondInd);
                  var.nr2Boot = length(BootInd);
                  var.Booting = (nr2Condition+1:nr2Condition+var.nr2Boot);
                % Partitioning Variable used to define the Folds
                  G = ones(n,1);if ~isempty(PartInd), G = IndVar(:,PartInd);end
                % Output
                  out.IndVar = IndVar;out.DepVar = DepVar;
                % examining variables to Boot and create reference faces, the lattter is
                % needed to rescale RIP variables onto common scales
                  for i=1:1:var.nr2Boot
                     v = IndVar(:,var.Booting(i));
                     var.Info{i} = BaseContainerv5.getVarInfo(v);    
                  end
                  out.Var = var;
                % Outer partitioning of data
                  Fouter = DOB_SCV(obj.OuterFold,DepVar,G); 
                  FoldResults = cell(1,obj.OuterFold);
                % avoid passing complete object to parfor  
                  RegRuns = obj.RegRuns;MaxIterations = obj.MaxIterations;Htest = obj.Htest;RIPNormalize = obj.RIPNormalize;
                  StopCorr = obj.StopCorr; RIPStatPerm = obj.RIPStatPerm;OuterFold = obj.OuterFold;InnerFold = obj.InnerFold;
                  [path,ID] = setupParForProgress(OuterFold); 
                  parfor fo=1:OuterFold
                        FoTestInd = find(Fouter==fo);
                        FoldResults{fo}.TestInd = FoTestInd;
                        FoTrInd = setdiff(Ind,FoTestInd);
                        FoldResults{fo}.TrInd = FoTrInd;
                        FoTrIndVar = IndVar(FoTrInd,:); %#ok<*PFBNS>
                        FoTrDepVar = DepVar(FoTrInd,:);
                        FoTestIndVar = IndVar(FoTestInd,:);
                        FoTestDepVar = DepVar(FoTestInd,:);
                        FoldResults{fo}.DepVar = FoTestDepVar;
                        FoTrOrigIndVar = FoTrIndVar;
                        FoTrNewIndVar = FoTrIndVar;
                        BootProgress = zeros(var.nr2Boot,MaxIterations,2);
                        ContBoot = true;counter = 0;
                        while ContBoot && Bootstrap>0
                           % seperate data into inner folds
                           Finner = DOB_SCV(InnerFold,FoTrDepVar,G(FoTrInd));
                           counter = counter + 1;
                           for fi=1:InnerFold
                                FiTestInd = find(Finner==fi);
                                FiTrInd = setdiff(1:length(FoTrInd),FiTestInd);
                                [M,H] = BaseContainerv5.getBootSampleRegression(FoTrIndVar(FiTrInd,:),FoTrDepVar(FiTrInd,:),var,RegRuns,obj.HtestP);
                                if Htest, M = H.*M;end% setting close to zero partial regression coefficients to zero
                                rip = BaseContainerv5.updateRIP(FoTrDepVar(FiTestInd,:),M)';
                                if RIPNormalize, rip = BaseContainerv5.normalizeRIP(rip,M,var);end
                                FoTrNewIndVar(FiTestInd,var.Booting) = rip;
                           end
                           % monitoring progress of Booting
                           ToStop = zeros(1,var.nr2Boot);
                           for b=1:1:var.nr2Boot
                               V1 = FoTrIndVar(:,var.Booting(b));
                               V2 = FoTrNewIndVar(:,var.Booting(b));
                               index =  ~isnan(V1);
                               tmp = corrcoef(V1(index),V2(index));
                               BootProgress(b,counter,2) = tmp(1,2);
                               if tmp(1,2)>=StopCorr, ToStop(b) = 1;end
                               V1 = FoTrOrigIndVar(:,var.Booting(b));
                               index =  ~isnan(V1);
                               tmp = corrcoef(V1(index),V2(index));
                               BootProgress(b,counter,1) = tmp(1,2);
                           end
                           if (sum(ToStop)/var.nr2Boot==1), ContBoot = false; end
                           if counter >= MaxIterations, ContBoot = false; end
                           FoTrIndVar = FoTrNewIndVar;
                        end
                        if counter>0, FoldResults{fo}.BootProgress = BootProgress(:,1:counter,:);end
                        FoldResults{fo}.BootIter = counter;
                        % getting the final M from outer training data
                        [M,H] = BaseContainerv5.getBootSampleRegression(FoTrIndVar,FoTrDepVar,var,RegRuns,obj.HtestP);
                        if Htest, M = H.*M;end
                        FoldResults{fo}.M = M;FoldResults{fo}.H = H;
                        % estimating rip values outer test data
                        rip = BaseContainerv5.updateRIP(FoTestDepVar,M)';
                        if RIPNormalize, rip = BaseContainerv5.normalizeRIP(rip,M,var);end
                        FoldResults{fo}.RIP = rip;
                        FoldResults{fo}.IndVar = FoTestIndVar(:,var.Booting);
                        % perform within outer fold statistics
                        [FoldResults{fo}.stat,FoldResults{fo}.statType,FoldResults{fo}.statP] = BaseContainerv5.performRIPStat(FoldResults{fo}.IndVar,rip,var,RIPStatPerm);
                        parfor_progress;
                  end
                  closeParForProgress(path,ID);
                % crossfold results
                  Stat = zeros(length(var.Booting),OuterFold);
                  StatP = zeros(length(var.Booting),OuterFold);
                  StatType = cell(length(var.Booting),OuterFold);
                  M = zeros(var.nr2Boot,size(DepVar,2),OuterFold);
                  for fo=1:OuterFold
                      Stat(:,fo) = FoldResults{fo}.stat;StatType(:,fo) = FoldResults{fo}.statType;
                      StatP(:,fo) = FoldResults{fo}.statP;M(:,:,fo) = FoldResults{fo}.M; 
                  end
                  out.Stat = Stat;out.StatType = StatType;out.StatP = StatP;
                  out.StatP(out.StatP==0) = 0.1/RIPStatPerm;
                  out.pfast = ones(1,var.nr2Boot);
                  for i=1:1:var.nr2Boot
                      out.pfast(i) = BaseContainerv5.pfast(out.StatP(i,:));
                  end
                % Extracting final M 
                  out.M = median(M,3);
                  out.Var = var;
                  if ~Htest, return; end
                  if OuterFold>=8, 
                     out.H = BaseContainerv5.wilcoxonM(M,var,obj.HtestP);
                  else
                     out.H = ones(size(out.M)); 
                  end
                  out.M = out.H.*out.M;
                  RIP = BaseContainerv5.updateRIP(DepVar,M)';
                  if obj.RIPNormalize, RIP = BaseContainerv5.normalizeRIP(RIP,M,var);end
                  out.RIP = RIP;
                  [out.RIPStat] = BaseContainerv5.performRIPStat(out.IndVar,RIP,var,0);
       end
       function out = BRIM(obj,DepVar,IndVar,BootInd,CondInd,PartInd)
                % Initialize
                  Bootstrap = true;if obj.MaxIterations == 0,Bootstrap = false;end
                  [n,~] = size(DepVar);Ind = (1:n);                    
                % Building Independent and Dependent Variables
                  nr2Condition = length(CondInd);
                  var.nr2Boot = length(BootInd);
                  var.Booting = (nr2Condition+1:nr2Condition+var.nr2Boot);
                % Partitioning Variable used to define the Folds
                  G = ones(n,1);if ~isempty(PartInd), G = IndVar(:,PartInd);end
                % Output
                  out.IndVar = IndVar;out.DepVar = DepVar;
                % examining variables to Boot and create reference faces, the lattter is
                % needed to rescale RIP variables onto common scales
                  for i=1:1:var.nr2Boot
                     v = IndVar(:,var.Booting(i));
                     var.Info{i} = BaseContainerv5.getVarInfo(v,DepVar);    
                  end
                  out.Var = var;
                  OrigIndVar = IndVar;
                  NewIndVar = IndVar;
                  BootProgress = zeros(var.nr2Boot,obj.MaxIterations,2);
                  ContBoot = true;counter = 0;
                  while ContBoot && Bootstrap>0
                           % seperate data into inner folds
                           Finner = DOB_SCV(obj.InnerFold,DepVar,G);
                           counter = counter + 1;
                           for fi=1:obj.InnerFold
                                FiTestInd = find(Finner==fi);
                                FiTrInd = setdiff(1:length(Ind),FiTestInd);
                                [M,H] = BaseContainerv5.getBootSampleRegression(IndVar(FiTrInd,:),DepVar(FiTrInd,:),var,obj.RegRuns,obj.SamplePercentage,obj.HtestP);
                                if obj.Htest, M = H.*M;end% setting close to zero partial regression coefficients to zero
                                rip = BaseContainerv5.updateRIP(DepVar(FiTestInd,:),M)';
                                if obj.RIPNormalize, rip = BaseContainerv5.normalizeRIP(rip,M,var);end
                                NewIndVar(FiTestInd,var.Booting) = rip;
                           end
                           % monitoring progress of Booting
                           ToStop = zeros(1,var.nr2Boot);
                           for b=1:1:var.nr2Boot
                               V1 = IndVar(:,var.Booting(b));
                               V2 = NewIndVar(:,var.Booting(b));
                               index =  ~isnan(V1);
                               tmp = corrcoef(V1(index),V2(index));
                               BootProgress(b,counter,2) = tmp(1,2);
                               if tmp(1,2)>=obj.StopCorr, ToStop(b) = 1;end
                               V1 = OrigIndVar(:,var.Booting(b));
                               index =  ~isnan(V1);
                               tmp = corrcoef(V1(index),V2(index));
                               BootProgress(b,counter,1) = tmp(1,2);
                           end
                           if (sum(ToStop)/var.nr2Boot==1), ContBoot = false; end
                           if counter >= obj.MaxIterations, ContBoot = false; end
                           IndVar = NewIndVar;
                  end
                  % getting the final M from outer training data
                  [M,H] = BaseContainerv5.getBootSampleRegression(IndVar,DepVar,var,obj.RegRuns,obj.SamplePercentage,obj.HtestP);
                  if obj.Htest, M = H.*M;end
                  out.M = M;out.H = H;
                  RIP = BaseContainerv5.updateRIP(DepVar,M)';
                  if obj.RIPNormalize, RIP = BaseContainerv5.normalizeRIP(RIP,M,var);end
                  out.RIP = RIP;
                  [out.RIPStat] = BaseContainerv5.performRIPStat(out.IndVar,RIP,var,0);
       end
       function out = PLSR(obj,DepVar,IndVar,BootInd,CondInd)                  
                % Building Independent and Dependent Variables
                  nr2Condition = length(CondInd);
                  var.nr2Boot = length(BootInd);
                  var.Booting = (nr2Condition+1:nr2Condition+var.nr2Boot);
                % Output
                  out.IndVar = IndVar;out.DepVar = DepVar;
                % examining variables to Boot and create reference faces, the lattter is
                % needed to rescale RIP variables onto common scales
                  for i=1:1:var.nr2Boot
                     v = IndVar(:,var.Booting(i));
                     var.Info{i} = BaseContainerv5.getVarInfo(v,DepVar);    
                  end
                  out.Var = var;
                  % getting the final M from outer training data
                  [M,H] = BaseContainerv5.getBootSampleRegression(IndVar,DepVar,var,obj.RegRuns,obj.SamplePercentage,obj.HtestP);
                  if obj.Htest, M = H.*M;end
                  out.M = M;out.H = H;
                  RIP = BaseContainerv5.updateRIP(DepVar,M)';
                  if obj.RIPNormalize, RIP = BaseContainerv5.normalizeRIP(RIP,M,var);end
                  out.RIP = RIP;
                  [out.RIPStat] = BaseContainerv5.performRIPStat(out.IndVar,RIP,var,0);
       end
       function obj = reduceSamples(obj,index)
          if nargout==1, obj = clone(obj); end
          obj.ShapeSpace.Tcoeff = obj.ShapeSpace.Tcoeff(index,:);
          obj.COV = obj.COV(index,:);
          obj.GB = obj.GB(index,:);
          if ~isempty(obj.RedShapeSpace),obj.RedShapeSpace.Tcoeff = obj.RedShapeSpace.Tcoeff(index,:);end
          %if~isempty(obj.MOD), obj.MOD.RIP = obj.MOD.RIP(index,:); end
          obj.COVMOD = [];
          obj.GBMOD = [];
       end
       function createReducedSpace(obj)
           if isempty(obj.ShapeSpace), return; end
           switch obj.RedShapeSpaceType
               case 'X'
                   if isempty(obj.COV), return; end
                   if isempty(obj.GB), return; end
                   A = [obj.COV, obj.GB];
               case 'RIP'
                   if isempty(obj.COVRIP), return; end
                   if isempty(obj.GBRIP), return; end
                   A = [obj.COVRIP, obj.GBRIP];
               otherwise
                   return;
           end
           % setting nanvalues to average
           for i=1:1:size(A,2)
               A(isnan(A(:,i)),i) = nanmean(A(:,i));
           end
           % Building regression model
           rm = PLSRShapeModel;rm.Model = clone(obj.ShapeSpace);rm.X = A;update(rm);
           % changing faces
           RedDepVar = obj.ShapeSpace.Tcoeff;avgA = mean(A);
           for i=1:1:obj.n
               RedDepVar(i,:) = changeX(rm,avgA,A(i,:),obj.ShapeSpace.Tcoeff(i,:));
           end
           tmpSpace = clone(obj.ShapeSpace);
           tmpSpace.Tcoeff = RedDepVar;
           shape = reconstructTraining(tmpSpace);
           RedShapeSpace = shapePCA;
           RedShapeSpace.RefScan = clone(obj.ShapeSpace.Average);
           getAverage(RedShapeSpace,shape);
           getModel(RedShapeSpace,shape);
           reduceNrPc(RedShapeSpace,obj.ShapeSpace.nrEV-size(A,2));
           obj.RedShapeSpace = RedShapeSpace;
       end
       function reduceGB(obj,index)
           obj.GB = obj.GB(:,index);
           obj.GBindex = obj.GBindex(index);
       end
   end  
   methods % MODEL FUNCTION
       function [Stat,pStat] = COVRIPStatistics(obj,testdepvar,testcov,t)
           % a cross test to see if the rip values in an external test set
           % are correlated with the original cov values in the test set
                RIP = BaseContainerv5.updateRIP(testdepvar,obj.COVMOD.M)';
                if obj.RIPNormalize, RIP = BaseContainerv5.normalizeRIP(RIP,obj.COVMOD.M,obj.COVMOD.Var);end
                [Stat,~,pStat] = BaseContainerv5.performRIPStat(testcov,RIP,obj.COVMOD.Var,t);
       end
       function [EER,R] = COVSubBiometrics(obj,testdepvar,testcov)
                RIP = BaseContainerv5.updateRIP(testdepvar,obj.COVMOD.M)';
                if obj.RIPNormalize, RIP = BaseContainerv5.normalizeRIP(RIP,obj.COVMOD.M,obj.COVMOD.Var);end
                EER = zeros(1,obj.nrCOV);
                R = zeros(obj.nrCOV,4);
                for i=1:1:obj.nrCOV
                    [EER(i),R(i,:)] = subBiometrics(obj,obj.COVMOD.Var.Info{i},obj.COV(:,i),obj.COVRIP(:,i),testcov(:,i),RIP(:,i));
                    
                end
       end
       function [EER,R,X,Y] = subBiometrics(obj,varinfo,trX,trrip,testX,testrip)
           % in subBiometrics a single variable is investigated where data
           % is seperated in opposite groups to form true and false matches
           % defining opposite groups
           switch varinfo.Type
               case 'Categorical'
                   posclass = -1;% males
                   Tind = find(testX==posclass);
                   Find = find(testX==-1*posclass);
               case 'Continuous'
                   Tind = find(testX<=nanmean(trX));
                   Find = find(testX>nanmean(trX));
           end
           nrT = length(Tind);TMatches = nan*zeros(1,nrT);
           nrF = length(Find);FMatches = nan*zeros(1,nrT,nrF);
           for t=1:1:nrT
               distr = BaseContainerv5.X2RIP(trX,trrip,testX(Tind(t)),varinfo.Type,6,200);
               TMatches(1,t) = getMatch(obj,distr,testrip(Tind(t)));
               FMatches(1,t,:) = getMatch(obj,distr,testrip(Find));
           end
           TScores = match2score(obj,TMatches);
           FScores = match2score(obj,FMatches);
           [R,EER,~,~,X,Y] = BaseContainerv5.getVerificationIdentification(TScores,FScores,ones(1,nrT),1);
       end
       function out = getMatch(obj,distr,rip)
                pT = normpdf(rip,distr.Tmu,distr.Tsigma)';
                pF = normpdf(rip,distr.Fmu,distr.Fsigma)';
                switch obj.Match
                    case 'QD'
                        out = pT./pF;
                    case 'INLIER'
                        out = pT./(pT+pF);
                end  
       end
       function out = match2class(obj,matches)
          switch obj.Match
              case 'QD'
                  T = 1;
              case 'INLIER'
                  T = 0.5;
          end
          out = BaseContainerv5.matchthresholding(matches,T);
       end
       function out = match2score(obj,matches)
          switch obj.Match
              case 'QD'
                  out = -log(matches);
                  %out = 1./matches;
              case 'INLIER'
                  out = 1-matches;
          end
       end
   end
   methods (Static = true)
       function E = getResiduals(A,B)
                [~,~,~,~,~,~,~,stats] = plsregress(A,B,min(size(A,2),size(B,2)));
                E = stats.Yresiduals;
       end
       function Stat = getLocalStatistic(R,B)
                n = size(B,1);P = B-R;A = repmat(mean(B),n,1);
                Stat = sum((P-A).^2)./sum((B-A).^2);
       end
       function M = getRegression(A,B,Booting)
                 if nargin < 3, Booting = (1:size(A,2));end
                 [A,B] = eliminateNAN(A,B);
                 [~,~,~,~,M] = plsregress(A,B,size(A,2));
                 M = M(1+Booting,:);
        end
       function [out] = updateRIP(in,M)
                    % in this implementation I take the reference as the origin
                    n2 = size(in,1); n1 = size(M,1);% determine input size
                    in = in'; out = nan*zeros(n1,n2);% allocate memory
                    for j=1:n1
                        out(j,:) = dot(in,repmat(M(j,:)'/norm(M(j,:)'),1,n2));
                    end
       end
       function [stat,statType,statP] = performRIPStat(vIn,vripIn,var,t)
                stat = zeros(1,length(var.Booting));
                statType = cell(1,length(var.Booting));
                statP = zeros(1,length(var.Booting));
                for s=1:1:length(var.Booting)
                   v = vIn(:,s);
                   vrip = vripIn(:,s);
                   index = find(~isnan(v));
                   switch var.Info{s}.Type
                       case 'Categorical'
                           statType{s} = 'F';
                           XG = cell(size(v(index)));
                           for l=1:length(var.Info{s}.El)
                               XG((v(index)==var.Info{s}.El(l))) = {num2str(l)};
                           end
                           [stat(s),~,statP(s)] = BaseContainerv5.myPermAnova(XG,vrip(index),t);
                       case 'Continuous'
                           statType{s} = 'Cor';
                           [stat(s),statP(s)] = BaseContainerv5.permCorr(v(index),vrip(index),t);
                   end
                end  
       end
       function [out] = normalizeRIP(rip,M,var)
             out = rip;
             for b= 1:1:var.nr2Boot % testing reference faces and rescale accordingly
                Mrip = updateRIP(var.Info{b}.MDepVar,M)';
                Prip = updateRIP(var.Info{b}.PDepVar,M)';
                tmp = rip(:,b);
                range = Prip(b)-Mrip(b);
                tmp = (tmp-Mrip(b))/range;           
                out(:,b) = tmp*var.Info{b}.Range+var.Info{b}.Mel;  
             end
       end
       function [M,H] = getBootSampleRegression(IndVar,DepVar,var,runs,perc,Tp)
                if nargin < 5, perc = 1; end
                nrS = size(IndVar,1);
                M = zeros(var.nr2Boot,size(DepVar,2),runs);
                for s=1:1:runs
                    if runs==1% no random sampling
                       SampleFi = 1:nrS; 
                    else % random sampling with replacement
                       SampleFi = randsample(nrS,round(perc*nrS),true);
                    end
                    M(:,:,s) = BaseContainerv5.getRegression(IndVar(SampleFi,:),DepVar(SampleFi,:),var.Booting);
                end
                if runs==1,M = squeeze(M);H = ones(size(M));return;end
                H = BaseContainerv5.wilcoxonM(M,var,Tp);
                M = median(M,3);       
       end
       function H = wilcoxonM(M,var,Tp)
            P = zeros(var.nr2Boot,size(M,2));
            for i=1:1:var.nr2Boot
               Mfor = squeeze(M(i,:,:)); 
               for j=1:1:size(Mfor,1)
                  P(i,j) = signrank(Mfor(j,:));    
               end
            end
            H = P<=Tp;
       end
       function p = pfast(p)
                product=prod(p);
                n=length(p);
                if n<=0
                   error('pfast was passed an empty array of p-values')
                elseif n==1
                   p = product;
                   return;
                elseif product == 0
                   p = 0;
                   return
                else
                   x = -log(product);
                   t=product;
                   p=product;
                   for i = 1:n-1
                       t = t * x / i;
                       p = p + t;
                   end
                end  
       end
       function [AVG,STD,O] = extractRobustRIPStatistics(rip,kappa,W2)
                nrO = length(rip);
                if nargin<3, W2 = ones(1,nrO); end 
                O.avg = sum(W2.*rip)/sum(W2);
                O.std = sqrt(sum(W2.*((rip-O.avg).^2))/sum(W2));
                W1 = ones(1,nrO);  
                for i=1:1:10
                    W = W1.*W2;
                    AVG = sum(W.*rip)/sum(W); 
                    STD = sqrt(sum(W.*((rip-AVG).^2))/sum(W));
                    IP = normpdf(rip,AVG,STD);
                    %L = (1/sqrt(((2*pi)^2)*STD))*exp(-0.5*kappa^2);
                    L = (1/(sqrt(2*pi)*STD))*exp(-0.5*kappa^2);
                    W1 = IP./(IP+L);
                end
       end
       function out = getVarInfo(v,DepVar)
                     el = unique(v(~isnan(v)));
                     if length(el)<6%Categorical parameter, currently allowing up to 6 categories  
                        out.Type = 'Categorical';
                        out.El = el;out.Mel = min(el);out.Pel = max(el);
                        out.Range = out.Pel-out.Mel;
                        % getting categorical reference faces
                        out.Mindex = find(v==out.Mel);
                        out.MDepVar = mean(DepVar(out.Mindex,:)); %#ok<*FNDSB>
                        out.Pindex = find(v==out.Pel);
                        out.PDepVar = mean(DepVar(out.Pindex,:));
                     else% variable is Continuous
                        out.Type = 'Continuous';
                        index = find(~isnan(v));
                        [sortv,ind] = sort(v(index),'ascend');
                        % lower 25%
                        indl = round(0.25*length(ind));
                        out.Mb = sortv(indl);
                        out.Mindex = find(v<=out.Mb);
                        out.Mel = mean(v(out.Mindex));
                        out.Mrange = [min(v(out.Mindex)) max(v(out.Mindex))];
                        out.MDepVar = mean(DepVar(out.Mindex,:));
                        % upper 25%
                        indu = round(0.75*length(ind));
                        out.Pb = sortv(indu);
                        out.Pindex = find(v>=out.Pb);
                        out.Pel = mean(v(out.Pindex));
                        out.Prange = [min(v(out.Pindex)) max(v(out.Pindex))];
                        out.PDepVar = mean(DepVar(out.Pindex,:));
                        out.Range = out.Pel-out.Mel;
                     end
       end
       function [FA,pFA,ppermFA,T] = myPermAnova(XG,Y,t)  
         [~,T] = anova1(Y,XG,'off');
         FA = T{2,5};
         pFA  = T{2,6};
         % permuting
         if t==0, ppermFA = nan; return;  end
         n = length(XG);
         FAnull = zeros(1,t);
         parfor i=1:t
             ind = randperm(n);
             Yfor = Y;
             Yfor = Yfor(ind,:);
             [~,Tfor] = anova1(Yfor,XG,'off');
             FAnull(i) = Tfor{2,5};       
         end
         ppermFA = sum(FAnull>=FA)/t;
       end
       function [c,pc]=permCorr(x,y,t)
                if nargin < 3, t = 1000; end
                n=length(x);
                corrobs=corrcoef(x,y);
                c = corrobs(2,1);
                if t==0, pc = nan; return; end
                corrCount=zeros(1,t);
                parfor i=1:t,
                    ind=randperm(n);
                    xfor=x(ind);
                    r=corrcoef(xfor,y);
                    cfor=r(2,1);
                    switch sign(c)
                        case 1
                            corrCount(i) = cfor>=c;
                        case -1
                            corrCount(i) = cfor<=c;
                    end
                end;
                pc = sum(corrCount)/t;
       end
       function out = matchthresholding(matches,T)
                out = matches>=T; 
       end
       function out = class2score(classes)
                out = 1-classes; 
       end   
       function [R,EER,ST,SF,X,Y] = getVerificationIdentification(TScores,FScores,wM,W)
                 nF = size(FScores,3);
                 nT = size(TScores,2);
                 tmpW = repmat(W,1,nT).*wM;
                 ST = BaseContainerv5.accScores(TScores,tmpW)';
                 SF = BaseContainerv5.accScores(FScores,repmat(tmpW,1,1,nF));
                 [R,EER,X,Y] = BaseContainerv5.evalScores(ST,SF,nT,nF); 
       end
       function Score = accScores(Scores,tmpW)
                Score = squeeze(nansum(tmpW.*Scores,1)./nansum(tmpW,1));
       end
       function [R,EER,x,y] = evalScores(ST,SF,nT,nF)
                 % Rank analysis
                 R = (sum(SF<repmat(ST,1,nF),2)+1)/(nF+1);
                 RM = median(R);
                 R1 = sum(R<=0.01)/length(R);
                 R10 = sum(R<=0.10)/length(R);
                 R20 = sum(R<=0.20)/length(R);
                 R = [R1 R10 R20 RM]; 
                 % EER analysis
                 g = [ones(1,nT), zeros(1,nT*nF)];
                 [~,order] = sort([ST;SF(:)]);
                 g = g(order);
                 true_neg = g == 0;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
                 true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
                 y = tpf(1:end-1)+dy./2;
                 x = fpf(1:end-1)+dx./2;
                 yn = 1-y;d = abs(x-yn);
                 [~,indmin] = min(d);
                 EER = ((x(indmin)+yn(indmin))/2);  
       end
       function [Distr,RIP] = X2RIP(trX,trrip,X,type,kappa,K)
           if nargin<5, kappa = 6; end
           if nargin<6, K = 300; end
           switch type
               case 'Categorical'
                  Tind = find(trX==X);% Inlier Distribution
                  Find = find(trX==-1*X);% Outlier Distribution
                  [Distr.Tmu,Distr.Tsigma,Distr.TO] = BaseContainerv5.extractRobustRIPStatistics(trrip(Tind)',kappa);
                  RIP = Distr.Tmu;Distr.TpMax = normpdf(Distr.Tmu,Distr.Tmu,Distr.Tsigma);
                  [Distr.Fmu,Distr.Fsigma,Distr.FO] = BaseContainerv5.extractRobustRIPStatistics(trrip(Find)',kappa);
                  Distr.FpMax = normpdf(Distr.Fmu,Distr.Fmu,Distr.Fsigma);      
               case 'Continuous'
                  ind = find(trX>0);% eliminate nans and wrong data entries (typically coded as 0);
                  trX = trX(ind);
                  W = abs(trX-X);
                  % select closest samples
                  [~,ind2] = sort(W,'ascend');
                  % TRUE DISTRIBUTION
                  ind2T = ind2(1:K);indT = ind(ind2T);
                  rip = trrip(indT)';
                  Distr.IndT = indT;
                  WT = W(ind2T);WT = WT-min(WT);WT = WT/max(WT);WT = 1-WT';
                  [Distr.Tmu,Distr.Tsigma,Distr.TO] = BaseContainerv5.extractRobustRIPStatistics(rip,kappa,WT);
                  Distr.TpMax = normpdf(Distr.Tmu,Distr.Tmu,Distr.Tsigma);
                  RIP = Distr.Tmu;
                  % FALSE DISTRIBUTION
                  ind2F = ind2(end-K+1:end);indF = ind(ind2F);
                  rip = trrip(indF)';
                  Distr.IndF = indF;
                  WF = W(ind2F);WF = WF-min(WF);WF = WF/max(WF);WF = WF';
                  [Distr.Fmu,Distr.Fsigma,Distr.FO] = BaseContainerv5.extractRobustRIPStatistics(rip,kappa,WF);
                  Distr.FpMax = normpdf(Distr.Fmu,Distr.Fmu,Distr.Fsigma); 
           end
       end
   end
end



% function out = extTest(obj,ExtDepVar,ExtCOV,ExtGB,t)
%            RIP = BaseContainerv5.updateRIP(ExtDepVar,obj.MOD.M)';
%            if obj.RIPNormalize, RIP = BaseContainerv5.normalizeRIP(RIP,obj.MOD.M,obj.MOD.Var);end
%            [out.Stat,~,out.pStat] = performRIPStat([ExtCOV, ExtGB],RIP,obj.MOD.Var,t);
% end


% function out = biometrics(obj,tDepVar,tCOV,tGB,rDepVar,vis)
%            if isscalar(rDepVar), rDepVar = randomFaces(obj,rDepVar);end
%            nT = size(tDepVar,1);nR = size(rDepVar,1);
%            tRIP = updateRIP(tDepVar,obj.MOD.M)';if obj.RIPNormalize, tRIP = BaseContainerv5.normalizeRIP(tRIP,obj.MOD.M,obj.MOD.Var);end
%            rRIP = updateRIP(rDepVar,obj.MOD.M)';if obj.RIPNormalize, rRIP = BaseContainerv5.normalizeRIP(rRIP,obj.MOD.M,obj.MOD.Var);end
%            tCOVMatches = nan*zeros(obj.nrCOV,nT,nT-1);rCOVMatches = nan*zeros(obj.nrCOV,nT,nR);COVMatches = nan*zeros(obj.nrCOV,nT);
%            tGBMatches = nan*zeros(obj.nrGB,nT,nT-1);rGBMatches = nan*zeros(obj.nrGB,nT,nR);GBMatches = nan*zeros(obj.nrGB,nT);
%            parfor t=1:nT
%                 tInd = setdiff(1:nT,t);
%                 [~,COVD] = COVX2RIP(obj,tCOV(t,:),6,300);
%                 [~,gbD] = GBX2RIP(obj,tGB(t,:),6,300);
%                 tmp = matchCOVRIP(obj,tRIP(:,1:obj.nrCOV),COVD);
%                 tCOVMatches(:,t,:) = tmp(:,tInd);COVMatches(:,t) = tmp(:,t);
%                 rCOVMatches(:,t,:) = matchCOVRIP(obj,rRIP(:,1:obj.nrCOV),COVD);
%                 tmp = matchGBRIP(obj,tRIP(:,obj.nrCOV+1:end),gbD);
%                 tGBMatches(:,t,:) = tmp(:,tInd);GBMatches(:,t) = tmp(:,t);
%                 rGBMatches(:,t,:) = matchGBRIP(obj,rRIP(:,obj.nrCOV+1:end),gbD);
%            end
%            Matches = [COVMatches;GBMatches];out.Matches = Matches;
%            tMatches = [tCOVMatches;tGBMatches];out.tMatches = tMatches;
%            rMatches = [rCOVMatches;rGBMatches];out.rMatchtes = rMatches;
%            switch obj.Classify
%                case 'SOFT'
%                    Scores = match2score(obj,Matches);out.Scores = Scores;
%                    tScores = match2score(obj,tMatches);out.tScores = tScores;
%                    rScores = match2score(obj,rMatches);out.rScores = rScores;
%                case 'HARD'
%                    Scores = BaseContainerv5.class2score(match2class(obj,Matches));out.Scores = Scores;
%                    tScores = BaseContainerv5.class2score(match2class(obj,tMatches));out.tScores = tScores;
%                    rScores = BaseContainerv5.class2score(match2class(obj,rMatches));out.rScores = rScores;
%            end    
%            out.ClassAcc = ((sum(match2class(obj,Matches),2)/nT)*100)';
%            nr = obj.nrCOVGB+4;
%            ind = cell(1,nr);
%            for i=1:1:obj.nrCOVGB % all variables individually
%                ind{i} = i; 
%            end
%            ind{obj.nrCOVGB+1} = 1:obj.nrCOV;% COVariates only
%            ind{obj.nrCOVGB+2} = obj.nrCOV+1:obj.nrCOVGB;% genetic background only
%            ind{obj.nrCOVGB+3} = 1:obj.nrCOVGB;% all
%            ind{obj.nrCOVGB+4} = [1 obj.nrCOV+1:obj.nrCOVGB];% sex and genetic background
%            ST = nan*zeros(nr,nT);SFr = nan*zeros(nr,nT,nR);SFt = nan*zeros(nr,nT,nT-1);
%            EER = nan*zeros(2,nr);Ranks = nan*zeros(2,nr,nT);
%            rX = nan*zeros(nr,(nT*(nR+1))-1);rY = nan*zeros(nr,(nT*(nR+1))-1);
%            tX = nan*zeros(nr,(nT*nT)-1);tY = nan*zeros(nr,(nT*nT)-1);
%            for r=1:1:nr
%                 [Ranks(1,r,:),EER(1,r),ST(r,:),SFr(r,:,:),rX(r,:),rY(r,:)] = BaseContainerv5.getBiometrics(Scores(ind{r},:),rScores(ind{r},:,:),ones(length(ind{r}),nT),ones(size(ones(length(ind{r}),nT),1),1));
%                 [Ranks(2,r,:),EER(2,r),~,SFt(r,:,:),tX(r,:),tY(r,:)] = BaseContainerv5.getBiometrics(Scores(ind{r},:),tScores(ind{r},:,:),ones(length(ind{r}),nT),ones(size(ones(length(ind{r}),nT),1),1));        
%            end
%            Ranks(1,:,:) = (Ranks(1,:,:)./(nR+1)).*100;
%            Ranks(2,:,:) = (Ranks(2,:,:)./nT).*100;
%            CumRanks = zeros(2,nr,100);
%            for cr=1:1:100
%                CumRanks(:,:,cr) = (sum(Ranks<=cr,3)/nT).*100;
%            end
%            out.EER = EER;
%            out.CumRanks = CumRanks;
%            if vis
%                str = {'b-' 'b-.' 'b:' 'b--' 'r-' 'r-.' 'r:' 'k-' 'm-' 'g-' 'c-' 'y-'};
%                figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
%                xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
%                title('Identification RANDOM');
%                plot(1:100,1:100,'k--');
%                for i=1:1:nr
%                    plot(1:100,squeeze(CumRanks(1,i,:))',str{i},'LineWidth',1);
%                end
%                figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
%                xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
%                title('Identification TEST');
%                plot(1:100,1:100,'k--');
%                for i=1:1:nr
%                    plot(1:100,squeeze(CumRanks(2,i,:)),str{i},'LineWidth',1);
%                end
%                figure;hold on;title('Verification RANDOM');
%                xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
%                plot(0:0.01:1,0:0.01:1,'k--');
%                for i=1:1:nr
%                    plot(rX(i,:),rY(i,:),str{i},'LineWidth',1);
%                end
%                figure;hold on;title('Verification TEST');
%                xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
%                plot(0:0.01:1,0:0.01:1,'k--');
%                for i=1:1:nr
%                    plot(tX(i,:),tY(i,:),str{i},'LineWidth',1);
%                end
%            end      
%end


% function out = getBiometricAcc(obj,tDepVar,tCOV,tGB)
%            nT = size(tDepVar,1);
%            tRIP = updateRIP(tDepVar,obj.MOD.M)';if obj.RIPNormalize, tRIP = BaseContainerv5.normalizeRIP(tRIP,obj.MOD.M,obj.MOD.Var);end
%            COVMatches = nan*zeros(obj.nrCOV,nT);
%            GBMatches = nan*zeros(obj.nrGB,nT);
%            parfor t=1:nT
%                 [~,COVD] = COVX2RIP(obj,tCOV(t,:),6,300);
%                 [~,gbD] = GBX2RIP(obj,tGB(t,:),6,300);
%                 COVMatches(:,t) = matchCOVRIP(obj,tRIP(t,1:obj.nrCOV),COVD);
%                 GBMatches(:,t) = matchGBRIP(obj,tRIP(t,obj.nrCOV+1:end),gbD);
%            end
%            Matches = [COVMatches;GBMatches];
%            out = ((sum(match2class(obj,Matches),2)/nT)*100)';      
% end

% function [scan,Cout,Atyp] = getBaseFace(obj,COVin,gbin,type,adjust,TrimBF)
%                 if nargin < 6, TrimBF = 2; end
%                 if nargin < 5, adjust = 'rescale'; end
%                 if nargin < 4, type = 'X'; end
%                 scan = [];C = []; Atyp = [];
%                 X = [COVin, gbin];
%                 if isempty(obj.MOD), return; end
%                 if isempty(obj.REF), return; end
%                 switch type
%                     case 'X'
%                         Xref = [obj.REF.AvgCOV, obj.REF.AvgGB];
%                         M = obj.MOD.MX;MI = obj.MOD.MIX;
%                     case 'RIP'
%                         Xref = [obj.REF.COVRIP, obj.REF.GBRIP];
%                         M = obj.MOD.MRIP;MI = obj.MOD.MIRIP;
%                 end
%                 X(isnan(X)) = Xref(isnan(X));
%                 C = [1 X]*[MI;M];
%                 Atyp = sqrt(sum(C.^2));
%                 Cout = C;
%                 switch lower(adjust);
%                     case 'rescale'
%                         C = C.*obj.ShapeSpace.EigStd';
%                     case 'trim'
%                         B = TrimBF*obj.ShapeSpace.EigStd';
%                         index = find(abs(C)>B);
%                         C(index) = B(index).*sign(C(index));
%                     otherwise
%                 end
%                 scan = getScan(obj.ShapeSpace,C);        
%        end
%        function [RIP,Distr] = COVX2RIP(obj,COVin,kappa,K)
%                 RIP = nan*zeros(1,obj.nrCOV);Distr = cell(1,obj.nrCOV);
%                 if isempty(obj.MOD), return; end
%                 if nargin < 4, K = 300; end
%                 if nargin < 3, kappa = 3; end 
%                 COVin(isnan(COVin)) = obj.AvgCOV(isnan(COVin));
%                 % Sex first as the only categorical case
%                 if ~(abs(COVin(1))-1==0), COVin(1) = 1; end% default choice is female
%                 Tind = find(obj.COV(:,1)==COVin(1));%Distr{1}.Tind = Tind;% Inlier Distribution
%                 Find = find(obj.COV(:,1)==-1*COVin(1));%Distr{1}.Find = Find;% Outlier Distribution
%                 [Distr{1}.Tmu,Distr{1}.Tsigma,Distr{1}.TO] = BaseContainerv5.extractRobustRIPStatistics(obj.COVRIP(Tind,1)',kappa);
%                 RIP(1) = Distr{1}.Tmu;Distr{1}.TpMax = normpdf(Distr{1}.Tmu,Distr{1}.Tmu,Distr{1}.Tsigma);
%                 [Distr{1}.Fmu,Distr{1}.Fsigma,Distr{1}.FO] = BaseContainerv5.extractRobustRIPStatistics(obj.COVRIP(Find,1)',kappa);
%                 Distr{1}.FpMax = normpdf(Distr{1}.Fmu,Distr{1}.Fmu,Distr{1}.Fsigma);
%                 % Age Weight and Height are Continuous
%                 for i=2:1:obj.nrCOV
%                     ind = find(obj.COV(:,i)>0);% eliminate nans and wrong data entries (typically coded as 0);
%                     X = obj.COV(ind,i);
%                     W = abs(X-COVin(i));
%                     % select closest samples
%                     [~,ind2] = sort(W,'ascend');
%                     % TRUE DISTRIBUTION
%                     ind2T = ind2(1:K);indT = ind(ind2T);
%                     rip = obj.COVRIP(indT,i)';
%                     Distr{i}.IndT = indT;
%                     WT = W(ind2T);WT = WT-min(WT);WT = WT/max(WT);WT = 1-WT';
%                     [Distr{i}.Tmu,Distr{i}.Tsigma,Distr{i}.TO] = BaseContainerv5.extractRobustRIPStatistics(rip,kappa,WT);
%                     Distr{i}.TpMax = normpdf(Distr{i}.Tmu,Distr{i}.Tmu,Distr{i}.Tsigma);
%                     RIP(i) = Distr{i}.Tmu;
%                     % FALSE DISTRIBUTION
%                     ind2F = ind2(end-K+1:end);indF = ind(ind2F);
%                     rip = obj.COVRIP(indF,i)';
%                     Distr{i}.IndF = indF;
%                     WF = W(ind2F);WF = WF-min(WF);WF = WF/max(WF);WF = WF';
%                     [Distr{i}.Fmu,Distr{i}.Fsigma,Distr{i}.FO] = BaseContainerv5.extractRobustRIPStatistics(rip,kappa,WF);
%                     Distr{i}.FpMax = normpdf(Distr{i}.Fmu,Distr{i}.Fmu,Distr{i}.Fsigma);
%                 end
%        end
%        function [RIP,Distr] = GBX2RIP(obj,gbin,kappa,K)
%                 RIP = nan*zeros(1,obj.nrGB);Distr = cell(1,obj.nrGB);
%                 if isempty(obj.MOD), return; end
%                 if nargin < 4, K = 300; end
%                 if nargin < 3, kappa = 3; end 
%                 gbin(isnan(gbin)) = obj.AvgGB(isnan(gbin));
%                 % Genetic Background axis are Continuous
%                 for i=1:1:obj.nrGB
%                     ind = find(~isnan(obj.GB(:,i)));% eliminate nans and wrong data entries (typically coded as 0);
%                     X = obj.GB(ind,i);
%                     W = abs(X-gbin(i));
%                     % select closest samples
%                     [~,ind2] = sort(W,'ascend');
%                     % TRUE DISTRIBUTION
%                     ind2T = ind2(1:K);indT = ind(ind2T);
%                     rip = obj.GBRIP(indT,i)';
%                     Distr{i}.IndT = indT;
%                     WT = W(ind2T);WT = WT-min(WT);WT = WT/max(WT);WT = 1-WT';
%                     [Distr{i}.Tmu,Distr{i}.Tsigma,Distr{i}.TO] = BaseContainerv5.extractRobustRIPStatistics(rip,kappa,WT);
%                     Distr{i}.TpMax = normpdf(Distr{i}.Tmu,Distr{i}.Tmu,Distr{i}.Tsigma);
%                     RIP(i) = Distr{i}.Tmu;
%                     % FALSE DISTRIBUTION
%                     ind2F = ind2(end-K+1:end);indF = ind(ind2F);
%                     rip = obj.GBRIP(indF,i)';
%                     Distr{i}.IndF = indF;
%                     WF = W(ind2F);WF = WF-min(WF);WF = WF/max(WF);WF = WF';
%                     [Distr{i}.Fmu,Distr{i}.Fsigma,Distr{i}.FO] = BaseContainerv5.extractRobustRIPStatistics(rip,kappa,WF);
%                     Distr{i}.FpMax = normpdf(Distr{i}.Fmu,Distr{i}.Fmu,Distr{i}.Fsigma);
%                 end
%        end
%        function [IB,pTIB,pT] = matchFaces(obj,faces,COVdistr,gbdistr)
%                 rip = updateRIP(faces,obj.MOD.M)';
%                 rip = BaseContainerv5.normalizeRIP(rip,obj.MOD.M,obj.MOD.Var);
%                 [COVmatches1,COVmatches2] = matchCOVRIP(obj,rip(:,1:obj.nrCOV,:),COVdistr);
%                 [gbmatches1,gbmatches2] = matchGBRIP(obj,rip(:,obj.nrCOV+1:end),gbdistr);
%                 pT = [COVmatches1;gbmatches1];
%                 IB = [COVmatches2;gbmatches2];
%                 pTIB = pT.*IB;
%        end
%        function [matches,matches1,matches2,matches3] = matchCOVRIP(obj,rip,distr)
%                 % initialize
%                   n = size(rip,1);
%                   matches1 = nan*zeros(obj.nrCOV,n);
%                   matches2 = nan*zeros(obj.nrCOV,n);
%                   matches3 = nan*zeros(obj.nrCOV,n);
%                   for i=1:1:obj.nrCOV
%                       pT = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)';
%                       pF = normpdf(rip(:,i),distr{i}.Fmu,distr{i}.Fsigma)';
%                       matches1(i,:) = (pT./distr{i}.TpMax);
%                       matches2(i,:) = (pT./(pT+pF));
%                       matches3(i,:) = pT./pF;
%                   end
%                   switch obj.Match
%                       case 'QD'
%                           matches = matches3;
%                       case 'INLIER'
%                           matches = matches2;
%                       otherwise
%                           matches = matches1;
%                   end
%        end
%        function [matches,matches1,matches2,matches3] = matchGBRIP(obj,rip,distr)
%                 % initialize
%                   n = size(rip,1);
%                   matches1 = nan*zeros(obj.nrGB,n);
%                   matches2 = nan*zeros(obj.nrGB,n);
%                   matches3 = nan*zeros(obj.nrGB,n);
%                   for i=1:1:obj.nrGB
%                       pT = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)';
%                       pF = normpdf(rip(:,i),distr{i}.Fmu,distr{i}.Fsigma)';
%                       matches1(i,:) = (pT./distr{i}.TpMax);
%                       matches2(i,:) = (pT./(pT+pF));
%                       matches3(i,:) = pT./pF;
%                   end
%                   switch obj.Match
%                       case 'QD'
%                           matches = matches3;
%                       case 'INLIER'
%                           matches = matches2;
%                       otherwise
%                           matches = matches1;
%                   end
%        end
%        function [matches1,matches2,matches3] = inversematchCOVRIP(obj,rip,distr)
%                 % initialize
%                   n = size(rip,1);
%                   matches1 = nan*zeros(obj.nrCOV,n);
%                   matches2 = nan*zeros(obj.nrCOV,n);
%                   matches3 = nan*zeros(obj.nrCOV,n);
%                   for i=1:1:obj.nrCOV
%                       pT = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)';
%                       pF = normpdf(rip(:,i),distr{i}.Fmu,distr{i}.Fsigma)';
%                       matches1(i,:) = (pF./distr{i}.FpMax);
%                       matches2(i,:) = (pF./(pT+pF));
%                       matches3(i,:) = pF./pT;
%                   end
%        end
%        function [matches1,matches2,matches3] = inversematchGBRIP(obj,rip,distr)
%                 % initialize
%                   n = size(rip,1);
%                   matches1 = nan*zeros(obj.nrGB,n);
%                   matches2 = nan*zeros(obj.nrGB,n);
%                   matches3 = nan*zeros(obj.nrGB,n);
%                   for i=1:1:obj.nrGB
%                       pT = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)';
%                       pF = normpdf(rip(:,i),distr{i}.Fmu,distr{i}.Fsigma)';
%                       matches1(i,:) = (pF./distr{i}.FpMax);
%                       matches2(i,:) = (pF./(pT+pF));
%                       matches3(i,:) = pF./pT;
%                   end
%        end
%        function out = randomFaces(obj,nr,type,norm,R)
%             if nargin < 5, R = 2.5; end
%             if nargin < 4, norm = true; end
%             if nargin < 3, type = 'Normal'; end
%             if nargin < 2, nr = 1000; end
%             switch norm
%                 case 1
%                     Std = ones(size(obj.ShapeSpace.EigStd));
%                 case 0
%                     Std = obj.ShapeSpace.EigStd;
%             end
%             switch type
%                 case 'Uniform'
%                     R = [obj.ShapeSpace.AvgCoeff-R*Std obj.ShapeSpace.AvgCoeff+R*Std];
%                     out = (repmat(R(:,1),1,nr) + repmat(R(:,2)-R(:,1),1,nr).*rand(obj.ShapeSpace.nrEV,nr))';
%                 case 'Normal'
%                     out = (repmat(obj.ShapeSpace.AvgCoeff,1,nr) + repmat(Std,1,nr).*randn(obj.ShapeSpace.nrEV,nr))';
%                 case 'Training'
%                     ind = randsample(obj.ShapeSpace.n,nr,true);
%                     out = obj.ShapeSpace.Tcoeff(ind,:);
%                     if norm, out = out./repmat(obj.ShapeSpace.EigStd',length(ind),1); end
%                 otherwise
%                     error('WRONG INITIALIZATION PROCEDURE');
%             end
%        end
%        function [Mscan,Pscan] = createRIPMorphsCOV(obj,type,nr)
%                 COVref = obj.REF.COVRIP;
%                 GBref = obj.REF.GBRIP;
%                 Yref = BaseCont.REF.RedCoeff;
%                 COVM = COVref;COVP = COVref;
%                 GBM = GBref;GBP = GBref;
%                 switch type
%                     case 'COV'
%                         if nr==1
%                         else
%                         end
%                     case 'gb'
%                 end
%                 
%                 maxcounter = 100;counter = 0;
%                 while Atypicality>=maxAtypicality && counter<=maxcounter
%                       counter = counter+1;
%                       X = Xref;X(end) = GTInfo.AVG+GTInfo.Sign*BF*GTInfo.STD;
%                       deltaX = X-Xref;dY = deltaX*obj.MOD.MMorphs;
%                       C = Yref+dY;
%                       Atypicality = sqrt(sum(C.^2));
%                       BF = BF-0.01;
%                 end
%                 if obj.ST4RedShape
%                    ShapeSpace = BaseCont.RedShapeSpace;
%                 else
%                    ShapeSpace = BaseCont.ShapeSpace;
%                 end
%                 scan = getScan(ShapeSpace,C);
%        end

