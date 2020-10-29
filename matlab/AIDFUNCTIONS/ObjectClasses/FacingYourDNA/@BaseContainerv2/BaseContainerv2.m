classdef BaseContainer < superClass
   % A container class containing all the parameter values
   properties
      TrackID = 'X';
      ShapeSpace = [];
      RedShapeSpace = [];
      RedShapeSpaceType = 'RIP'; % RIP (using ripped values) or X (using original values)
      Cov = [];
      CovNames = {};
      GB = [];% genetic background axes
   end
   properties (Dependent = true)
       GBNames;
       CovRIP;
       GBRIP;
   end
   properties (Dependent = true, Hidden = true)
      nrCov;
      nrGB;
      DepVar;
      OrigDepVar;
      RedDepVar;
      OrigRedDepVar;
      AvgCov;
      AvgGB;
      AvgCovRIP;
      AvgGBRIP;
      n;
      ShapeDim;
      RedShapeDim;
      indCov;
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
       ST3PT = 0.05;
   end
   properties (Hidden = true, Dependent = true)
       ST3Code;
   end
   properties
       ST4MOD = [];
       ST4Ref = [];
   end
   properties (Hidden = true)% (NGBRIM)
      RegRuns = 50; % Number of Regression runs (NGBRIM)
      OuterFold = 16; % number of Outer Folds (NGBRIM)
      InnerFold = 10; % number of Inner Folds (NGBRIM)
      MaxIterations = 5; % Maximum number of BRIM iterations (NGBRIM)
      Htest = true; % Perform Wilcoxon test om partial regression coefficients (NGBRIM)
      RIPNormalize = true; % normalize RIP values based on group averages (NGBRIM)
      StopCorr = 0.98; % Stopping correlation between subsequent iterations (NGBRIM)
      RIPStatPerm = 1000; % Number of permutations in testing significance of RIP values (NGBRIM)
      ParOuter = true;% execute parallel computation on outer fold or not (inner fold then) (NGBRIM)
   end
   methods % Constructor
        function obj = BaseContainer(varargin)
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
       function out = get.nrCov(obj)
                out = size(obj.Cov,2); 
       end
       function out = get.indCov(obj)
                if isempty(obj.Cov), out = []; return; end
                index = (1:size(obj.Cov,1));
                [i,~] = find(isnan(obj.Cov));
                i = unique(i);out = setdiff(index,i);
       end
       function out = get.nrGB(obj)
                out = size(obj.GB,2); 
       end
       function out = get.AvgCov(obj)
           if isempty(obj.Cov), out = []; return; end
           out = nanmean(obj.Cov);
       end
       function out = get.AvgGB(obj)
           if isempty(obj.GB), out = []; return; end
           out = nanmean(obj.GB);
       end
       function out = get.AvgCovRIP(obj)
           if isempty(obj.CovRIP), out = []; return; end
           out = nanmean(obj.CovRIP);
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
       function out = get.CovRIP(obj)
           if isempty(obj.ST4MOD), out = []; return; end
           out = obj.ST4MOD.RIP(:,1:obj.nrCov);
       end
       function out = get.GBRIP(obj)
           if isempty(obj.ST4MOD), out = []; return; end
           out = obj.ST4MOD.RIP(:,obj.nrCov+1:end);
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
   methods % GENERAL INTERFACE FUNCTIONS
       function createReducedSpace(obj)
           if isempty(obj.ShapeSpace), return; end
           switch obj.RedShapeSpaceType
               case 'X'
                   if isempty(obj.Cov), return; end
                   if isempty(obj.GB), return; end
                   A = [obj.Cov, obj.GB];
               case 'RIP'
                   if isempty(obj.CovRIP), return; end
                   if isempty(obj.GBRIP), return; end
                   A = [obj.CovRIP, obj.GBRIP];
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
       function out = NGBRIM(obj,DepVar,IndVar,BootInd,CondInd,PartInd)
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
                     el = unique(v(~isnan(v)));
                     if length(el)<6%Categorical parameter, currently allowing up to 6 categories  
                        var.Info{i}.Type = 'Categorical';
                        var.Info{i}.El = el;var.Info{i}.Mel = min(el);var.Info{i}.Pel = max(el);
                        var.Info{i}.Range = var.Info{i}.Pel-var.Info{i}.Mel;
                        % getting categorical reference faces
                        var.Info{i}.Mindex = find(v==var.Info{i}.Mel);
                        var.Info{i}.MDepVar = mean(DepVar(var.Info{i}.Mindex,:)); %#ok<*FNDSB>
                        var.Info{i}.Pindex = find(v==var.Info{i}.Pel);
                        var.Info{i}.PDepVar = mean(DepVar(var.Info{i}.Pindex,:));
                     else% variable is continous
                        var.Info{i}.Type = 'Continous';
                        index = find(~isnan(v));
                        [sortv,ind] = sort(v(index),'ascend');
                        % lower 25%
                        indl = round(0.25*length(ind));
                        var.Info{i}.Mb = sortv(indl);
                        var.Info{i}.Mindex = find(v<=var.Info{i}.Mb);
                        var.Info{i}.Mel = mean(v(var.Info{i}.Mindex));
                        var.Info{i}.MDepVar = mean(DepVar(var.Info{i}.Mindex,:));
                        % upper 25%
                        indu = round(0.75*length(ind));
                        var.Info{i}.Pb = sortv(indu);
                        var.Info{i}.Pindex = find(v>=var.Info{i}.Pb);
                        var.Info{i}.Pel = mean(v(var.Info{i}.Pindex));
                        var.Info{i}.PDepVar = mean(DepVar(var.Info{i}.Pindex,:));
                        var.Info{i}.Range = var.Info{i}.Pel-var.Info{i}.Mel;
                     end     
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
                                [M,H] = BaseContainer.getBootSampleRegression(FoTrIndVar(FiTrInd,:),FoTrDepVar(FiTrInd,:),var,RegRuns);
                                if Htest, M = H.*M;end% setting close to zero partial regression coefficients to zero
                                rip = BaseContainer.updateRIP(FoTrDepVar(FiTestInd,:),M)';
                                if RIPNormalize, rip = BaseContainer.normalizeRIP(rip,M,var);end
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
                        [M,H] = BaseContainer.getBootSampleRegression(FoTrIndVar,FoTrDepVar,var,RegRuns);
                        if Htest, M = H.*M;end
                        FoldResults{fo}.M = M;FoldResults{fo}.H = H;
                        % estimating rip values outer test data
                        rip = BaseContainer.updateRIP(FoTestDepVar,M)';
                        if RIPNormalize, rip = BaseContainer.normalizeRIP(rip,M,var);end
                        FoldResults{fo}.RIP = rip;
                        FoldResults{fo}.IndVar = FoTestIndVar(:,var.Booting);
                        % perform within outer fold statistics
                        [FoldResults{fo}.stat,FoldResults{fo}.statType,FoldResults{fo}.statP] = BaseContainer.performStat(FoldResults{fo}.IndVar,rip,var,RIPStatPerm);
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
                      out.pfast(i) = BaseContainer.pfast(out.StatP(i,:));
                  end
                % Extracting final M 
                  out.M = median(M,3);
                  out.Var = var;
                  if ~Htest, return; end
                  if OuterFold>=8, 
                     out.H = BaseContainer.wilcoxonM(M,var);
                  else
                     out.H = ones(size(out.M)); 
                  end
                  out.M = out.H.*out.M;
       end
       function obj = reduceSamples(obj,index)
          if nargout==1, obj = clone(obj); end
          obj.ShapeSpace.Tcoeff = obj.ShapeSpace.Tcoeff(index,:);
          obj.Cov = obj.Cov(index,:);
          obj.GB = obj.GB(index,:);
          if ~isempty(obj.RedShapeSpace),obj.RedShapeSpace.Tcoeff = obj.RedShapeSpace.Tcoeff(index,:);end
          %if~isempty(obj.ST4MOD), obj.ST4MOD.RIP = obj.ST4MOD.RIP(index,:); end
          obj.ST4MOD = [];
       end
   end
   methods % STAGE INTERFACE FUNCTIONS    
       function runStage1(obj,t)
           E = BaseContainer.getResiduals(obj.Cov(obj.indCov,:),obj.ShapeSpace.Tcoeff(obj.indCov,:));
           ind = obj.indCov;Cov = obj.Cov(ind,:);gb = obj.GB(ind,:);nrGB = obj.nrGB; %#ok<*PROP>
           [path,ID] = setupParForProgress(nrGB);
           r = nan*zeros(1,nrGB);rp = nan*zeros(1,nrGB);n = length(ind);
           parfor i=1:1:nrGB
               T = gb(:,i);T = BaseContainer.getResiduals(Cov,T);
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
           index = find(obj.ST1Code);
           reduceGB(obj,index);
       end
       function runStage2(obj,t)
           ind = obj.indCov;Cov = obj.Cov(ind,:);E = obj.ShapeSpace.Tcoeff(ind,:);
           gb = obj.GB(ind,:);nrGB = obj.nrGB;all = 1:nrGB; %#ok<*PROP>
           [path,ID] = setupParForProgress(nrGB);
           r = nan*zeros(1,nrGB);rp = nan*zeros(1,nrGB);n = length(ind);
           parfor i=1:1:nrGB
               tr = setdiff(all,i);
               T = BaseContainer.getResiduals([Cov, gb(:,tr)],gb(:,i));
               ET = BaseContainer.getResiduals([Cov, gb(:,tr)],E);
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
       function runStage3(obj,t)
           obj.RIPStatPerm = t;obj.OuterFold = 2;
           out = NGBRIM(obj,obj.DepVar,[obj.Cov, obj.GB],1:obj.nrCov+obj.nrGB,[],1);
           obj.ST3P = out.pfast(obj.nrCov+1:end);
           obj.ST3Stat = mean(out.Stat(obj.nrCov+1:end,:),2)';
           index = find(obj.ST3Code);
           if length(index)==obj.nrGB; disp('All Survived');return;end
           % reduce and do it again, untill all survive
           disp('Reducing and Running Again');
           reduceGB(obj,index);runStage3(obj,t);
       end
       function runStage4(obj)
            obj.RIPStatPerm = 10;obj.OuterFold = 10;IndVar = [obj.Cov, obj.GB];
            out = NGBRIM(obj,obj.DepVar,IndVar,1:obj.nrCov+obj.nrGB,[],1);
            obj.ST4MOD.P = out.pfast;obj.ST4MOD.M = out.M;obj.ST4MOD.Var = out.Var;
          % Setting up reference scan 
            obj.ST4Ref.Scan = clone(obj.ShapeSpace.Average);obj.ST4Ref.Coeff = (obj.ShapeSpace.AvgCoeff./obj.ShapeSpace.EigStd)';
            obj.ST4Ref.AvgCov = obj.AvgCov;obj.ST4Ref.AvgGB = obj.AvgGB;
            RIPref = BaseContainer.updateRIP(obj.ST4Ref.Coeff,obj.ST4MOD.M)';
            RIPref = BaseContainer.normalizeRIP(RIPref,obj.ST4MOD.M,obj.ST4MOD.Var);
            obj.ST4Ref.CovRIP = RIPref(1:obj.nrCov);obj.ST4Ref.GBRIP = RIPref(obj.nrCov+1:end);
          % Population used to model things
            RIP = BaseContainer.updateRIP(obj.DepVar,obj.ST4MOD.M)';
            RIP = BaseContainer.normalizeRIP(RIP,obj.ST4MOD.M,obj.ST4MOD.Var);
            [obj.ST4MOD.RIPStat] = BaseContainer.performStat(IndVar,RIP,obj.ST4MOD.Var,10);
            obj.ST4MOD.RIP = RIP;
           % Obtaining regressions to contruct morphs
            [obj.ST4MOD.MX,obj.ST4MOD.MIX] = getRegression(IndVar,obj.DepVar,(1:size(IndVar,2)));
            [obj.ST4MOD.MRIP,obj.ST4MOD.MIRIP] = getRegression(obj.ST4MOD.RIP,obj.DepVar,(1:size(obj.ST4MOD.RIP,2)));
            %[obj.ST4MOD.MXOrig,obj.ST4MOD.MIXOrig] = getRegression(IndVar,obj.OrigDepVar,(1:size(IndVar,2)));
            %[obj.ST4MOD.MRIPOrig,obj.ST4MOD.MIRIPOrig] = getRegression(obj.ST4MOD.RIP,obj.OrigDepVar,(1:size(obj.ST4MOD.RIP,2)));
           % Creating Reduced Space
            createReducedSpace(obj);
            obj.ST4Ref.RedScan = clone(obj.RedShapeSpace.Average);obj.ST4Ref.RedCoeff = (obj.RedShapeSpace.AvgCoeff./obj.RedShapeSpace.EigStd)';
       end
   end
   methods % MODEL FUNCTION
       function [scan,Cout,Atyp] = getBaseFace(obj,covin,gbin,type,adjust,TrimBF)
                if nargin < 6, TrimBF = 2; end
                if nargin < 5, adjust = 'rescale'; end
                if nargin < 4, type = 'X'; end
                scan = [];C = []; Atyp = [];
                X = [covin, gbin];
                if isempty(obj.ST4MOD), return; end
                if isempty(obj.ST4Ref), return; end
                switch type
                    case 'X'
                        Xref = [obj.ST4Ref.AvgCov, obj.ST4Ref.AvgGB];
                        M = obj.ST4MOD.MX;MI = obj.ST4MOD.MIX;
                    case 'RIP'
                        Xref = [obj.ST4Ref.CovRIP, obj.ST4Ref.GBRIP];
                        M = obj.ST4MOD.MRIP;MI = obj.ST4MOD.MIRIP;
                end
                X(isnan(X)) = Xref(isnan(X));
                C = [1 X]*[MI;M];
                Atyp = sqrt(sum(C.^2));
                Cout = C;
                switch lower(adjust);
                    case 'rescale'
                        C = C.*obj.ShapeSpace.EigStd';
                    case 'trim'
                        B = TrimBF*obj.ShapeSpace.EigStd';
                        index = find(abs(C)>B);
                        C(index) = B(index).*sign(C(index));
                    otherwise
                end
                scan = getScan(obj.ShapeSpace,C);        
       end
       function [RIP,Distr] = CovX2RIP(obj,covin,kappa,K)
                RIP = nan*zeros(1,obj.nrCov);Distr = cell(1,obj.nrCov);
                if isempty(obj.ST4MOD), return; end
                if nargin < 4, K = 300; end
                if nargin < 3, kappa = 3; end 
                covin(isnan(covin)) = obj.AvgCov(isnan(covin));
                % Sex first as the only categorical case
                if ~(abs(covin(1))-1==0), covin(1) = 1; end% default choice is female
                Tind = find(obj.Cov(:,1)==covin(1));%Distr{1}.Tind = Tind;% Inlier Distribution
                Find = find(obj.Cov(:,1)==-1*covin(1));%Distr{1}.Find = Find;% Outlier Distribution
                [Distr{1}.Tmu,Distr{1}.Tsigma,Distr{1}.TO] = BaseContainer.extractRobustRIPStatistics(obj.CovRIP(Tind,1)',kappa);
                RIP(1) = Distr{1}.Tmu;Distr{1}.TpMax = normpdf(Distr{1}.Tmu,Distr{1}.Tmu,Distr{1}.Tsigma);
                [Distr{1}.Fmu,Distr{1}.Fsigma,Distr{1}.FO] = BaseContainer.extractRobustRIPStatistics(obj.CovRIP(Find,1)',kappa);
                Distr{1}.FpMax = normpdf(Distr{1}.Fmu,Distr{1}.Fmu,Distr{1}.Fsigma);
                % Age Weight and Height are continous
                for i=2:1:obj.nrCov
                    ind = find(obj.Cov(:,i)>0);% eliminate nans and wrong data entries (typically coded as 0);
                    X = obj.Cov(ind,i);
                    W = abs(X-covin(i));
                    % select closest samples
                    [~,ind2] = sort(W,'ascend');
                    %ind2 = ind2(1:K);ind = ind(ind2);
                    %rip = obj.CovRIP(ind,i)';
                    %Distr{i}.Ind = ind;
                    %W = W(ind2);W = W-min(W);W = W/max(W);W = 1-W';
                    %[Distr{i}.Tmu,Distr{i}.Tsigma,Distr{i}.TO] = BaseContainer.extractRobustRIPStatistics(rip,kappa,W);
                    %Distr{i}.TpMax = normpdf(Distr{i}.Tmu,Distr{i}.Tmu,Distr{i}.Tsigma);
                    %RIP(i) = Distr{i}.Tmu;
                    % TRUE DISTRIBUTION
                    ind2T = ind2(1:K);indT = ind(ind2T);
                    rip = obj.CovRIP(indT,i)';
                    Distr{i}.IndT = indT;
                    WT = W(ind2T);WT = WT-min(WT);WT = WT/max(WT);WT = 1-WT';
                    [Distr{i}.Tmu,Distr{i}.Tsigma,Distr{i}.TO] = BaseContainer.extractRobustRIPStatistics(rip,kappa,WT);
                    Distr{i}.TpMax = normpdf(Distr{i}.Tmu,Distr{i}.Tmu,Distr{i}.Tsigma);
                    RIP(i) = Distr{i}.Tmu;
                    % FALSE DISTRIBUTION
                    ind2F = ind2(end-K+1:end);indF = ind(ind2F);
                    rip = obj.CovRIP(indF,i)';
                    Distr{i}.IndF = indF;
                    WF = W(ind2F);WF = WF-min(WF);WF = WF/max(WF);WF = WF';
                    [Distr{i}.Fmu,Distr{i}.Fsigma,Distr{i}.FO] = BaseContainer.extractRobustRIPStatistics(rip,kappa,WF);
                    Distr{i}.FpMax = normpdf(Distr{i}.Fmu,Distr{i}.Fmu,Distr{i}.Fsigma);
                end
       end
       function [RIP,Distr] = GBX2RIP(obj,gbin,kappa,K)
                RIP = nan*zeros(1,obj.nrGB);Distr = cell(1,obj.nrGB);
                if isempty(obj.ST4MOD), return; end
                if nargin < 4, K = 300; end
                if nargin < 3, kappa = 3; end 
                gbin(isnan(gbin)) = obj.AvgGB(isnan(gbin));
                % Genetic Background axis are continous
                for i=1:1:obj.nrGB
                    ind = find(~isnan(obj.GB(:,i)));% eliminate nans and wrong data entries (typically coded as 0);
                    X = obj.GB(ind,i);
                    W = abs(X-gbin(i));
                    % select closest samples
                    [~,ind2] = sort(W,'ascend');
                    % TRUE DISTRIBUTION
                    ind2T = ind2(1:K);indT = ind(ind2T);
                    rip = obj.GBRIP(indT,i)';
                    Distr{i}.IndT = indT;
                    WT = W(ind2T);WT = WT-min(WT);WT = WT/max(WT);WT = 1-WT';
                    [Distr{i}.Tmu,Distr{i}.Tsigma,Distr{i}.TO] = BaseContainer.extractRobustRIPStatistics(rip,kappa,WT);
                    Distr{i}.TpMax = normpdf(Distr{i}.Tmu,Distr{i}.Tmu,Distr{i}.Tsigma);
                    RIP(i) = Distr{i}.Tmu;
                    % FALSE DISTRIBUTION
                    ind2F = ind2(end-K+1:end);indF = ind(ind2F);
                    rip = obj.GBRIP(indF,i)';
                    Distr{i}.IndF = indF;
                    WF = W(ind2F);WF = WF-min(WF);WF = WF/max(WF);WF = WF';
                    [Distr{i}.Fmu,Distr{i}.Fsigma,Distr{i}.FO] = BaseContainer.extractRobustRIPStatistics(rip,kappa,WF);
                    Distr{i}.FpMax = normpdf(Distr{i}.Fmu,Distr{i}.Fmu,Distr{i}.Fsigma);
                end
       end
       function [IB,pTIB,pT] = matchFaces(obj,faces,covdistr,gbdistr)
                rip = updateRIP(faces,obj.ST4MOD.M)';
                rip = BaseContainer.normalizeRIP(rip,obj.ST4MOD.M,obj.ST4MOD.Var);
                [covmatches1,covmatches2] = matchCovRIP(obj,rip(:,1:obj.nrCov,:),covdistr);
                [gbmatches1,gbmatches2] = matchGBRIP(obj,rip(:,obj.nrCov+1:end),gbdistr);
                pT = [covmatches1;gbmatches1];
                IB = [covmatches2;gbmatches2];
                pTIB = pT.*IB;
       end
       function [matches1,matches2,matches3] = matchCovRIP(obj,rip,distr)
                % initialize
                  n = size(rip,1);
                  matches1 = nan*zeros(obj.nrCov,n);
                  matches2 = nan*zeros(obj.nrCov,n);
                  matches3 = nan*zeros(obj.nrCov,n);
                  for i=1:1:obj.nrCov
                      pT = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)';
                      pF = normpdf(rip(:,i),distr{i}.Fmu,distr{i}.Fsigma)';
                      matches1(i,:) = (pT./distr{i}.TpMax);
                      matches2(i,:) = (pT./(pT+pF));
                      matches3(i,:) = pT./pF;
                  end
                % Gender first
                  %pT = normpdf(rip(:,1),distr{1}.Tmu,distr{1}.Tsigma)';
                  %pF = normpdf(rip(:,1),distr{1}.Fmu,distr{1}.Fsigma)';
                  %matches1(1,:) = (pT./(pT+pF));
                  %matches2(1,:) = (pT./(pT+pF));
                  %
                  %matches(1,:) = (pT./distr{1}.TpMax).*(pT./(pT+pF));
                % Now the rest
                  %for i=2:1:obj.nrCov   
                  %    matches1(i,:) = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)'./distr{i}.TpMax;
                  %end
       end
       function [matches1,matches2,matches3] = matchGBRIP(obj,rip,distr)
                % initialize
                  n = size(rip,1);
                  matches1 = nan*zeros(obj.nrGB,n);
                  matches2 = nan*zeros(obj.nrGB,n);
                  matches3 = nan*zeros(obj.nrGB,n);
                  for i=1:1:obj.nrGB
                      pT = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)';
                      pF = normpdf(rip(:,i),distr{i}.Fmu,distr{i}.Fsigma)';
                      matches1(i,:) = (pT./distr{i}.TpMax);
                      matches2(i,:) = (pT./(pT+pF));
                      matches3(i,:) = pT./pF;
                  end
                  %matches = nan*zeros(obj.nrGB,n);
                % Now the rest
                  %for i=1:1:obj.nrGB
                  %    matches(i,:) = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)'./distr{i}.TpMax;
                  %end
       end
       function [matches1,matches2,matches3] = inversematchCovRIP(obj,rip,distr)
                % initialize
                  n = size(rip,1);
                  matches1 = nan*zeros(obj.nrCov,n);
                  matches2 = nan*zeros(obj.nrCov,n);
                  matches3 = nan*zeros(obj.nrCov,n);
                  for i=1:1:obj.nrCov
                      pT = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)';
                      pF = normpdf(rip(:,i),distr{i}.Fmu,distr{i}.Fsigma)';
                      matches1(i,:) = (pF./distr{i}.FpMax);
                      matches2(i,:) = (pF./(pT+pF));
                      matches3(i,:) = pF./pT;
                  end
       end
       function [matches1,matches2,matches3] = inversematchGBRIP(obj,rip,distr)
                % initialize
                  n = size(rip,1);
                  matches1 = nan*zeros(obj.nrGB,n);
                  matches2 = nan*zeros(obj.nrGB,n);
                  matches3 = nan*zeros(obj.nrGB,n);
                  for i=1:1:obj.nrGB
                      pT = normpdf(rip(:,i),distr{i}.Tmu,distr{i}.Tsigma)';
                      pF = normpdf(rip(:,i),distr{i}.Fmu,distr{i}.Fsigma)';
                      matches1(i,:) = (pF./distr{i}.FpMax);
                      matches2(i,:) = (pF./(pT+pF));
                      matches3(i,:) = pF./pT;
                  end
       end
       function out = randomFaces(obj,nr,type,norm,R)
            if nargin < 5, R = 2.5; end
            if nargin < 4, norm = true; end
            if nargin < 3, type = 'Normal'; end
            if nargin < 2, nr = 1000; end
            switch norm
                case 1
                    Std = ones(size(obj.ShapeSpace.EigStd));
                case 0
                    Std = obj.ShapeSpace.EigStd;
            end
            switch type
                case 'Uniform'
                    R = [obj.ShapeSpace.AvgCoeff-R*Std obj.ShapeSpace.AvgCoeff+R*Std];
                    out = (repmat(R(:,1),1,nr) + repmat(R(:,2)-R(:,1),1,nr).*rand(obj.ShapeSpace.nrEV,nr))';
                case 'Normal'
                    out = (repmat(obj.ShapeSpace.AvgCoeff,1,nr) + repmat(Std,1,nr).*randn(obj.ShapeSpace.nrEV,nr))';
                case 'Training'
                    ind = randsample(obj.ShapeSpace.n,nr,true);
                    out = obj.ShapeSpace.Tcoeff(ind,:);
                    if norm, out = out./repmat(obj.ShapeSpace.EigStd',length(ind),1); end
                otherwise
                    error('WRONG INITIALIZATION PROCEDURE');
            end
       end
       function [Mscan,Pscan] = createRIPMorphsCov(obj,type,nr)
                Covref = obj.ST4Ref.CovRIP;
                GBref = obj.ST4Ref.GBRIP;
                Yref = BaseCont.ST4Ref.RedCoeff;
                CovM = Covref;CovP = Covref;
                GBM = GBref;GBP = GBref;
                switch type
                    case 'cov'
                        if nr==1
                        else
                        end
                    case 'gb'
                end
                
                maxcounter = 100;counter = 0;
                while Atypicality>=maxAtypicality && counter<=maxcounter
                      counter = counter+1;
                      X = Xref;X(end) = GTInfo.AVG+GTInfo.Sign*BF*GTInfo.STD;
                      deltaX = X-Xref;dY = deltaX*obj.ST4MOD.MMorphs;
                      C = Yref+dY;
                      Atypicality = sqrt(sum(C.^2));
                      BF = BF-0.01;
                end
                if obj.ST4RedShape
                   ShapeSpace = BaseCont.RedShapeSpace;
                else
                   ShapeSpace = BaseCont.ShapeSpace;
                end
                scan = getScan(ShapeSpace,C);
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
       function [stat,statType,statP] = performStat(vIn,vripIn,var,t)
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
                           [stat(s),~,statP(s)] = myPermAnova(XG,vrip(index),t);
                       case 'Continous'
                           statType{s} = 'Cor';
                           [stat(s),statP(s)] = permCorr(v(index),vrip(index),t);
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
       function [M,H] = getBootSampleRegression(IndVar,DepVar,var,runs)
                nrS = size(IndVar,1);
                M = zeros(var.nr2Boot,size(DepVar,2),runs);
                for s=1:1:runs
                    if runs==1% no random sampling
                       SampleFi = 1:nrS; 
                    else % random sampling with replacement
                       SampleFi = randsample(nrS,nrS,true);
                    end
                    M(:,:,s) = BaseContainer.getRegression(IndVar(SampleFi,:),DepVar(SampleFi,:),var.Booting);
                end
                if runs==1,M = squeeze(M);H = ones(size(M));return;end
                H = BaseContainer.wilcoxonM(M,var);
                M = median(M,3);       
       end
       function H = wilcoxonM(M,var)
            H = zeros(var.nr2Boot,size(M,2));
            for i=1:1:var.nr2Boot
               Mfor = squeeze(M(i,:,:)); 
               for j=1:1:size(Mfor,1)
                  [~,H(i,j)] = signrank(Mfor(j,:));    
               end
            end
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
   end
end

% for i=1:1:obj.nrGB
%                     ind = find(~isnan(obj.GB(:,i)));% eliminate nans and wrong data entries (typically coded as 0);
%                     X = obj.GB(ind,i);
%                     W = abs(X-gbin(i));
%                     % select closest samples
%                     [~,ind2] = sort(W,'ascend');
%                     ind2 = ind2(1:K);ind = ind(ind2);
%                     rip = obj.GBRIP(ind,i)';
%                     Distr{i}.Ind = ind;
%                     W = W(ind2);W = W-min(W);W = W/max(W);W = 1-W';
%                     [Distr{i}.Tmu,Distr{i}.Tsigma,Distr{i}.TO] = BaseContainer.extractRobustRIPStatistics(rip,kappa,W);
%                     Distr{i}.TpMax = normpdf(Distr{i}.Tmu,Distr{i}.Tmu,Distr{i}.Tsigma);
%                     RIP(i) = Distr{i}.Tmu;
%                 end

% for i=2:1:obj.nrCov
%                     ind = find(obj.Cov(:,i)>0);% eliminate nans and wrong data entries (typically coded as 0);
%                     X = obj.Cov(ind,i);
%                     W = abs(X-covin(i));
%                     % select closest samples
%                     [~,ind2] = sort(W,'ascend');
%                     ind2 = ind2(1:K);ind = ind(ind2);
%                     rip = obj.CovRIP(ind,i)';
%                     Distr{i}.Ind = ind;
%                     W = W(ind2);W = W-min(W);W = W/max(W);W = 1-W';
%                     [Distr{i}.Tmu,Distr{i}.Tsigma,Distr{i}.TO] = BaseContainer.extractRobustRIPStatistics(rip,kappa,W);
%                     Distr{i}.TpMax = normpdf(Distr{i}.Tmu,Distr{i}.Tmu,Distr{i}.Tsigma);
%                     RIP(i) = Distr{i}.Tmu;
%                 end



%Yref = obj.ST4Ref.Coeff;
                %deltaX = X-Xref;
                %dY = deltaX*M;
                %C = Yref+dY;
                %Atyp = sqrt(sum(C.^2));
                %scan = getScan(obj.ShapeSpace,C);
