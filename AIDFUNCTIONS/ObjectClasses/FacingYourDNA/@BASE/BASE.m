classdef BASE < superClass
    properties
       TrackID = 'X';
       ShapeSpace = [];
       ShapeIndex = [];
       COV = [];
       COVNames = {'Sex' 'Age' 'Weight' 'Height'};
       COVMOD = [];
       GB = [];
       GBindex = [];
       GBMOD = {};
       REF = [];
       GBMultiple = false;
       COVMultiple = true;
    end
    properties (Dependent = true)
       GBNames;
       nCOV;
       nGB;
       nCOVGB;
       nS;
       DepVar;
       AvgCOV;
       AvgGB;
       ShapeDim;
    end
    properties
        BuildMethod = 'PLSR';
        RegRuns = 100;
        SamplePerc = 1;
        HTest = true;
        HTestp = 0.0001;
        RIPNorm = true;
    end
    properties % TESTING
        COVR2Test = [];
        COVATest = [];
        COVBioTest = [];
        GBR2Test = [];
        GBATest = [];
        GBBioTest = [];
    end
    properties % BIOMETRICS
        Match = 'QD';
        KN = 100;
        Kappa = 6;
        Classify = 'SOFT';
    end
    methods % Constructor
        function obj = BASE(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % GENERAL GETTING/SETTING
       function out = get.GBNames(obj)
           out = cell(1,obj.nGB);
           for i=1:1:obj.nGB
               out{i} = ['PC' num2str(obj.GBindex(i))];
           end
       end
       function out = get.nCOV(obj)
                out = size(obj.COV,2); 
       end
       function out = get.nGB(obj)
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
       function out = get.ShapeSpace(obj)
           out = obj.ShapeSpace;
           if ~superClass.isH(out), out = []; end
       end
       function out = get.nS(obj)
          if isempty(obj.ShapeSpace), out = 0; return; end
          out = length(obj.ShapeIndex);
       end 
       function out = get.DepVar(obj)
            if isempty(obj.ShapeSpace), out = [];return; end
            out = obj.ShapeSpace.Tcoeff(obj.ShapeIndex,:)./repmat(obj.ShapeSpace.EigStd',obj.nS,1);
       end
       function out = get.ShapeDim(obj)
           if isempty(obj.ShapeSpace), out = 0; return; end
           out = obj.ShapeSpace.nrEV;
       end
       function out = get.nCOVGB(obj)
           out = obj.nCOV+obj.nGB;
       end
   end
   methods % COVARIATE ANALYSES
       function out = runR2TestCOV(obj,GBindex,t)
           R2 = nan*zeros(1,obj.nCOV);
           pR2 = nan*zeros(1,obj.nCOV);
           [path,ID] = setupParForProgress(obj.nCOV);
           all = 1:obj.nCOV;
           if isempty(GBindex)
              GB = []; %#ok<*PROP,*PROPLC>
           else
              GB = obj.GB(:,index);
           end
           parfor i=1:1:obj.nCOV
               switch obj.COVMultiple
                   case false
                       [R2(i),pR2(i)] = BASE.testPartialPLSR(obj.COV(:,i),obj.DepVar,GB,t);
                   case true
                       rm = setdiff(all,i);
                       [R2(i),pR2(i)] = BASE.testPartialPLSR(obj.COV(:,i),obj.DepVar,[obj.COV(:,rm) GB],t);
               end
               parfor_progress;
           end
           closeParForProgress(path,ID);
           out.R2 = R2;out.pR2 = pR2;
           obj.COVR2Test = out;
       end
       function out = runAngleTestCOV(obj,t,maxM,Aval,GBindex)
           M = nan*zeros(obj.ShapeDim,2,maxM,obj.nCOV);
           A = nan*zeros(t,obj.nCOV);
           pA = nan*zeros(length(Aval),obj.nCOV);
           avgA = nan*zeros(1,obj.nCOV);
           stdA = nan*zeros(1,obj.nCOV);
           medA = nan*zeros(1,obj.nCOV);
           madA = nan*zeros(1,obj.nCOV);
           lowerciA = nan*zeros(1,obj.nCOV);
           upperciA = nan*zeros(1,obj.nCOV);
           [path,ID] = setupParForProgress(obj.nCOV);
           if isempty(GBindex)
              GB = []; %#ok<*PROPLC>
           else
              GB = obj.GB(:,index);
           end
           all = 1:obj.nCOV;
           parfor a=1:obj.nCOV
               diff = setdiff(all,a);
               [tmpA,tmpM,STAT] = AngleTest(obj,[obj.COV(:,diff) GB obj.COV(:,a)],obj.DepVar,Aval,maxM,t);
               A(:,a) = tmpA;M(:,:,:,a) = tmpM;pA(:,a) = STAT.pA;avgA(a) = STAT.avgA;
               medA(a) = STAT.medA;stdA(a) = STAT.stdA;madA(a) = STAT.madA;
               upperciA(a) = STAT.upperciA;lowerciA(a) = STAT.lowerciA; 
               parfor_progress;
           end
           closeParForProgress(path,ID);
           out.A = A;out.M = M;
           out.AvgA = avgA;out.MedA = medA;out.StdA = stdA;out.MadA = madA;out.LowA = lowerciA;out.UpA = upperciA;out.pA = pA;
           out.Aval = Aval;
           obj.COVATest = out;
       end
       function buildFromAngleTestCOV(obj,test)
           out = cell(1,obj.nCOV);
           for i=1:obj.nCOV
               out{i}.IndVar = obj.COV(:,i);
               forM = squeeze(test.M(:,:,:,i));
             % VARIABLE INFORMATION  
               out{i}.Var.nr2Boot = 1;
               out{i}.Var.Booting = obj.nCOV+1;
               out{i}.Var.Info{1} = BASE.getVarInfo(out{i}.IndVar,obj.DepVar);
               fullM = [squeeze(forM(:,1,:)) squeeze(forM(:,2,:))];
               out{i}.M = nanmedian(fullM,2)';
               out{i}.H = BASE.wilcoxonMSimple(fullM',obj.HTestp);
               out{i}.HM = out{i}.M.*out{i}.H;
              % RIP COMPUTATIONS
               out{i}.MRIP = BASE.updateRIP(obj.DepVar,out{i}.M)';
               out{i}.NormMRIP = BASE.normalizeRIP(out{i}.MRIP,out{i}.M,out{i}.Var);
               out{i}.HMRIP = BASE.updateRIP(obj.DepVar,out{i}.HM)';
               out{i}.NormHMRIP = BASE.normalizeRIP(out{i}.HMRIP,out{i}.HM,out{i}.Var);
           end
           obj.COVMOD = out;
       end
       function out = biometricTestCOV(obj,TestDepVar,TestCOV,nrF)
            EER = nan*zeros(1,obj.nCOV);
            CumRanks = nan*zeros(100,obj.nCOV);
            ClassACC = nan*zeros(1,obj.nCOV);
            nrT = size(TestDepVar,1);
            if nargin<5, nrF = 200; end
            [path,ID] = setupParForProgress(obj.nCOV);
            parfor i=1:obj.nCOV
                MOD = obj.COVMOD{i};
                testcov = TestCOV(:,i);
                RIP = BASE.updateRIP(TestDepVar,getMODM(obj,MOD))';
                if obj.RIPNorm, RIP = BASE.normalizeRIP(RIP,getMODM(obj,MOD),MOD.Var);end
                TMatches = nan*zeros(1,nrT);
                FMatches = nan*zeros(1,nrT,nrF);
                for t=1:nrT
                    diff = abs(testcov-testcov(t));
                    [~,Find] = sort(diff);
                    Find = Find(end-nrF+1:end);
                    Distr = BASE.X2RIP(MOD.IndVar,getMODRIP(obj,MOD),testcov(t),MOD.Var.Info{1}.Type,obj.Kappa,obj.KN);
                    m = getMatch(obj,Distr,RIP);
                    TMatches(1,t) = m(t);
                    FMatches(1,t,:) = m(Find);
                end
                TScores = TMatches;FScores = FMatches;
                switch obj.Classify
                   case 'SOFT'
                       TScores = match2score(obj,TMatches);
                       FScores = match2score(obj,FMatches);
                   case 'HARD'
                       TScores = BASE.class2score(match2class(obj,TMatches));
                       FScores = BASE.class2score(match2class(obj,FMatches));
                end
                ClassACC(i) = ((sum(match2class(obj,TMatches),2)/nrT)*100);
                [forRanks,EER(i)] = BASE.getBiometrics(TScores,FScores,ones(1,nrT),ones(size(ones(1,nrT),1),1));
                forRanks = (forRanks./(nrF+1)).*100;
                CumRank = zeros(1,100);
                for cr=1:1:100
                    CumRank(1,cr) = (sum(forRanks<=cr)/(nrT)).*100;
                end
                CumRanks(:,i) = CumRank';
                parfor_progress;
            end
            closeParForProgress(path,ID);
            out.EER = EER;
            out.CumRanks = CumRanks;
            out.ClassACC = ClassACC;
       end
   end
   methods % GENETIC BACKGROUND ANALYSES
       function out = runR2TestGB(obj,t)
           R2 = nan*zeros(1,obj.nGB);
           pR2 = nan*zeros(1,obj.nGB);
           [path,ID] = setupParForProgress(obj.nGB);
           all = 1:obj.nGB;
           parfor i=1:1:obj.nGB
               switch obj.GBMultiple
                   case false
                       [R2(i),pR2(i)] = BASE.testPartialPLSR(obj.GB(:,i),obj.DepVar,obj.COV,t);
                   case true
                       rm = setdiff(all,i);
                       [R2(i),pR2(i)] = BASE.testPartialPLSR(obj.GB(:,i),obj.DepVar,[obj.COV obj.GB(:,rm)],t);
               end
               parfor_progress;
           end
           closeParForProgress(path,ID);
           out.R2 = R2;out.pR2 = pR2;
           obj.GBR2Test = out;
       end  
       function out = runAngleTestGB(obj,t,maxM,Aval)
           if nargin<4, Aval = 0; end
           if nargin<3, maxM = t; end
           M = nan*zeros(obj.ShapeDim,2,maxM,obj.nGB);
           A = nan*zeros(t,obj.nGB);
           pA = nan*zeros(length(Aval),obj.nGB);
           avgA = nan*zeros(1,obj.nGB);
           stdA = nan*zeros(1,obj.nGB);
           medA = nan*zeros(1,obj.nGB);
           madA = nan*zeros(1,obj.nGB);
           lowerciA = nan*zeros(1,obj.nGB);
           upperciA = nan*zeros(1,obj.nGB);
           [path,ID] = setupParForProgress(obj.nGB);
           parfor a=1:obj.nGB
               [tmpA,tmpM,STAT] = AngleTest(obj,[obj.COV obj.GB(:,a)],obj.DepVar,Aval,maxM,t);
               A(:,a) = tmpA;M(:,:,:,a) = tmpM;pA(:,a) = STAT.pA;avgA(a) = STAT.avgA;
               medA(a) = STAT.medA;stdA(a) = STAT.stdA;madA(a) = STAT.madA;
               upperciA(a) = STAT.upperciA;lowerciA(a) = STAT.lowerciA;
               parfor_progress;
           end
           closeParForProgress(path,ID);
           out.A = A;out.M = M;
           out.AvgA = avgA;out.MedA = medA;out.StdA = stdA;out.MadA = madA;out.LowA = lowerciA;out.UpA = upperciA;out.pA = pA;
           out.Aval = Aval;
           obj.GBATest = out;
       end
       function buildFromAngleTestGB(obj,test)
           out = cell(1,obj.nGB);
           for i=1:obj.nGB
               out{i}.IndVar = obj.GB(:,i);
               forM = squeeze(test.M(:,:,:,i));
             % VARIABLE INFORMATION  
               out{i}.Var.nr2Boot = 1;
               out{i}.Var.Booting = obj.nCOV+1;
               out{i}.Var.Info{1} = BASE.getVarInfo(out{i}.IndVar,obj.DepVar);
               fullM = [squeeze(forM(:,1,:)) squeeze(forM(:,2,:))];
               out{i}.M = nanmedian(fullM,2)';
               out{i}.H = BASE.wilcoxonMSimple(fullM',obj.HTestp);
               out{i}.HM = out{i}.M.*out{i}.H;
              % RIP COMPUTATIONS
               out{i}.MRIP = BASE.updateRIP(obj.DepVar,out{i}.M)';
               out{i}.NormMRIP = BASE.normalizeRIP(out{i}.MRIP,out{i}.M,out{i}.Var);
               out{i}.HMRIP = BASE.updateRIP(obj.DepVar,out{i}.HM)';
               out{i}.NormHMRIP = BASE.normalizeRIP(out{i}.HMRIP,out{i}.HM,out{i}.Var);
           end
           obj.GBMOD = out;
       end
       function out = buildGB(obj)
           % TO BE UPDATED
           [path,ID] = setupParForProgress(obj.nGB);
           out = cell(1,obj.nGB);
           allGB = 1:obj.nGB;
           parfor i=1:obj.nGB
               switch obj.GBMultiple
                   case false
                       out{i} = build(obj,obj.DepVar,[obj.COV obj.GB(:,i)],obj.nCOV+1,1:obj.nCOV,obj.COV(:,1));
                   case true
                       rmGB = setdiff(allGB,i);
                       out{i} = build(obj,obj.DepVar,[obj.COV obj.GB(:,rmGB) obj.GB(:,i)],obj.nCOV+obj.nGB,1:(obj.nCOV+obj.nGB-1),obj.COV(:,1));
               end
               parfor_progress;
           end
           closeParForProgress(path,ID);
           obj.GBMOD = out;
       end
       function out = biometricTestGB(obj,TestDepVar,TestGB,nrF)
            EER = nan*zeros(1,obj.nGB);
            CumRanks = nan*zeros(100,obj.nGB);
            ClassACC = nan*zeros(1,obj.nGB);
            nrT = size(TestDepVar,1);
            if nargin<5, nrF = 200; end
            [path,ID] = setupParForProgress(obj.nGB);
            parfor i=1:obj.nGB
                MOD = obj.GBMOD{i};
                testgb = TestGB(:,i);
                RIP = BASE.updateRIP(TestDepVar,getMODM(obj,MOD))';
                if obj.RIPNorm, RIP = BASE.normalizeRIP(RIP,getMODM(obj,MOD),MOD.Var);end
                %RIP = BASE.updateRIP(TestDepVar,MOD.M)';
                %if obj.RIPNorm, RIP = BASE.normalizeRIP(RIP,MOD.M,MOD.Var);end
                TMatches = nan*zeros(1,nrT);
                FMatches = nan*zeros(1,nrT,nrF);
                for t=1:nrT
                    diff = abs(testgb-testgb(t));
                    [~,Find] = sort(diff);
                    Find = Find(end-nrF+1:end);
                    Distr = BASE.X2RIP(MOD.IndVar,getMODRIP(obj,MOD),testgb(t),'Continuous',obj.Kappa,obj.KN);
                    %Distr = BASE.X2RIP(MOD.IndVar,MOD.MRIP,testgb(t),'Continuous',obj.Kappa,obj.KN);
                    m = getMatch(obj,Distr,RIP);
                    TMatches(1,t) = m(t);
                    FMatches(1,t,:) = m(Find);
                end
                TScores = TMatches;FScores = FMatches;
                switch obj.Classify
                   case 'SOFT'
                       TScores = match2score(obj,TMatches);
                       FScores = match2score(obj,FMatches);
                   case 'HARD'
                       TScores = BASE.class2score(match2class(obj,TMatches));
                       FScores = BASE.class2score(match2class(obj,FMatches));
                end
                ClassACC(i) = ((sum(match2class(obj,TMatches),2)/nrT)*100);
                [forRanks,EER(i)] = BASE.getBiometrics(TScores,FScores,ones(1,nrT),ones(size(ones(1,nrT),1),1));
                forRanks = (forRanks./(nrF+1)).*100;
                CumRank = zeros(1,100);
                for cr=1:1:100
                    CumRank(1,cr) = (sum(forRanks<=cr)/(nrT)).*100;
                end
                CumRanks(:,i) = CumRank';
                parfor_progress;
            end
            closeParForProgress(path,ID);
            out.EER = EER;
            out.CumRanks = CumRanks;
            out.ClassACC = ClassACC;
            obj.GBBioTest = out;
       end
   end
   methods % INTERFACING FUNCTIONS
       function obj = reduceSamples(obj,index)
          if nargout==1, obj = clone(obj); end
          obj.ShapeIndex = obj.ShapeIndex(index);
          obj.COV = obj.COV(index,:);
          obj.GB = obj.GB(index,:);
          obj.COVMOD = [];
          obj.GBMOD = [];
       end
       function out = getShapeDM(obj)
           if isempty(obj.ShapeSpace),out = [];return;end
           out = squareform(pdist(obj.DepVar,'euclidean'));% Shape Similarity Matrix
       end
       function out = getGBDM(obj)
           if isempty(obj.GB),out = [];return;end
           out = squareform(pdist(obj.GB,'euclidean'));
       end
       function out = reduceGB(obj,index)
           if nargout==1, obj = clone(obj);end
           obj.GB = obj.GB(:,index);
           obj.GBindex = obj.GBindex(index);
           if ~isempty(obj.GBMOD), obj.GBMOD = obj.GBMOD(index);end
           out = obj;
       end
       function out = build(obj,DepVar,IndVar,BootInd,CondInd,PartInd)
           switch obj.BuildMethod
               case 'PLSR'
                   out = PLSR(obj,DepVar,IndVar,BootInd,CondInd);
               case 'BRIM'
                   out = BRIM(obj,DepVar,IndVar,BootInd,CondInd,PartInd);
               case 'FOLDEDBRIM'
                   out = FOLDEDBRIM(obj,DepVar,IndVar,BootInd,CondInd,PartInd);
           end
       end
       function out = PLSR(obj,DepVar,IndVar,BootInd,CondInd)                  
                % Building Independent and Dependent Variables
                  nr2Condition = length(CondInd);
                  var.nr2Boot = length(BootInd);
                  var.Booting = (nr2Condition+1:nr2Condition+var.nr2Boot);
                % Output
                  out.IndVar = IndVar(:,var.Booting);out.DepVar = DepVar;
                % examining variables to Boot and create reference faces, the lattter is
                % needed to rescale RIP variables onto common scales
                  for i=1:1:var.nr2Boot
                     v = IndVar(:,var.Booting(i));
                     var.Info{i} = BASE.getVarInfo(v,DepVar);    
                  end
                  out.Var = var;
                  % getting the final M from outer training data
                  [M,H] = BASE.getBootSampleRegression(IndVar,DepVar,var,obj.RegRuns,obj.SamplePerc,obj.Htestp);
                  if obj.Htest, M = H.*M;end
                  out.M = M;out.H = H;
                  RIP = BASE.updateRIP(DepVar,M)';
                  if obj.RIPNorm, RIP = BASE.normalizeRIP(RIP,M,var);end
                  out.RIP = RIP;
                  [out.RIPStat] = BASE.performRIPStat(out.IndVar,RIP,var,0);
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
                     var.Info{i} = BASE.getVarInfo(v,DepVar);    
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
                                [M,H] = BASE.getBootSampleRegression(IndVar(FiTrInd,:),DepVar(FiTrInd,:),var,obj.RegRuns,obj.SamplePercentage,obj.HtestP);
                                if obj.Htest, M = H.*M;end% setting close to zero partial regression coefficients to zero
                                rip = BASE.updateRIP(DepVar(FiTestInd,:),M)';
                                if obj.RIPNorm, rip = BASE.normalizeRIP(rip,M,var);end
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
                  [M,H] = BASE.getBootSampleRegression(IndVar,DepVar,var,obj.RegRuns,obj.SamplePercentage,obj.HtestP);
                  if obj.Htest, M = H.*M;end
                  out.M = M;out.H = H;
                  RIP = BASE.updateRIP(DepVar,M)';
                  if obj.RIPNorm, RIP = BASE.normalizeRIP(RIP,M,var);end
                  out.RIP = RIP;
                  %[out.RIPStat] = BASE.performRIPStat(out.IndVar,RIP,var,0);
       end
       function illustrateGB(obj,index)
           disp(['Axis Name: ' obj.GBNames{index}]);
           mod = obj.GBMOD{index};
           %disp(['Correlation: ' num2str(mod.RIPStat)]);
           figure;plot(mod.IndVar,getMODRIP(obj,mod),'b.');
           C = mod.M.*obj.ShapeSpace.EigStd';
           vec = obj.ShapeSpace.EigVec*C';
           vec = reshape(vec,3,length(vec)/3);
           vec = sqrt(sum(vec.^2,1));
           scan = clone(obj.ShapeSpace.Average);
           scan.Value = vec;
           scan.ColorMode = 'Indexed';
           viewer(scan);
       end
       function [MCOR,GBCOR] = illustrateGBCOR(obj)
           MCOR = nan*zeros(obj.nGB,obj.nGB);
           GBCOR = nan*zeros(obj.nGB,obj.nGB);
           for i=1:1:obj.nGB
               for j=1:1:obj.nGB
                   MCOR(i,j) = angle(obj.GBMOD{i}.M',obj.GBMOD{j}.M');
                   tmpc = corrcoef(obj.GB(:,i),obj.GB(:,j));
                   GBCOR(i,j) = tmpc(1,2);
               end
           end
           figure;imagesc(abs(MCOR));set(gca,'clim',[0 1]);
           figure;imagesc(abs(GBCOR));set(gca,'clim',[0 1]);
       end    
       function out = getPartitionsCOVR2(obj,F)
           [~,t] = size(F);
           K = length(unique(F(:,1)));
           out = nan*zeros(K,t);
           parfor i=1:t
              tmp = zeros(K,1);
              for k=1:K
                  ind = find(F(:,i)==k);
                  [~,tmp(k)] = BASE.getMR2(obj.COV(ind,:),obj.DepVar(ind,:));
              end
              out(:,i) = tmp;
           end
       end
       function out = getMODRIP(obj,MOD)
          if obj.HTest, 
             if obj.RIPNorm, out = MOD.NormHMRIP; return; end
             out = MOD.HMRIP; return;
          end
          if obj.RIPNorm, out = MOD.NormMRIP; return; end
          out = MOD.MRIP;
       end
       function out = getMODM(obj,MOD)
           if obj.HTest, out = MOD.HM; return; end
           out = MOD.M;
       end
       function [A,M,STAT] = AngleTest(obj,IndVar,DepVar,Aval,maxM,t)
               A = nan*zeros(t,1);
               M = nan*zeros(obj.ShapeDim,2,maxM);
               [nS,vid] = size(IndVar);
               el = 0;
               for i=1:t
                   el = el+1;
                   F = crossvalind('Kfold',nS,2);
                   forM = nan*zeros(obj.ShapeDim,2);
                   for k=1:1:2
                       ind = find(F==k);
                       forM(:,k) = BASE.getRegression(IndVar(ind,:),DepVar(ind,:),vid);
                   end
                   A(i) = angle(forM(:,1),forM(:,2));
                   M(:,:,el) = forM;
                   if mod(i,maxM)==0, el = 0; end
                   if mod(i,100)==0
                      tmppA = (length(find(A(1:i)<=0))+1)/(i+1); 
                      acc = 10/i;
                      if tmppA>acc, break;end
                   end
               end
               STAT.Aval = Aval;
               STAT.pA = zeros(length(Aval),1);
               for j=1:1:length(Aval)
                   STAT.pA(j) = (length(find(A(1:i)<=STAT.Aval(j)))+1)/(i+1);  
               end
               STAT.avgA = mean(A(1:i));
               STAT.medA = median(A(1:i));
               STAT.stdA = std(A(1:i));
               STAT.madA = mad(A(1:i));
               sortA = sort(A(1:i));
               STAT.upperciA = sortA(round(0.975*i));
               STAT.lowerciA = sortA(round(0.025*i+1));    
       end
   end
   methods % BIOMETRICS
       function out = biometrics(obj,TestDepVar,TestCOV,TestGB,indCOV,short,vis)
           % TO BE UPDATED
           % Initialize
            nrT = size(TestDepVar,1);
           % extracting COV RIP values
            RIPCOV = nan*zeros(nrT,obj.nCOV);
            wCOVA = nan*zeros(obj.nCOV,1);
            wCOVP = nan*zeros(obj.nCOV,1);
            for i=1:1:obj.nCOV
                tmp = BASE.updateRIP(TestDepVar,getMODM(obj,obj.COVMOD{i}))';
                if obj.RIPNorm, tmp = BASE.normalizeRIP(tmp,getMODM(obj,obj.COVMOD{i}),obj.COVMOD{i}.Var);end
                RIPCOV(:,i) = tmp;
                wCOVA(i) = max(obj.COVATest.AvgA(i),0);
                wCOVP(i) = -log(obj.COVATest.pA(1,i));
                %figure;plot(TestCOV(:,i),tmp,'b.');
            end 
           % extracting GB RIP values
            RIPGB = nan*zeros(nrT,obj.nGB);
            wGBA = nan*zeros(obj.nGB,1);
            wGBP = nan*zeros(obj.nGB,1);
            for i=1:1:obj.nGB
                tmp = BASE.updateRIP(TestDepVar,getMODM(obj,obj.GBMOD{i}))';
                if obj.RIPNorm, tmp = BASE.normalizeRIP(tmp,getMODM(obj,obj.GBMOD{i}),obj.GBMOD{i}.Var);end
                RIPGB(:,i) = tmp;
                wGBA(i) = max(obj.GBATest.AvgA(i),0);
                wGBP(i) = -log(obj.GBATest.pA(1,i));
                %figure;plot(TestGB(:,i),tmp,'b.');
            end
           % extracting matching scores 
            TMatchesGB = nan*zeros(obj.nGB,nrT);
            FMatchesGB = nan*zeros(obj.nGB,nrT,nrT-1);
            TMatchesCOV = nan*zeros(obj.nCOV,nrT);
            FMatchesCOV = nan*zeros(obj.nCOV,nrT,nrT-1);
            parfor t=1:nrT
                % t=1;
                tInd = setdiff(1:nrT,t);
                % GENETIC BACKGROUND MATCHES
                tmpFMatches = nan*zeros(obj.nGB,nrT-1);
                tmpTMatches = nan*zeros(obj.nGB,1);
                for i=1:1:obj.nGB
                    Distr = BASE.X2RIP(obj.GBMOD{i}.IndVar,getMODRIP(obj,obj.GBMOD{i}),TestGB(t,i),obj.GBMOD{i}.Var.Info{1}.Type,obj.Kappa,obj.KN);
                    m = getMatch(obj,Distr,RIPGB(:,i));
                    tmpFMatches(i,:) = m(tInd);
                    tmpTMatches(i,:) = m(t);
                end
                TMatchesGB(:,t) = tmpTMatches;
                FMatchesGB(:,t,:) = tmpFMatches;
                % COVARIATE MATCHES
                tmpFMatches = nan*zeros(obj.nCOV,nrT-1);
                tmpTMatches = nan*zeros(obj.nCOV,1);
                for i=1:1:obj.nCOV
                    Distr = BASE.X2RIP(obj.COVMOD{i}.IndVar,getMODRIP(obj,obj.COVMOD{i}),TestCOV(t,i),obj.COVMOD{i}.Var.Info{1}.Type,obj.Kappa,obj.KN);
                    m = getMatch(obj,Distr,RIPCOV(:,i));
                    tmpFMatches(i,:) = m(tInd);
                    tmpTMatches(i,:) = m(t);
                end
                TMatchesCOV(:,t) = tmpTMatches;
                FMatchesCOV(:,t,:) = tmpFMatches;
            end
            switch obj.Classify
               case 'SOFT'
                   TScoresGB = match2score(obj,TMatchesGB);out.TScoresGB = TScoresGB;
                   FScoresGB = match2score(obj,FMatchesGB);out.FScoresGB = FScoresGB;
                   TScoresCOV = match2score(obj,TMatchesCOV);out.TScoresCOV = TScoresCOV;
                   FScoresCOV = match2score(obj,FMatchesCOV);out.FScoresCOV = FScoresCOV;
               case 'HARD'
                   TScoresGB = BASE.class2score(match2class(obj,TMatchesGB));out.TScoresGB = TScoresGB;
                   FScoresGB = BASE.class2score(match2class(obj,FMatchesGB));out.FScoresGB = FScoresGB;
                   TScoresCOV = BASE.class2score(match2class(obj,TMatchesCOV));out.TScoresCOV = TScoresCOV;
                   FScoresCOV = BASE.class2score(match2class(obj,FMatchesCOV));out.FScoresCOV = FScoresCOV;
            end 
            out.ClassAccGB = ((sum(match2class(obj,TMatchesGB),2)/nrT)*100)';
            out.ClassAccCOV = ((sum(match2class(obj,TMatchesCOV),2)/nrT)*100)';
            out.wCOVA = wCOVA;out.wCOVP = wCOVP;
            out.wGBA = wGBA;out.wGBP = wGBP;
            
            if short, return; end
            %indCOV = 1;% SELECTION SEX ONLY
            %indCOV = 1:4;
            nCOV = length(indCOV);
            
            wCOV = wCOVA(indCOV);
            wGB = wGBA;
            
            %wCOV = wCOVP(indCOV);
            %wGB = wGBP;
            
            
            nr = 6;
            EER = nan*zeros(1,nr);Ranks = nan*zeros(1,nr,nrT);
            X = nan*zeros(nr,(nrT*nrT)-1);Y = nan*zeros(nr,(nrT*nrT)-1);
            r = 0;
            r=r+1;% COVARIATES UNWEIGHTED
            [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BASE.getBiometrics(TScoresCOV(indCOV,:),FScoresCOV(indCOV,:,:),ones(nCOV,nrT),ones(nCOV,1));
            r=r+1;% COVARIATES WEIGHTED
            [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BASE.getBiometrics(TScoresCOV(indCOV,:),FScoresCOV(indCOV,:,:),ones(nCOV,nrT),wCOV);
            r=r+1;% GENETIC BACKGROUND UNWEIGHTED
            [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BASE.getBiometrics(TScoresGB,FScoresGB,ones(obj.nGB,nrT),ones(obj.nGB,1));
            r=r+1;% GENETIC BACKGROUND WEIGHTED
            [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BASE.getBiometrics(TScoresGB,FScoresGB,ones(obj.nGB,nrT),wGB);
            r=r+1;% COMBINED UNWEIGHTED
            [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BASE.getBiometrics([TScoresCOV(indCOV,:); TScoresGB],[FScoresCOV(indCOV,:,:); FScoresGB],ones(obj.nGB+nCOV,nrT),ones(nCOV+obj.nGB,1));
             r=r+1;% COMBINED WEIGHTED
            [Ranks(:,r,:),EER(:,r),~,~,X(r,:),Y(r,:)] = BASE.getBiometrics([TScoresCOV(indCOV,:); TScoresGB],[FScoresCOV(indCOV,:,:); FScoresGB],ones(obj.nGB+nCOV,nrT),[wCOV;wGB]);
            Ranks(1,:,:) = (Ranks(1,:,:)./nrT).*100;
            CumRanks = zeros(1,nr,100);
            for cr=1:1:100
                CumRanks(:,:,cr) = (sum(Ranks<=cr,3)/nrT).*100;
            end
            out.EER = EER;
            out.CumRanks = CumRanks;
            if vis
               str = {'g-.' 'g-' 'b-.' 'b-' 'k-.' 'k-'};
               figure;set(gca,'ylim',[0 100],'xlim',[1 100]);hold on;
               xlabel('Rank (%)');ylabel('Identified (%)' );grid on;
               title('Identification BASE');
               plot(1:100,1:100,'k--');
               for i=1:1:nr
                   plot(1:100,squeeze(CumRanks(1,i,:)),str{i},'LineWidth',1);
               end
               figure;hold on;title('Verification BASE');
               xlabel('false positive fraction');ylabel('true positive fraction' );grid on;
               plot(0:0.01:1,0:0.01:1,'k--');
               plot(0:0.01:1,1:-0.01:0,'k-');
               for i=1:1:nr
                   plot(X(i,:),Y(i,:),str{i},'LineWidth',1);
               end
            end      
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
       function out = match2score(obj,matches)
          switch obj.Match
              case 'QD'
                  out = -log(matches);
                  %out = 1./matches;
              case 'INLIER'
                  out = 1-matches;
          end
       end
       function out = match2class(obj,matches)
          switch obj.Match
              case 'QD'
                  T = 1;
              case 'INLIER'
                  T = 0.5;
          end
          out = BASE.matchthresholding(matches,T);
       end
   end
   methods (Static = true)
       function [R2,pR2,M] = testPartialPLSR(X,Y,C,t)
           index = intersect(BASE.notNAN(C),BASE.notNAN(X));
           if ~isempty(C)
            E = BASE.getResiduals(C(index,:),Y(index,:));
            X = BASE.getResiduals(C(index,:),X(index,:));
           else
            E = Y(index,:);
            X = X(index,:);
           end
           [~,~,~,~,M,var] = plsregress(X,E,1);R2 = var(2);
           M = M(end,:);
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
       function out = getResiduals(X,Y)
                [~,~,~,~,~,~,~,stats] = plsregress(X,Y,min(size(X,2),size(Y,2)));
                out = stats.Yresiduals;
       end
       function out = notNAN(in)
           index = (1:size(in,1));
           [i,~] = find(isnan(in));
           out = setdiff(index,unique(i));
       end
       function [A,pA] = regressionAngle(M1,M2,t)
           dim = size(M1,2);
           A = angle(M1',M2');
           if t<1, pA = nan; return;end
           ACount = false(t,1);
           parfor p=1:t
               ind = randperm(dim);
               Afor = angle(M1(ind)',M2');%#ok<*PFBNS>
               ACount(p) = Afor>=A;
           end
           pA = (sum(ACount)+1)/(t+1);
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
       function [out] = updateRIP(in,M)
                    % in this implementation I take the reference as the origin
                    n2 = size(in,1); n1 = size(M,1);% determine input size
                    in = in'; out = nan*zeros(n1,n2);% allocate memory
                    for j=1:n1
                        out(j,:) = dot(in,repmat(M(j,:)'/norm(M(j,:)'),1,n2));
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
                    M(:,:,s) = BASE.getRegression(IndVar(SampleFi,:),DepVar(SampleFi,:),var.Booting);
                end
                if runs==1,M = squeeze(M);H = ones(size(M));return;end
                H = BASE.wilcoxonM(M,var,Tp);
                M = median(M,3);       
       end
       function M = getRegression(A,B,Booting)
                 if nargin < 3, Booting = (1:size(A,2));end
                 [A,B] = eliminateNAN(A,B);
                 [~,~,~,~,M] = plsregress(A,B,size(A,2));
                 M = M(1+Booting,:);
       end
       function [M,R2] = getMR2(A,B)
                [A,B] = eliminateNAN(A,B);
                [~,~,~,~,M,var] = plsregress(A,B,size(A,2));
                R2 = var(2);
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
       function H = wilcoxonMSimple(M,Tp)
            nrD = size(M,2);
            %H = zeros(1,nrD);
            P = zeros(1,nrD);
            for j=1:1:nrD
               %[P(j),H(j)] = signrank(M(:,j));% wilcoxon test for median  
               P(j) = signrank(M(:,j));% wilcoxon test for median  
            end
            H = P<=Tp;
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
       function [Distr,RIP] = X2RIP(trX,trrip,X,type,kappa,K)
           if nargin<5, kappa = 6; end
           if nargin<6, K = 300; end
           switch type
               case 'Categorical'
                  Tind = find(trX==X);% Inlier Distribution
                  Find = find(trX==-1*X);% Outlier Distribution
                  [Distr.Tmu,Distr.Tsigma,Distr.TO] = BASE.extractRobustRIPStatistics(trrip(Tind)',kappa);
                  RIP = Distr.Tmu;Distr.TpMax = normpdf(Distr.Tmu,Distr.Tmu,Distr.Tsigma);
                  [Distr.Fmu,Distr.Fsigma,Distr.FO] = BASE.extractRobustRIPStatistics(trrip(Find)',kappa);
                  Distr.FpMax = normpdf(Distr.Fmu,Distr.Fmu,Distr.Fsigma);      
               case 'Continuous'
                  ind = find(~isnan(trX));% eliminate nans as wrong data entries (typically coded as 0);
                  trX = trX(ind);
                  W = abs(trX-X);
                  % select closest samples
                  [~,ind2] = sort(W,'ascend');
                  % TRUE DISTRIBUTION
                  ind2T = ind2(1:K);indT = ind(ind2T);
                  rip = trrip(indT)';
                  Distr.IndT = indT;
                  WT = W(ind2T);WT = WT-min(WT);WT = WT/max(WT);WT = 1-WT';
                  [Distr.Tmu,Distr.Tsigma,Distr.TO] = BASE.extractRobustRIPStatistics(rip,kappa,WT);
                  Distr.TpMax = normpdf(Distr.Tmu,Distr.Tmu,Distr.Tsigma);
                  RIP = Distr.Tmu;
                  % FALSE DISTRIBUTION
                  ind2F = ind2(end-K+1:end);indF = ind(ind2F);
                  rip = trrip(indF)';
                  Distr.IndF = indF;
                  WF = W(ind2F);WF = WF-min(WF);WF = WF/max(WF);WF = WF';
                  [Distr.Fmu,Distr.Fsigma,Distr.FO] = BASE.extractRobustRIPStatistics(rip,kappa,WF);
                  Distr.FpMax = normpdf(Distr.Fmu,Distr.Fmu,Distr.Fsigma); 
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
                           [stat(s),~,statP(s)] = BASE.myPermAnova(XG,vrip(index),t);
                       case 'Continuous'
                           statType{s} = 'Cor';
                           [stat(s),statP(s)] = BASE.permCorr(v(index),vrip(index),t);
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
       function out = matchthresholding(matches,T)
                out = matches>=T; 
       end
       function out = class2score(classes)
                out = 1-classes; 
       end
       function [R,EER,ST,SF,X,Y] = getBiometrics(TScores,FScores,wM,W)
                 nF = size(FScores,3);
                 nT = size(TScores,2);
                 tmpW = repmat(W,1,nT).*wM;
                 ST = BASE.accScores(TScores,tmpW)';
                 SF = BASE.accScores(FScores,repmat(tmpW,1,1,nF));
                 [R,EER,X,Y] = BASE.evalScores(ST,SF,nT,nF); 
       end
       function Score = accScores(Scores,tmpW)
                Score = squeeze(nansum(tmpW.*Scores,1)./nansum(tmpW,1));
       end
       function [R,EER,x,y] = evalScores(ST,SF,nT,nF)
                 % Rank analysis
                 R = (sum(SF<repmat(ST,1,nF),2)+1);%/(nF+1);
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
   end
end

% function out = runAngleTestGB(obj,t,maxM)
%            if nargin<3, maxM = t; end
%            K=2;
%            M = nan*zeros(obj.ShapeDim,2,maxM,obj.nGB);
%            Aval = 0:0.025:0.4;
%            A = nan*zeros(t,obj.nGB);
%            %pA = nan*zeros(1,obj.nGB);
%            pA = nan*zeros(length(Aval),obj.nGB);
%            avgA = nan*zeros(1,obj.nGB);
%            stdA = nan*zeros(1,obj.nGB);
%            medA = nan*zeros(1,obj.nGB);
%            madA = nan*zeros(1,obj.nGB);
%            lowerciA = nan*zeros(1,obj.nGB);
%            upperciA = nan*zeros(1,obj.nGB);
%            vid = obj.nCOV+1;
%            [path,ID] = setupParForProgress(obj.nGB);
%            parfor a=1:obj.nGB
%                tmpA = nan*zeros(t,1);
%                tmpM = nan*zeros(obj.ShapeDim,2,maxM);
%                el = 0;
%                for i=1:t
%                    el = el+1;
%                    F = crossvalind('Kfold',obj.nS,2);
%                    forM = nan*zeros(obj.ShapeDim,K);
%                    for k=1:1:K
%                        ind = find(F==k);
%                        forM(:,k) = BASE.getRegression([obj.COV(ind,:) obj.GB(ind,a)],obj.DepVar(ind,:),vid);
%                    end
%                    tmpA(i) = angle(forM(:,1),forM(:,2));
%                    tmpM(:,:,el) = forM;
%                    if mod(i,maxM)==0, el = 0; end
%                    if mod(i,100)==0
%                       tmppAB = (length(find(tmpA(1:i)<=0))+1)/(i+1); 
%                       acc = 10/i;
%                       if tmppAB>acc, break;end
%                    end
%                end
%                A(:,a) = tmpA;
%                M(:,:,:,a) = tmpM;
%                %pA(a) = (length(find(tmpA(1:i)<=0))+1)/(i+1);
%                tmppA = zeros(length(Aval),1);
%                for j=1:1:length(Aval)
%                    tmppA(j) = (length(find(tmpA(1:i)<=Aval(j)))+1)/(i+1);  
%                end
%                pA(:,a) = tmppA;
%                avgA(a) = mean(tmpA(1:i));
%                medA(a) = median(tmpA(1:i));
%                stdA(a) = std(tmpA(1:i));
%                madA(a) = mad(tmpA(1:i));
%                sortA = sort(tmpA(1:i));
%                upperciA(a) = sortA(round(0.975*i));
%                lowerciA(a) = sortA(round(0.025*i+1));
%                parfor_progress;
%            end
%            closeParForProgress(path,ID);
%            out.A = A;out.M = M;
%            out.AvgA = avgA;out.MedA = medA;out.StdA = stdA;out.MadA = madA;out.LowA = lowerciA;out.UpA = upperciA;out.pA = pA;
%            out.Aval = Aval;
%            obj.GBATest = out;
%        end

% function [out,A,R2] = partitionTest(obj,t)
%            R2 = nan*zeros(t,2,obj.nGB);
%            A = nan*zeros(t,obj.nGB);
%            F = DOB_SCV_DM(2,obj.nS,[]);
%            [path,ID] = setupParForProgress(t);
%            parfor i=1:t
%                ind = randperm(obj.nS);
%                Ffor = F(ind);
%                M = nan*zeros(2,obj.nGB,obj.ShapeDim);
%                tmpR2 = nan*zeros(2,obj.nGB);
%                for f=1:1:2
%                   SInd = find(Ffor==f);
%                   for a=1:obj.nGB
%                     [tmpR2(f,a),~,M(f,a,:)] = BASE.testPartialPLSR(obj.GB(SInd,a),obj.DepVar(SInd,:),obj.COV(SInd,:),0);
%                   end
%                end
%                tmpA = nan*zeros(1,obj.nGB);
%                for a=1:1:obj.nGB
%                    tmpA(a) = angle(squeeze(M(1,a,:)),squeeze(M(2,a,:)));
%                end
%                A(i,:) = tmpA;
%                R2(i,:,:) = tmpR2;
%                parfor_progress;
%            end
%            closeParForProgress(path,ID);
%            R2 = [squeeze(R2(:,1,:));squeeze(R2(:,2,:))];
%            pA = nan*zeros(1,obj.nGB);
%            pAB = nan*zeros(1,obj.nGB);
%            avgA = nan*zeros(1,obj.nGB);
%            stdA = nan*zeros(1,obj.nGB);
%            lowciA = nan*zeros(1,obj.nGB);
%            upperciA = nan*zeros(1,obj.nGB);
%            avgR2 = nan*zeros(1,obj.nGB);
%            stdR2 = nan*zeros(1,obj.nGB);
%            lowciR2 = nan*zeros(1,obj.nGB);
%            upperciR2 = nan*zeros(1,obj.nGB);
%            for i=1:obj.nGB
%               pA(i) = signrank(A(:,i));
%               pAB(i) = (length(find(A(:,i)<=0))+1)/(t+1);
%               avgA(i) = mean(A(:,i));
%               stdA(i) = std(A(:,i));
%               sortA = sort(A(:,i));
%               lowciA(i) = sortA(round(0.975*t));
%               upperciA(i) = sortA(round(0.025*t+1));
%               avgR2(i) = mean(R2(:,i));
%               stdR2(i) = std(R2(:,i));
%               sortR2 = sort(R2(:,i));
%               lowciR2(i) = sortR2(round(0.975*t));
%               upperciR2(i) = sortR2(round(0.025*t+1));
%            end
%            out.AvgA = avgA;out.StdA = stdA;out.LowA = lowciA;out.UpA = upperciA;out.pA = pA;out.pAB = pAB;
%            out.AvgR2 = avgR2;out.StdR2 = stdR2;out.LowR2 = lowciR2;out.UpR2 = upperciR2;
%        end

