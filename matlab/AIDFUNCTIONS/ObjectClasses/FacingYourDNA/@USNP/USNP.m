classdef USNP < superClassSNP
   % General Properties
   properties
      RS = [];
      GT = [];
   end
   properties (Dependent = true)
      GG;
      GGGTBal;
      AL;
      mAl;
      mAlF;
      MAl;
      MAlF;
      nm;
      nM;
      nmm;
      nMm;
      nMM;
      fmm;
      fMm;
      fMM;
      Label;
      nS;
   end
   properties
      TestInd = 0;
      BuildMethod = 'PLSR';
      RIPNorm = true;
      HTest = true;
      MOD = [];
   end
   properties (Dependent = true)
      MODGG;
      MODM;
      MODRIP;
   end
   properties
       Match = 'QD';
       Classify = 'SOFT';
   end
   properties
       R2 = 0;
       pR2 = 1;
       A = 0;
       pA = 1;
       TestBio = [];
   end
   properties (Hidden = true, Dependent = true)
      Balances;
      mmGT;
      MmGT;
      MMGT;
      mmind;
      Mmind;
      MMind;
      nAA;
      nAB;
      nBB;
      AAind;
      ABind;
      BBind;
      fAA;
      fAB;
      fBB;
      nA;
      nB;
      fA;
      fB;
      nAl;
      nGT;
   end
   properties (Hidden = true)
      RegRuns = 100; % Number of Regression runs (SNPBRIM)
      SamplePercentage = 1;
      SampleWR = true;
      RegSampling = true;% Bootstrap/None (SNPBRIM)
      RegBalance = true; % balance the data yes or no (SNPBRIM)
      RegBalanceTH = 0.4; % desired balance factor (SNPBRIM)
      RegBalanceMethod = 'UpSample'; % UpSample/DownSample/ADASYN (SNPBRIM)
      OuterFold = 8; % number of Outer Folds (SNPBRIM)
      InnerFold = 10; % number of Inner Folds (SNPBRIM)
      UseRedShape = false; % Use reduced space or covariates (SNPBRIM)
      MaxIterations = 5; % Maximum number of BRIM iterations (SNPBRIM)
      StopCorr = 0.98; % Stopping correlation between subsequent iterations (SNPBRIM)
      RIPStatPerm = 1000; % Number of permutations in testing significance of RIP values (SNPBRIM)
      GroupDef = 'GT';% Definition of groups, GG = according to Test performed, GT = according to orginal genotypes (SNPBRIM)
      RedShape = false;
   end
   methods % CONSTRUCTOR
        function obj = USNP(varargin)
            obj = obj@superClassSNP(varargin{:});         
        end
   end
   methods % GETTING/SETTING
       function out = get.GG(obj)
           if isempty(obj.GT), out = []; return; end
           out = GT2GG(obj,obj.GT);
       end
       function out = get.GGGTBal(obj)
          if isempty(obj.GT), out = []; return; end
          out = USNP.getGGGTBalance(obj.GG,obj.GT);
       end
       function out = get.AL(obj)
           if isempty(obj.GT), out = [];return;end
           out = nan*obj.GT;
           out(obj.ABind) = 1;
           switch obj.mAl
               case 'A'
                  out(obj.AAind) = 2;
                  out(obj.BBind) = 0;
               case 'B'
                  out(obj.AAind) = 0;
                  out(obj.BBind) = 2;
           end
       end
       function out = get.AAind(obj)
           if isempty(obj.GT), out = []; return; end
           out = find(obj.GT==-1);
       end
       function out = get.ABind(obj)
           if isempty(obj.GT), out = []; return; end
           out = find(obj.GT==0);
       end
       function out = get.BBind(obj)
           if isempty(obj.GT), out = []; return; end
           out = find(obj.GT==1);
       end
       function out = get.nAA(obj)
           if isempty(obj.GT), out = 0; return; end
           out = length(obj.AAind);
       end
       function out = get.nAB(obj)
           if isempty(obj.GT), out = 0; return; end
           out = length(obj.ABind);
       end
       function out = get.nBB(obj)
           if isempty(obj.GT), out = 0; return; end
           out = length(obj.BBind);
       end
       function out = get.nA(obj)
                out = 2*obj.nAA+obj.nAB; 
       end
       function out = get.nB(obj)
                out = 2*obj.nBB+obj.nAB;
       end
       function out = get.nm(obj)
           switch obj.mAl
               case 'A'
                   out = obj.nA;
               case 'B'
                   out = obj.nB;
           end
       end
       function out = get.nM(obj)
           switch obj.mAl
               case 'A'
                   out = obj.nB;
               case 'B'
                   out = obj.nA;
           end
       end
       function out = get.fAA(obj)
                out = sqrt(obj.nAA/(obj.nAA+obj.nAB+obj.nBB));
       end
       function out = get.fAB(obj)
                out = sqrt(obj.nAB/(obj.nAA+obj.nAB+obj.nBB));
       end
       function out = get.fBB(obj)
                out = sqrt(obj.nBB/(obj.nAA+obj.nAB+obj.nBB));
       end
       function out = get.nGT(obj)
           out = length(obj.GT);
       end
       function out = get.nAl(obj)
           out = obj.nGT*2;
       end
       function out = get.fA(obj)
           if isempty(obj.GT), out = 0; return; end
           out = obj.nA/obj.nAl;
       end
       function out = get.fB(obj)
           if isempty(obj.GT), out = 0; return; end
           out = obj.nB/obj.nAl;
       end
       function out = get.mAl(obj)
           if obj.nA<=obj.nB
               out = 'A';
           else
               out = 'B';
           end
       end
       function out = get.MAl(obj)
           if obj.nA>=obj.nB
               out = 'A';
           else
               out = 'B';
           end
       end
       function out = get.Balances(obj)
           out = zeros(1,6);
           if isempty(obj.GT), return; end
           out(1) = USNP.getBalance(obj.nAA+obj.nAB,obj.nBB);% Recessive
           out(2) = USNP.getBalance(obj.nBB+obj.nAB,obj.nAA);% Dominant
           out(3) = USNP.getBalance(obj.nAA+obj.nBB,obj.nAB);% Over Dominant
           out(4) = USNP.getBalance(obj.nAA,obj.nBB);% AA versus BB
           out(5) = USNP.getBalance(obj.nAA,obj.nAB);% AA versus AB
           out(6) = USNP.getBalance(obj.nAB,obj.nBB);% AB versus BB
       end
       function out = get.Label(obj)
           switch obj.TestInd
               case 0
                   out = 'NO EFFECT';
               case 1
                   out = 'RECESSIVE';
               case 2
                   out = 'DOMINANT';
               case 3
                   out = 'OVER DOMINANT';
               case 4
                   out = 'ADDITIVE';
               case 5
                   out = 'AA BB';
               case 6
                   out = 'AA AB';
               case 7
                   out = 'AB BB';
               case 8
                   out = 'NON ADDITIVE';
               case 9
                   out = 'ADDITIVE AA AB';
               case 10
                   out = 'ADDITIEVE AB BB';    
           end
       end
       function out = get.mAlF(obj)
           out = min([obj.fA obj.fB]); 
       end
       function out = get.MAlF(obj)
           out = max([obj.fA obj.fB]);
       end
       function out = get.mmind(obj)
          if isempty(obj.AL), out = [];return;end 
          out = find(obj.AL==2); 
       end
       function out = get.Mmind(obj)
          if isempty(obj.AL), out = [];return;end
          out = find(obj.AL==1); 
       end
       function out = get.MMind(obj)
          if isempty(obj.AL), out = [];return;end 
          out = find(obj.AL==0); 
       end
       function out = get.nmm(obj)
           out = length(obj.mmind);
       end
       function out = get.nMm(obj)
           out = length(obj.Mmind);
       end
       function out = get.nMM(obj)
           out = length(obj.MMind);
       end
       function out = get.fmm(obj)
                out = sqrt(obj.nmm/(obj.nmm+obj.nMm+obj.nMM));
       end
       function out = get.fMm(obj)
                out = sqrt(obj.nMm/(obj.nmm+obj.nMm+obj.nMM));
       end
       function out = get.fMM(obj)
                out = sqrt(obj.nMM/(obj.nmm+obj.nMm+obj.nMM));
       end
       function out = get.mmGT(obj)
                switch obj.mAl
                    case 'A'
                        out = -1;
                    case 'B'
                        out = 1;
                end 
       end
       function out = get.MmGT(obj) %#ok<*MANU>
                out = 0; 
       end
       function out = get.MMGT(obj)
                switch obj.mAl
                    case 'A'
                        out = 1;
                    case 'B'
                        out = -1;
                end 
       end
       function out = get.nS(obj)
          out = length(obj.GT); 
       end
       function out = get.MODM(obj)
          if isempty(obj.MOD), out = []; return; end
          if obj.HTest, out = obj.MOD.HM; return; end
          out = obj.MOD.M;
       end
       function out = get.MODRIP(obj)
          if isempty(obj.MOD), out = []; return; end
          if obj.HTest, 
             if obj.RIPNorm, out = obj.MOD.NormHMRIP; return; end
             out = obj.MOD.HMRIP; return;
          end
          if obj.RIPNorm, out = obj.MOD.NormMRIP; return; end
          out = obj.MOD.MRIP;
       end
       function out = get.MODGG(obj)
           if isempty(obj.MOD), out = []; return; end
           out = obj.MOD.GG;
       end
   end
   methods % DATA CHECKING
       function out = dataCheckGT(obj,MinSampleSize)
           if nargin<2, MinSampleSize = 10;end
           mmOk = obj.nmm>=MinSampleSize;
           MmOk = obj.nMm>=MinSampleSize;
           MMOk = obj.nMM>=MinSampleSize;
           switch obj.TestInd
               case 0
                   out = false;
               case {1 2 3 4} % All three groups required
                   out = mmOk&MmOk&MMOk;
               case 5 % test 4: AA BB
                   out = mmOk&MMOk;
               case 6 % test 5: AA AB
                   out = mmOk&MmOk;
               case 7 % test 6: AB BB
                   out = MmOk&MMOk;
               otherwise
                   error('Wrong Test Index');
           end
       end
       function out = dataCheckGTFreq(obj,MinSampleFreq)
           if nargin<2, MinSampleFreq = 0.1;end
           mmOk = obj.fmm>=MinSampleFreq;
           MmOk = obj.fMm>=MinSampleFreq;
           MMOk = obj.fMM>=MinSampleFreq;
           switch obj.TestInd
               case 0
                   out = false;
               case {1 2 3 4} % All three groups required
                   out = mmOk&MmOk&MMOk;
               case 5 % test 4: AA BB
                   out = mmOk&MMOk;
               case 6 % test 5: AA AB
                   out = mmOk&MmOk;
               case 7 % test 6: AB BB
                   out = MmOk&MMOk;
               otherwise
                   error('Wrong Test Index');
           end
       end
       function out = dataCheckmAlF(obj,MinAlFreq)
           if nargin<2, MinAlFreq = 0.1;end
           out = obj.mAlF>=MinAlFreq;
       end
   end
   methods % SNP ANALYSIS
       function out = runR2Test(obj,DepVar,COV,t)
           [out.R2,out.pR2] = USNP.testPartialPLSR(GT2GG(obj,obj.GT),DepVar,COV,t);
           obj.R2 = out.R2;obj.pR2 = out.pR2;
       end
       function out = runAngleTest(obj,DepVar,COV,t,maxM)
           if nargin < 5, maxM = t; end
           K=2;ShapeDim = size(DepVar,2);
           M = nan*zeros(ShapeDim,2,maxM);% only keep track of maxM regressions for further use
           A = nan*zeros(t,1);
           GG = GT2GG(obj,obj.GT);
           el = 0;
           for i=1:t
               el = el+1;
               F = crossvalind('Kfold',obj.nS,2);
               forM = nan*zeros(ShapeDim,2);
               for k=1:1:K
                   ind = find(F==k);
                   forM(:,k) = USNP.getRegression([COV(ind,:) GG(ind)],DepVar(ind,:));
               end
               A(i) = angle(forM(:,1),forM(:,2));
               M(:,:,el) = forM;
               if mod(i,maxM)==0, el = 0; end
               if mod(i,100)==0
                  tmppAB = (length(find(A(1:i)<=0))+1)/(i+1); 
                  acc = 10/i;
                  if tmppAB>acc, break;end
               end
           end
           Aval = 0:0.025:0.4;
           pA = zeros(length(Aval),1);
           for j=1:1:length(Aval)
               pA(j) = (length(find(A(1:i)<=Aval(j)))+1)/(i+1);  
           end
           %pA = (length(find(A(1:i)<=0))+1)/(i+1); %#ok<*PROP>
           avgA = mean(A(1:i));
           medA = median(A(1:i));
           stdA = std(A(1:i));
           madA = mad(A(1:i));
           sortA = sort(A(1:i));
           upperciA = sortA(round(0.975*i));
           lowerciA = sortA(round(0.025*i+1));
           if i<maxM,M = M(:,:,1:i);end
           A = A(1:i);
           out.A = A;out.M = M;
           out.AvgA = avgA;out.MedA = medA;out.StdA = stdA;out.MadA = madA;out.LowA = lowerciA;out.UpA = upperciA;out.pA = pA;
           obj.A = out.AvgA;
           obj.pA = out.pA(1);
       end
       function out = runAngleTest_DOB_SCV(obj,DepVar,COV,t,maxM)
           if nargin < 5, maxM = t; end
           K=2;ShapeDim = size(DepVar,2);
           M = nan*zeros(ShapeDim,2,maxM);% only keep track of maxM regressions for further use
           A = nan*zeros(t,1);
           GG = GT2GG(obj,obj.GT);
           el = 0;
           SDM = squareform(pdist(DepVar,'euclidean'));
           for i=1:t
               el = el+1;
               %F = crossvalind('Kfold',obj.nS,2);
               F = DOB_SCV_DM(2,SDM,obj.GT);
               forM = nan*zeros(ShapeDim,2);
               for k=1:1:K
                   ind = find(F==k);
                   forM(:,k) = USNP.getRegression([COV(ind,:) GG(ind)],DepVar(ind,:));
               end
               A(i) = angle(forM(:,1),forM(:,2));
               M(:,:,el) = forM;
               if mod(i,maxM)==0, el = 0; end
               if mod(i,100)==0
                  tmppAB = (length(find(A(1:i)<=0))+1)/(i+1); 
                  acc = 10/i;
                  if tmppAB>acc, break;end
               end
           end
           Aval = 0:0.025:0.4;
           pA = zeros(length(Aval),1);
           for j=1:1:length(Aval)
               pA(j) = (length(find(A(1:i)<=Aval(j)))+1)/(i+1);  
           end
           %pA = (length(find(A(1:i)<=0))+1)/(i+1); %#ok<*PROP>
           avgA = mean(A(1:i));
           medA = median(A(1:i));
           stdA = std(A(1:i));
           madA = mad(A(1:i));
           sortA = sort(A(1:i));
           upperciA = sortA(round(0.975*i));
           lowerciA = sortA(round(0.025*i+1));
           if i<maxM,M = M(:,:,1:i);end
           A = A(1:i);
           out.A = A;out.M = M;
           out.AvgA = avgA;out.MedA = medA;out.StdA = stdA;out.MadA = madA;out.LowA = lowerciA;out.UpA = upperciA;out.pA = pA;
           obj.A = out.AvgA;
           obj.pA = out.pA(1);
       end
       function out = runBiometricTest(obj,t,TestDepVar,TestGT)
                % Initializing
                 TrGG = obj.MODGG;TrRIP = obj.MODRIP;PosClass = obj.MOD.Var.PosClass;
                 if nargin<3
                    TestGG = TrGG;
                    TestRIP = TrRIP;
                 else
                    TestGG = GT2GG(obj,TestGT);
                    TestRIP = SNP.updateRIP(TestDepVar,obj.MODM);
                    if obj.RIPNorm, TestRIP = SNP.normalizeRIP(TestRIP,obj.MODM,obj.MOD.Var); end
                    index = find(~isnan(TestGG));
                    TestGG = TestGG(index);TestRIP = TestRIP(index)';
                 end
                % homozygotes only under an additive model
                 if obj.TestInd==4
                    index = find(abs(TrGG));
                    TrGG = TrGG(index);
                    TrRIP = TrRIP(index);
                    if sum(TrGG==-1)<=sum(TrGG==1)
                       PosClass = -1;
                    else
                       PosClass = 1;
                    end
                    index = find(abs(TestGG));
                    TestGG = TestGG(index);
                    TestRIP = TestRIP(index);
                 end
                 [out.EER,out.G,out.AUC,out.R] = oppGGBiometric(obj,TrGG,TrRIP,TestGG,TestRIP,PosClass);
                 if t==0; return; end
                 N = length(TestGG);
                 EER = nan*ones(1,t);
                 for i=1:1:t
                     ind = randsample(N,N,true);
                     EER(i) = oppGGBiometric(obj,TrGG,TrRIP,TestGG(ind),TestRIP(ind),PosClass);
                     if mod(i,100)==0
                        index = find(~isnan(EER)); 
                        tmppEER = (sum(EER(index)>=0.5)+1)/(length(index)+1); 
                        acc = 10/length(index);
                        if tmppEER>acc, break;end
                     end
                 end
                 index = find(~isnan(EER));
                 EER = EER(index);
                 EER = sort(EER);
                 out.AvgEER = nanmean(EER);
                 out.UpperEER = EER(round(0.975*length(index)));
                 out.LowerEER = EER(round(0.025*length(index)+1));
                 out.pEER = (sum(EER>=0.5)+1)/(length(index)+1);
       end
       function [LocD,LocP] = locationTest(obj,t,DepVar)
           GG = GT2GG(obj,obj.GT);
           if ~(obj.TestInd==4)
              [LocD,LocP] = USNP.DStatTest(DepVar(GG==-1,:),DepVar(GG==1,:),t); 
           else
              [LocD,LocP] = USNP.FStatTest(DepVar(GG==-1,:),DepVar(GG==0,:),DepVar(GG==1,:),t);
           end
       end
   end
   methods % MODEL BUILDING
       function buildFromAngleTest(obj,DepVar,COV,TestA)
          % INITIALIZING 
           GG = GT2GG(obj,obj.GT);
           ind = find(~isnan(GG));
           out.UsedInd = ind;
           GG = GG(ind);
           out.GG = GG;out.GT = obj.GT(ind);out.GGGTBal = USNP.getGGGTBalance(out.GG,out.GT);
           DepVar = DepVar(ind,:);
          % VAR INFORMATION 
           var.El = unique(GG);var.nrEl = length(var.El);
           var.Booting = size(size(COV,2)+1,2);var.nr2Boot = 1;
           var.Mel = min(var.El);var.Pel = max(var.El);var.Range = var.Pel-var.Mel;
           var.Mindex = find(GG==var.Mel);var.MDepVar = mean(DepVar(var.Mindex,:)); %#ok<*FNDSB>
           var.Pindex = find(GG==var.Pel);var.PDepVar = mean(DepVar(var.Pindex,:));
           var.nrWEl = zeros(1,var.nrEl);
           for i=1:1:var.nrEl
               var.nrWEl(i) = length(find(GG==var.El(i)));
           end
           [~,tmp] = min(var.nrWEl);var.PosClass = var.El(tmp);% needed for ROC analysis (minority group definition)
           out.Var = var;
          % M CONSTRUCTION
           fullM = [squeeze(TestA.M(:,1,:)) squeeze(TestA.M(:,2,:))];
           out.M = nanmedian(fullM,2)';
           out.H = USNP.wilcoxonM(fullM');
           out.HM = out.M.*out.H;
          % RIP COMPUTATIONS
           out.MRIP = USNP.updateRIP(DepVar,out.M)';
           out.NormMRIP = USNP.normalizeRIP(out.MRIP,out.M,out.Var);
           out.HMRIP = USNP.updateRIP(DepVar,out.HM)';
           out.NormHMRIP = USNP.normalizeRIP(out.HMRIP,out.HM,out.Var);
          % RIP STAT
           %[F,FP,~,Gmean,~,~] = USNP.performRIPStat(out.GG,out.MRIP,out.Var,0);
          % FINAL OUTPUT
           %[out.GGInfo,out.Freq] = USNP.getGTInfo(GG,RIP);
           [out.Balance,out.MinBalance] = USNP.getGGBalance(GG);
           obj.MOD = out;
       end
       function build(obj,DepVar,COV)
                % Run SNPv4BRIM
                 obj.MOD = getMOD(obj,DepVar,COV);
                 obj.MOD = rmfield(obj.MOD,{'DepVar' 'IndVar','Grp'});
               % Setting up reference scan 
                 %RIPref = SNPv4.updateRIP(BaseCont.ST4Ref.Coeff,obj.MOD.M);
                 %if obj.RIPNormalize, RIPref = SNPv4.normalizeRIP(RIPref,obj.MOD.M,obj.MOD.Var);end
                 %obj.MOD.RIPRef = RIPref;
               % Retrieving GT Info and statistics
                 [obj.MOD.GGInfo,obj.MOD.Freq] = USNP.getGTInfo(obj.MOD.GG,obj.MOD.RIP);  
       end
       function out = FOLDEDBRIM(obj,DepVar,Cov)
           % To bootstrap or not to bootstrap that is the question
              Bootstrap = true;if obj.MaxIterations == 0, Bootstrap = false;end
           % Defining Group Memberhip
              GG = getGG(obj,obj.GT); 
              ind = find(~isnan(GG));GG = GG(ind);
              switch obj.GroupDef
                  case 'GG' % Group membership is defined by the test performed
                      GM = GG;
                  case 'GT' % Group membership is defined by the original genotypes
                      GM = obj.GT(ind);
              end
              out.UsedInd = ind;% store the samples used
              out.GG = GG;out.GM = GM;
           % Defining Independent and Dependent Variables
              DepVar = DepVar(ind,:);
              if isempty(Cov);
                 IndVar = GG;
              else
                 IndVar = [Cov(ind,:) GG];
              end
              out.DepVar = DepVar;out.IndVar = IndVar;
            % Defining test variable info
              var.El = unique(GG);var.nrEl = length(var.El);
              var.Booting = size(IndVar,2);var.nr2Boot = 1;
              var.Mel = min(var.El);var.Pel = max(var.El);var.Range = var.Pel-var.Mel;
              var.Mindex = find(GG==var.Mel);var.MDepVar = mean(DepVar(var.Mindex,:)); %#ok<*FNDSB>
              var.Pindex = find(GG==var.Pel);var.PDepVar = mean(DepVar(var.Pindex,:));
              var.nrWEl = zeros(1,var.nrEl);
              for i=1:1:var.nrEl
                  var.nrWEl(i) = length(find(GG==var.El(i)));
              end
              [~,tmp] = min(var.nrWEl);var.PosClass = var.El(tmp);% needed for ROC analysis (minority group definition)
              out.Var = var;
            % Defining group info
              grp.El = unique(GM);grp.nrEl = length(grp.El);
              out.Grp = grp;
              [n,dim] = size(DepVar);
            % Outer Partitioning of Data
              Fouter = DOB_SCV(obj.OuterFold,DepVar,GM);
              FoldResults = cell(1,obj.OuterFold);
              parfor fo=1:obj.OuterFold
                 % Extracting training and testing 
                  FoTestInd = find(Fouter==fo);FoTrInd = setdiff((1:n),FoTestInd);nrFoTr = length(FoTrInd);
                  FoldResults{fo}.TestInd = FoTestInd;
                  FoldResults{fo}.TrInd = FoTrInd;
                  FoTrIndVar = IndVar(FoTrInd,:);FoTrDepVar = DepVar(FoTrInd,:);FoTrGM = GM(FoTrInd);
                  FoldResults{fo}.GG = IndVar(FoTestInd,end);% Only keeping origninal GG
                  FoldResults{fo}.TrGG = IndVar(FoTrInd,end);
                  FoTestDepVar = DepVar(FoTestInd,:);
                  %FoldResults{fo}.DepVar = FoTestDepVar;
                  BootProgress = zeros(1,obj.MaxIterations);
                 % Initialize Bootstrapping
                  %FoldResults{fo}.UpSampled = zeros(1,grp.nrEl);
                  ContBoot = true;counter = 0;
                  while ContBoot && Bootstrap>0
                    % Keeping track of iterations  
                     counter = counter + 1; 
                   % seperate data into inner folds
                     Finner = DOB_SCV(obj.InnerFold,FoTrDepVar,FoTrGM);
                     TMPIndVar = FoTrIndVar(:,end);% Allocate Memory, only last collumn for SNP
                     for fi=1:obj.InnerFold % Fi...
                         FiTestInd = find(Finner==fi);
                         FiTrInd = setdiff(1:nrFoTr,FiTestInd);
                         [M,H] = robustPLSRegression(obj,FoTrIndVar(FiTrInd,:),FoTrDepVar(FiTrInd,:),FoTrGM(FiTrInd));% Compute BootSampled Regression
                         if obj.Htest, M = H.*M;end% Only keep significant regression coefficients
                         rip = USNP.updateRIP(FoTrDepVar(FiTestInd,:),M)';% Get RIP scores
                         if obj.RIPNormalize, rip = USNP.normalizeRIP(rip,M,var);end% rescale RIP scores
                         TMPIndVar(FiTestInd) = rip;% Store updated RIP scores, always in the last collumn
                     end
                   % monitoring progress of Booting
                     tmpC = corrcoef(FoTrIndVar(:,end),TMPIndVar);
                     BootProgress(counter) = tmpC(1,2);
                     if BootProgress(counter)>=obj.StopCorr, ContBoot = false; end
                     if counter >= obj.MaxIterations, ContBoot = false; end
                     FoTrIndVar(:,end) = TMPIndVar;% Update Training IndVar for next round
                  end
                  %if counter>0, FoldResults{fo}.BootProgress = BootProgress(1:counter);end
                  %FoldResults{fo}.BootIter = counter;
                 % Extracting final M from outer training data
                  [M,H] = robustPLSRegression(obj,FoTrIndVar,FoTrDepVar,FoTrGM);
                  if obj.Htest, M = H.*M;end
                  FoldResults{fo}.M = M;FoldResults{fo}.H = H;
                 % Evaluate outer test data
                  rip = updateRIP(FoTestDepVar,M)';% Get RIP Scores Test partition 
                  if obj.RIPNormalize, rip = USNP.normalizeRIP(rip,M,var);end% rescale RIP scores
                  FoldResults{fo}.RIP = rip;
                  [FoldResults{fo}.F,FoldResults{fo}.FP,FoldResults{fo}.FPar,FoldResults{fo}.G,FoldResults{fo}.GP,~] = USNP.performRIPStat(FoldResults{fo}.GG,FoldResults{fo}.RIP,var,obj.RIPStatPerm);
                  rip = updateRIP(FoTrDepVar,M)';% Get RIP Scores training partition
                  if obj.RIPNormalize, rip = USNP.normalizeRIP(rip,M,var);end% rescale RIP scores
                  FoldResults{fo}.TrRIP = rip;
              end
             % Gathering Fold Results
               FoldF = zeros(1,obj.OuterFold);
               FoldFP = zeros(1,obj.OuterFold);
               FoldFPar = zeros(1,obj.OuterFold);
               FoldG = zeros(1,obj.OuterFold);
               FoldGP = zeros(1,obj.OuterFold);
               Folds = cell(1,obj.OuterFold);
               M = zeros(obj.OuterFold,dim);
               for fo=1:obj.OuterFold
                   FoldF(fo) = FoldResults{fo}.F;FoldFP(fo) = FoldResults{fo}.FP;FoldFPar(fo) = FoldResults{fo}.FPar;
                   FoldG(fo) = FoldResults{fo}.G;FoldGP(fo) = FoldResults{fo}.GP;M(fo,:) = FoldResults{fo}.M;
                   Folds{fo}.M = FoldResults{fo}.M;Folds{fo}.GG = FoldResults{fo}.GG;Folds{fo}.RIP = FoldResults{fo}.RIP;
                   Folds{fo}.TestInd = FoldResults{fo}.TestInd;Folds{fo}.TrGG = FoldResults{fo}.TrGG;Folds{fo}.TrRIP = FoldResults{fo}.TrRIP;
               end
               out.F = nanmean(FoldF);
               out.G = nanmean(FoldG);
               out.FPar = nan;
               index = find(~isnan(FoldFPar)); %#ok<*EFIND>
               if ~isempty(index), out.FPar = USNP.pfast(FoldFPar); end
               out.FP = nan;
               out.GP = nan;
               if obj.RIPStatPerm>0
                  FoldFP(FoldFP==0) = 0.1/obj.RIPStatPerm;
                  index = find(~isnan(FoldFP));
                  if ~isempty(index), out.FP = USNP.pfast(FoldFP(index)); end
                  FoldGP(FoldGP==0) = 0.1/obj.RIPStatPerm;
                  index = find(~isnan(FoldGP));
                  if ~isempty(index), out.GP = USNP.pfast(FoldGP(index)); end
               end
               out.FoldF = FoldF;
               out.FoldG = FoldG;
               out.FoldFPar = FoldFPar;
               out.FoldFP = FoldFP;
               out.FoldGP = FoldGP;
               out.Folds = Folds;
              % Extracting final M 
               out.M = median(M);
               if ~obj.Htest, return; end
               if obj.OuterFold>=8, 
                  out.H = USNP.wilcoxonM(M);
               else
                  out.H = ones(1,dim); 
               end
               out.M = out.H.*out.M;
               RIP = USNP.updateRIP(DepVar,out.M);
               if obj.RIPNormalize, RIP = USNP.normalizeRIP(RIP,out.M,var); end
               out.RIP = RIP;
               [out.RIPF,out.RIPFP,out.RIPFPar,out.RIPG,out.RIPGP,~] = USNP.performRIPStat(out.GG,RIP,var,0);
              % Obtaining regressions to contruct morphs
               %IndVar = [Cov(out.UsedInd,:) RIP'];
               %[out.MMorphs,out.MIMorphs,~] = getRegression(IndVar,DepVar,(1:size(IndVar,2)));   
       end
       function out = BRIM(obj,DepVar,Cov)
           % To bootstrap or not to bootstrap that is the question
              Bootstrap = true;if obj.MaxIterations == 0, Bootstrap = false;end
           % Defining Group Memberhip
              GG = getGG(obj,obj.GT); 
              ind = find(~isnan(GG));GG = GG(ind);
              switch obj.GroupDef
                  case 'GG' % Group membership is defined by the test performed
                      GM = GG;
                  case 'GT' % Group membership is defined by the original genotypes
                      GM = obj.GT(ind);
              end
              out.UsedInd = ind;% store the samples used
              out.GG = GG;out.GM = GM;
           % Defining Independent and Dependent Variables
              DepVar = DepVar(ind,:);
              if isempty(Cov);
                 IndVar = GG;
              else
                 IndVar = [Cov(ind,:) GG];
              end
              out.DepVar = DepVar;out.IndVar = IndVar;
            % Defining test variable info
              var.El = unique(GG);var.nrEl = length(var.El);
              var.Booting = size(IndVar,2);var.nr2Boot = 1;
              var.Mel = min(var.El);var.Pel = max(var.El);var.Range = var.Pel-var.Mel;
              var.Mindex = find(GG==var.Mel);var.MDepVar = mean(DepVar(var.Mindex,:)); %#ok<*FNDSB>
              var.Pindex = find(GG==var.Pel);var.PDepVar = mean(DepVar(var.Pindex,:));
              var.nrWEl = zeros(1,var.nrEl);
              for i=1:1:var.nrEl
                  var.nrWEl(i) = length(find(GG==var.El(i)));
              end
              [~,tmp] = min(var.nrWEl);var.PosClass = var.El(tmp);% needed for ROC analysis (minority group definition)
              out.Var = var;
            % Defining group info
              grp.El = unique(GM);grp.nrEl = length(grp.El);
              out.Grp = grp;
              [n,dim] = size(DepVar); %#ok<ASGLU>
            % Initialize Bootstrapping           
              BootProgress = zeros(1,obj.MaxIterations);
              ContBoot = true;counter = 0;
              while ContBoot && Bootstrap>0
                    % Keeping track of iterations  
                      counter = counter + 1; 
                    % seperate data into inner folds
                     Finner = DOB_SCV(obj.InnerFold,DepVar,GM);
                     TMPIndVar = IndVar(:,end);% Allocate Memory, only last collumn for SNP
                     for fi=1:obj.InnerFold % Fi...
                         FiTestInd = find(Finner==fi);
                         FiTrInd = setdiff(1:n,FiTestInd);
                         [M,H] = robustPLSRegression(obj,IndVar(FiTrInd,:),DepVar(FiTrInd,:),GM(FiTrInd));% Compute BootSampled Regression
                         if obj.Htest, M = H.*M;end% Only keep significant regression coefficients
                         rip = USNP.updateRIP(DepVar(FiTestInd,:),M)';% Get RIP scores
                         if obj.RIPNormalize, rip = USNP.normalizeRIP(rip,M,var);end% rescale RIP scores
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
              [M,H] = robustPLSRegression(obj,IndVar,DepVar,GM);
              if obj.Htest, M = H.*M;end
              out.M = M;out.H = H;
              RIP = USNP.updateRIP(DepVar,out.M);
              if obj.RIPNormalize,RIP = USNP.normalizeRIP(RIP,out.M,var);end
              out.RIP = RIP;
              [out.RIPF,out.RIPFP,out.RIPFPar,out.RIPG,out.RIPGP,~] = USNP.performRIPStat(out.GG,RIP,var,0);  
       end
       function out = PLSR(obj,DepVar,Cov)
           % Defining Group Memberhip
              GG = GT2GG(obj,obj.GT); 
              ind = find(~isnan(GG));GG = GG(ind);
              switch obj.GroupDef
                  case 'GG' % Group membership is defined by the test performed
                      GM = GG;
                  case 'GT' % Group membership is defined by the original genotypes
                      GM = obj.GT(ind);
              end
              out.UsedInd = ind;% store the samples used
              out.GG = GG;out.GM = GM;
           % Defining Independent and Dependent Variables
              DepVar = DepVar(ind,:);
              if isempty(Cov);
                 IndVar = GG;
              else
                 IndVar = [Cov(ind,:) GG];
              end
              out.DepVar = DepVar;out.IndVar = IndVar;
            % Defining test variable info
              var.El = unique(GG);var.nrEl = length(var.El);
              var.Booting = size(IndVar,2);var.nr2Boot = 1;
              var.Mel = min(var.El);var.Pel = max(var.El);var.Range = var.Pel-var.Mel;
              var.Mindex = find(GG==var.Mel);var.MDepVar = mean(DepVar(var.Mindex,:)); %#ok<*FNDSB>
              var.Pindex = find(GG==var.Pel);var.PDepVar = mean(DepVar(var.Pindex,:));
              var.nrWEl = zeros(1,var.nrEl);
              for i=1:1:var.nrEl
                  var.nrWEl(i) = length(find(GG==var.El(i)));
              end
              [~,tmp] = min(var.nrWEl);var.PosClass = var.El(tmp);% needed for ROC analysis (minority group definition)
              out.Var = var;
            % Defining group info
              grp.El = unique(GM);grp.nrEl = length(grp.El);
              out.Grp = grp;
              [n,dim] = size(DepVar); %#ok<ASGLU>
            % Extracting final M from outer training data
              [M,H] = robustPLSRegression(obj,IndVar,DepVar,GM,obj.SamplePercentage);
              if obj.Htest, M = H.*M;end
              out.M = M;out.H = H;
              RIP = USNP.updateRIP(DepVar,out.M);
              if obj.RIPNorm, RIP = USNP.normalizeRIP(RIP,out.M,var);end
              out.RIP = RIP;
              [out.RIPF,out.RIPFP,out.RIPFPar,out.RIPG,out.RIPGP,~] = USNP.performRIPStat(out.GG,RIP,var,0);  
       end
       function [M,H,Mfor] = robustPLSRegression(obj,indvar,depvar,gg,perc)
             if nargin<5, perc = 1; end
           % initialize
             sampling = obj.RegSampling;runs = obj.RegRuns;balance = obj.RegBalance;% only read once
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
                                    distances = squareform(pdist(depvar(tmpindm,:)));
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
                            [DepVarNew,IndVarNew] = USNP.regUpSample(depvar(balInfo{i}.Ind,:),indvar(balInfo{i}.Ind,:),balInfo{i});
                            SampleDepVar = [SampleDepVar; DepVarNew]; %#ok<*AGROW>
                            SampleIndVar = [SampleIndVar; IndVarNew];
                        end
                    else % Balancing by Decreasing Majority Groups
                        SampleIndVar = [];SampleDepVar = [];
                        for i=1:1:nrgg
                            [DepVarNew,IndVarNew] = USNP.regDownSample(depvar(balInfo{i}.Ind,:),indvar(balInfo{i}.Ind,:),balInfo{i}.nrKeep);
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
                 Mfor(s,:) = USNP.getRegression(SampleIndVar(SampleInd,:),SampleDepVar(SampleInd,:));
             end
           % finalizing output
             M = median(Mfor,1); % extract the median of all regression coefficients
             if ~obj.Htest, H=ones(1,nrD);return; end
             if runs==1, H=ones(1,nrD); return; end
             if nargout<2, return; end
             H = USNP.wilcoxonM(Mfor);% wilcoxon test for median
       end       
       function out = getMOD(obj,DepVar,Cov)
          switch obj.BuildMethod
              case 'PLSR'
                  out = PLSR(obj,DepVar,Cov);
              case 'BRIM'
                  out = BRIM(obj,DepVar,Cov);
              case 'FOLDEDBRIM'
                  out = FOLDEDBRIM(obj,DepVar,Cov);
          end
       end
   end
   methods % BIOMETRICS
       function [EER,G,AUC,R] = oppGGBiometric(obj,trgg,trrip,gg,rip,posclass)
                n = length(rip);
                C = unique(trgg);nrC = length(C);allC = 1:nrC;
                AVG = nan*zeros(1,nrC);STD = nan*zeros(1,nrC);
                for c=1:1:nrC
                    trindex = find(trgg==C(c));
                    AVG(c) = mean(trrip(trindex));
                    STD(c) = std(trrip(trindex));
                end
                PDF = normpdf(repmat(rip,1,nrC),repmat(AVG,n,1),repmat(STD,n,1));
                % true matches only defined by minority class
                cp = find(C==posclass);
                indexp = find(gg==posclass);
                indexn = setdiff(1:n,indexp);
                remCp = setdiff(allC,cp);
                switch obj.Match
                    case 'QD'
                        ratio = (PDF(:,cp))./PDF(:,remCp);
                        switch obj.Classify
                            case 'SOFT'
                                ratio = -log(ratio);
                            case 'HARD'
                                ratio = ratio<1;
                        end
                    case 'INLIER'
                        ratio = (PDF(:,cp))./sum(PDF,2);
                        switch obj.Classify
                            case 'SOFT'
                                ratio = 1-ratio;
                            case 'HARD'
                                ratio = ratio<0.5;
                        end
                end
                [EER,G,AUC] = USNP.getEER(ratio(indexp),ratio(indexn));
                R = USNP.getRANK(ratio(indexp),ratio(indexn));
       end
   end
   methods % MODEL APPLICATION
   end
   methods % INTERFACING     
       function obj = reduceSamples(obj,index)
           if nargout==1, obj = clone(obj); end
           obj.GT = obj.GT(index);
           obj.MOD = [];
       end
       function saveToDisk(obj,suffix)
          save([obj.RS '_' suffix],'obj');           
       end
       function out = GT2GG(obj,in)
           out = in;% Initialize
           mmind = find(in==obj.mmGT);
           Mmind = find(in==obj.MmGT);
           MMind = find(in==obj.MMGT);
           switch obj.TestInd
               case 0 % NO EFFECT
                   out = nan*out;
               case 1 % test 1: AAAB BB (RECESSIVE)
                   out(union(Mmind,MMind)) = -1;
                   out(mmind) = 1;
               case 2 % test 2: AA ABBB (DOMINANT)
                   out(MMind) = -1;
                   out(union(mmind,Mmind)) = 1;
               case 3 % test 3: AABB AB (OVER DOMINANT)
                   out(union(mmind,MMind)) = -1;
                   out(Mmind) = 1;
               case 4 % test 4 AA AB BB (ADDITIVE)
                   out(mmind) = 1;
                   out(Mmind) = 0;
                   out(MMind) = -1;
               case 5 % test 5: AA BB
                   out(mmind) = 1;
                   out(Mmind) = nan;
                   out(MMind) = -1;
               case 6 % test 6: AA AB
                   out(mmind) = 1;
                   out(Mmind) = -1;
                   out(MMind) = nan;
               case 7 % test 7: AB BB
                   out(mmind) = nan;
                   out(Mmind) = 1;
                   out(MMind) = -1;
           end       
       end
       function [Distr,Freq,RIP] = GT2RIP(obj,in)
           Distr = {};Freq = nan;RIP = nan;
           if isempty(obj.MOD), return; end
           gg = GT2GG(obj,in);
           if isnan(gg),return; end
           if (obj.TestInd==4)&&(gg==0),return;end% eliminate heterozygotes for additive models for now
           [Distr,Freq,RIP] = USNP.GG2RIP(obj.MOD.GG,getMODRIP(obj,obj.MOD),gg);
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
   end
   methods (Static = true)
       function [R2,pR2,M] = testPartialPLSR(X,Y,C,t)
           index = intersect(BASE.notNAN(C),BASE.notNAN(X));
           E = BASE.getResiduals(C(index,:),Y(index,:));
           X = BASE.getResiduals(C(index,:),X(index,:));
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
       function [D,P] = DStatTest(X1,X2,t)
         % initializing  
           nX1 = size(X1,1);nX2 = size(X2,1);
         % Getting First Order moments & Residus
           AvgX1 = mean(X1);AvgX2 = mean(X2);
         % Getting Distance Based Test Statistic
           D =  sqrt((AvgX1-AvgX2)*(AvgX1-AvgX2)');
         % Permutation test  
           StatCount = false(1,t);
           nT = nX1+nX2;
           X = [X1; X2];
           parfor i=1:t
                  ind = randperm(nT);
                  X1perm = X(ind(1:nX1),:);X2perm = X(ind(nX1+1:end),:); %#ok<*PFBNS>
                  AvgX1perm = mean(X1perm);AvgX2perm = mean(X2perm);
                  Dperm = sqrt((AvgX1perm-AvgX2perm)*(AvgX1perm-AvgX2perm)');
                  StatCount(i) = Dperm>=D;
           end
           P = sum(StatCount)/t;           
       end
       function [F,P] = FStatTest(X1,X2,X3,t)
         % initializing  
           nX1 = size(X1,1);nX2 = size(X2,1);nX3 = size(X3,1);nX = nX1+nX2+nX3;
           X = [X1;X2;X3];
         % group averages
           AvgX1 = mean(X1);AvgX2 = mean(X2);AvgX3 = mean(X3);AvgX = mean(X);   
         % SSB
           SSB = sum(sum(([AvgX1;AvgX2;AvgX3] - repmat(AvgX,3,1)).^2,2).*[nX1;nX2;nX3])/2;
         % SSW
           SSW = sum([sum((X1-repmat(AvgX1,nX1,1)).^2,2);sum((X2-repmat(AvgX2,nX2,1)).^2,2);sum((X3-repmat(AvgX3,nX3,1)).^2,2)])/(nX-2);
         % F
           F = SSB/SSW;
         % Permutation test  
           StatCount = false(1,t);
           parfor i=1:t
                 % initialize permutation
                  ind = randperm(nX);
                  X1perm = X(ind(1:nX1),:);X2perm = X(ind(nX1+1:nX1+nX2),:);X3perm = X(ind(nX1+nX2+1:end),:); %#ok<*PFBNS>
                 % group averages 
                  AvgX1perm = mean(X1perm);AvgX2perm = mean(X2perm);AvgX3perm = mean(X3perm);AvgXperm = mean(X(ind,:));
                 % SSB
                  SSBperm = sum(sum(([AvgX1perm;AvgX2perm;AvgX3perm] - repmat(AvgXperm,3,1)).^2,2).*[nX1;nX2;nX3])/2;
                 % SSW 
                  SSWperm = sum([sum((X1perm-repmat(AvgX1perm,nX1,1)).^2,2);sum((X2perm-repmat(AvgX2perm,nX2,1)).^2,2);sum((X3perm-repmat(AvgX3perm,nX3,1)).^2,2)])/(nX-2);
                 % F
                  StatCount(i) = (SSBperm/SSWperm)>=F;
           end
           P = sum(StatCount)/t;
       end
       function out = getBalance(nr1,nr2)
                out = min([nr1 nr2])./max([nr1 nr2]);
       end
       function [GAdd,CAdd] = foldUpSample(GO,CO,nr)        
           nO = size(GO,1);
           GAdd = zeros(nr,size(GO,2));
           CAdd = zeros(nr,size(CO,2));
           s = randsample(nO,nr,true);
           alpha = rand(1,nr);
           for f=1:1:nr
               nb = randsample(setdiff(1:nO,s(f)),1);
               GAdd(f,:) = GO(s(f),:) + (GO(nb,:)-GO(s(f),:))*alpha(f);% create new Face
               CAdd(f,:) = CO(s(f),:) + (CO(nb,:)-CO(s(f),:))*alpha(f);% create new Covariates
           end
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
       function M = getRegression(A,B)
           [A,B] = eliminateNAN(A,B);% typically present in covariates
           [~,~,~,~,M] = plsregress(A,B,size(A,2));
           M = M(end,:);
       end
       function H = wilcoxonM(M)
            nrD = size(M,2);
            %H = zeros(1,nrD);
            P = zeros(1,nrD);
            for j=1:1:nrD
               %[P(j),H(j)] = signrank(M(:,j));% wilcoxon test for median  
               P(j) = signrank(M(:,j));% wilcoxon test for median  
            end
            H = P<=0.00001;
       end
       function [out] = normalizeRIP(rip,M,var)
             Mrip = dot(var.MDepVar',M'/norm(M'));
             Prip = dot(var.PDepVar',M'/norm(M'));
             out = ((rip-Mrip)/(Prip-Mrip))*var.Range+var.Mel;
       end
       function [out] = updateRIP(in,M)
            out = dot(in',repmat(M'/norm(M'),1,size(in,1)));
       end
       function [F,FP,FPar,Gmean,GmeanP,GT] = performRIPStat(g,rip,var,t)
                FP = nan;Gmean = nan;GmeanP = nan;GT = nan;
                if var.nrEl<3, doroc = true; else doroc = false;end
                index = find(~isnan(g));g = g(index);
                XG = cell(size(g));
                for l=1:var.nrEl % building grouping variable
                    XG((g==var.El(l))) = {num2str(l)};
                end
                rip = rip(index);n = length(index);
                [F,FPar] = USNP.ripAnova(XG,rip);
                if isnan(F), return; end % BAD TEST GROUP
                if doroc, [Gmean,GT] = USNP.ripROC(g,rip,var.PosClass); end
                if t == 0,return;end
                FCount = zeros(1,t);
                if doroc, GCount = zeros(1,t); end
                parfor i=1:t
                    ind = randperm(n);
                    Ffor = USNP.ripAnova(XG,rip(ind));
                    FCount(i) = Ffor>=F;
                    if doroc, Gmeanfor = USNP.ripROC(g,rip(ind),var.PosClass); GCount(i) = Gmeanfor>=Gmean; end
                end
                FP = sum(FCount)/t;
                if doroc, GmeanP = sum(GCount)/t; end
       end
       function [F,P] = ripAnova(XG,rip)
                [P,T] = anova1(rip,XG,'off');
                F = T{2,5};  
       end
       function [G,T] = ripROC(g,rip,posclass)
                g = g*posclass;rip = rip*posclass;% Minority group must be positive class
                [rip,order] = sort(rip);g = g(order);
                true_neg = g == -1;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
                true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
                y = tpf(1:end-1)+dy./2;x = fpf(1:end-1)+dx./2;
                [~,ind] = min(abs(x-(1-y)));
                FN = tpf(ind+1)*na;TN = fpf(ind+1)*nn;
                TP = na-FN;FP = nn-TN;
                G = sqrt((TP/(TP+FN))*(TN/(TN+FP))); % G-mean measure for imbalanced datasets
                if nargout>1, T=rip(ind+1)*posclass;end % classification is done -1 if rip<=T and 1 if rip > T     
       end
       function [G,T] = ripROCFull(g,rip,posclass,testg,testrip)
                g = g*posclass;rip = rip*posclass;% Minority group must be positive class
                testg = testg*posclass;testrip = testrip*posclass;
                [rip,order] = sort(rip);g = g(order);
                true_neg = g == -1;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
                true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
                y = tpf(1:end-1)+dy./2;x = fpf(1:end-1)+dx./2;
                [~,ind] = min(abs(x-(1-y)));
                FN = tpf(ind+1)*na;TN = fpf(ind+1)*nn;
                TP = na-FN;FP = nn-TN;
                G = sqrt((TP/(TP+FN))*(TN/(TN+FP))); % G-mean measure for imbalanced datasets
                nr = length(x)-1;
                GR = zeros(1,nr);
                AC = zeros(1,nr);
                PR = zeros(1,nr);
                RE = zeros(1,nr);
                for i=1:1:nr
                    FNR = tpf(i+1)*na;TNR = fpf(i+1)*nn;
                    TPR = na-FN;FPR = nn-TNR;
                    GR(i) = sqrt((TPR/(TPR+FNR))*(TNR/(TNR+FPR)));
                    AC(i) = (TPR+TNR)/(TPR+FNR+FPR+TNR);
                    PR(i) = TPR/(TPR+FPR);
                    RE(i) = TPR/(TPR+FNR);
                end
                figure;hold on;
                plot(1:nr,GR,'b-');plot(1:nr,PR,'g-');plot(1:nr,RE,'r-');plot(1:nr,AC,'k-');
                [MGR,Mind] = max(GR);
                if nargout>1, T=rip(ind+1)*posclass;end % classification is done -1 if rip<=T and 1 if rip > T     
       end
       function [Soft,Hard] = evalQD(trgg,trrip,gg,rip,posclass,t)
                if nargin < 6, t = 0; end
                n = length(rip);
                C = unique(trgg);nrC = length(C);allC = 1:nrC;
                AVG = nan*zeros(1,nrC);STD = nan*zeros(1,nrC);
                for c=1:1:nrC
                    trindex = find(trgg==C(c));
                    AVG(c) = mean(trrip(trindex));
                    STD(c) = std(trrip(trindex));
                end
                PDF = normpdf(repmat(rip,1,nrC),repmat(AVG,n,1),repmat(STD,n,1));
                % true matches only defined by minority class
                % SOFT ASSESSMENT
                cp = find(C==posclass);
                indexp = find(gg==posclass);
                indexn = setdiff(1:n,indexp);
                remCp = setdiff(allC,cp);
                ratio = (length(remCp)*PDF(:,cp))./sum(PDF(:,remCp),2);
                tmatches = -log(ratio(indexp));fmatches = -log(ratio(indexn));
                [Soft.EER,Soft.G,Soft.AUC] = USNP.getEER(tmatches,fmatches);
                Soft.R = USNP.getRANK(tmatches,fmatches);
                % HARD ASSESSMENT
                [Hard.PREC,Hard.REC,Hard.G,Hard.TNF] = USNP.getHardClass(ratio,indexp,indexn,1);
                if t==0, return; end
                if isnan(Soft.EER), return; end
                EERcount = zeros(1,t);
                AUCcount = zeros(1,t);
                Rcount = zeros(1,t);
                Gcount = zeros(1,t);
                parfor i=1:t
                   ind = randperm(n);
                   forPDF = normpdf(repmat(rip(ind),1,nrC),repmat(AVG,n,1),repmat(STD,n,1));
                   forratio = (length(remCp)*forPDF(:,cp))./sum(forPDF(:,remCp),2);
                   fortmatches = -log(forratio(indexp));forfmatches = -log(forratio(indexn)); 
                   [forEER,~,forAUC] = USNP.getEER(fortmatches,forfmatches);
                   tmp = USNP.getRANK(fortmatches,forfmatches);
                   forR = tmp(4);
                   EERcount(i) = forEER<=Soft.EER;
                   AUCcount(i) = forAUC>=Soft.AUC;
                   Rcount(i) = forR<=Soft.R(4);
                   [~,~,forG] = USNP.getHardClass(forratio,indexp,indexn,1);
                   Gcount(i) = forG>=Hard.G;
                end
                Soft.pAUC = sum(AUCcount)/t;
                Soft.pEER = sum(EERcount)/t;
                Hard.pG = sum(Gcount)/t;
                Soft.pR = sum(Rcount)/t;
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
       function [out,freq] = getGTInfo(gg,rip)
              % Obtaining group memberships and additional needed info
                 elGG = unique(gg);
                 nrGG = length(elGG);
                 out = cell(1,nrGG);
                 freq = [];
                 for i=1:1:nrGG
                     out{i}.GG = elGG(i);
                     out{i}.Ind = find(gg==out{i}.GG);
                     out{i}.Freq = length(out{i}.Ind)/length(find(~isnan(gg)));
                     out{i}.Nr = length(out{i}.Ind);
                     if out{i}.Freq==0, out{i}.Freq = nan; continue;end
                     [out{i}.AVG,out{i}.STD,out{i}.MED,out{i}.MAD] = USNP.extractRIPStatistics(rip(out{i}.Ind));
                     freq = [freq out{i}.Freq];
                 end
       end
       function [Distr,Freq,RIP] = GG2RIP(trGG,trrip,GG)
           Tind = find(trGG==GG);% Inlier Distribution
           Find = find(trGG==-1*GG);% Outlier Distribution
           [Distr.Tmu,Distr.Tsigma,Distr.TO] = USNP.extractRIPStatistics(trrip(Tind)');
           RIP = Distr.Tmu;Distr.TpMax = normpdf(Distr.Tmu,Distr.Tmu,Distr.Tsigma);
           [Distr.Fmu,Distr.Fsigma,Distr.FO] = USNP.extractRIPStatistics(trrip(Find)');
           Distr.FpMax = normpdf(Distr.Fmu,Distr.Fmu,Distr.Fsigma);
           Freq = length(Tind)/(length(Tind)+length(Find));
       end
       function [AVG,STD,MED,MAD,O] = extractRobustRIPStatistics(rip,kappa)
                % Initialize
                  O.avg = mean(rip);O.std = std(rip);
                  MED = median(rip);MAD = mad(rip);
                  nrO = length(rip);W = ones(1,nrO);
                  for i=1:1:10
                     AVG = sum(W.*rip)/sum(W); 
                     STD = sqrt(sum(W.*((rip-AVG).^2))/sum(W));
                     IP = normpdf(rip,AVG,STD);
                     %L = (1/sqrt(((2*pi)^2)*STD))*exp(-0.5*kappa^2);
                     L = (1/(sqrt(2*pi)*STD))*exp(-0.5*kappa^2);
                     W = IP./(IP+L);
                  end
       end
       function [AVG,STD,MED,MAD] = extractRIPStatistics(rip)
                % Initialize
                  AVG = mean(rip);STD = std(rip);
                  MED = median(rip);MAD = mad(rip);
       end
       function [EER,G,AUC,x,y] = getEER(tmatches,fmatches)
                g = [ones(1,length(tmatches)), -1*ones(1,length(fmatches))];
                [~,order] = sort([tmatches;fmatches],'ascend');
                g = g(order);
                true_neg = g == -1;nn = sum(true_neg);fpf = scale(cumsum(true_neg),nn);dx = diff(fpf);
                true_pos = g == 1;na = sum(true_pos);tpf = scale(cumsum(true_pos),na);dy = diff(tpf);
                y = tpf(1:end-1)+dy./2;
                x = fpf(1:end-1)+dx./2;
                yn = 1-y;d = abs(x-yn);
                [~,ind] = min(d);
                EER = ((x(ind)+yn(ind))/2);
                AUC = sum(dx.*y);
                FN = tpf(ind+1)*na;TN = fpf(ind+1)*nn;
                TP = na-FN;FP = nn-TN;
                G = 1-sqrt((TP/(TP+FN))*(TN/(TN+FP))); 
       end
       function out = getRANK(tmatches,fmatches)
                R = sum(repmat(tmatches,1,length(fmatches))>repmat(fmatches',length(tmatches),1),2);
                R = R./(length(fmatches)+1);
                RM = mean(R);
                R1 = sum(R<=0.01)/length(R);
                R10 = sum(R<=0.10)/length(R);
                R20 = sum(R<=0.20)/length(R);
                out = [R1 R10 R20 RM];
       end
       function [PREC,REC,G,TNF] = getHardClass(ratio,indexp,indexn,Thres)
                Pc = length(indexp);Nc = length(indexn);
                TP = sum(ratio(indexp)>Thres);
                TN = sum(ratio(indexn)<Thres);TNF = TN/Nc;
                FN = Pc-TP;FP = Nc-TN;
                PREC = TP/(TP+FP);
                REC = TP/(TP+FN);
                G = sqrt((TP/(TP+FN))*(TN/(TN+FP)));
       end
       function [bal,minbal] = getGGBalance(gg)
          elgg = unique(gg);
          nrgg = length(elgg);
          nS = length(gg);
          bal = nan*zeros(1,nrgg);
          for i=1:1:nrgg
              bal(i) = length(find(gg==elgg(i)))/nS;
          end
          minbal = min(bal);
       end
       function out = getGGGTBalance(gg,gt)
           elgg = unique(gg);
           nrelgg = length(elgg);
           out = zeros(1,nrelgg);
           for i=1:1:nrelgg
              ind = find(gg==elgg(i));
              forgt = gt(ind);
              elforgt = unique(forgt);
              nrelforgt = length(elforgt);
              if nrelforgt==1, out(i) = 1; continue; end
              forout = zeros(1,nrelforgt);
              for j=1:1:nrelforgt
                 forout(j) = length(find(forgt==elforgt(j))); 
              end
              out(i) = min(forout)/max(forout);
           end
       end
   end
end


% function buildFromAngleTest(obj,DepVar,COV,TestA,Htest)
%           % INITIALIZING 
%            GG = GT2GG(obj,obj.GT);
%            ind = find(~isnan(GG));
%            out.UsedInd = ind;
%            GG = GG(ind);
%            out.GG = GG;out.GT = obj.GT(ind);out.GGGTBal = USNP.getGGGTBalance(out.GG,out.GT);
%            DepVar = DepVar(ind,:);
%            forM = obj.TestA.M;
%            if Weigthed
%                forW = obj.TestA.A';
%                forW(forW<=0) = 0;forW(isnan(forW)) = 0;
%            else
%                forW = ones(1,length(obj.TestA.A));
%            end
%            forW = forW./sum(forW);
%            var.El = unique(GG);var.nrEl = length(var.El);
%            var.Booting = size(size(COV,2)+1,2);var.nr2Boot = 1;
%            var.Mel = min(var.El);var.Pel = max(var.El);var.Range = var.Pel-var.Mel;
%            var.Mindex = find(GG==var.Mel);var.MDepVar = mean(DepVar(var.Mindex,:)); %#ok<*FNDSB>
%            var.Pindex = find(GG==var.Pel);var.PDepVar = mean(DepVar(var.Pindex,:));
%            var.nrWEl = zeros(1,var.nrEl);
%            for i=1:1:var.nrEl
%                var.nrWEl(i) = length(find(GG==var.El(i)));
%            end
%            [~,tmp] = min(var.nrWEl);var.PosClass = var.El(tmp);% needed for ROC analysis (minority group definition)
%            out.Var = var;
%            forM = squeeze(mean(forM,2));
%            %[~,maxind] = max(forW);% select the best
%            %out{i}.M = forM(:,maxind)';
%            %out.M = nanmedian(forM,2)';
%            ShapeDim = size(forM,1);
%            M = zeros(1,ShapeDim);
%            for i=1:1:ShapeDim
%                M(i) = weightedMedian(forM(i,:),forW);
%            end
%            out.M = M;
%            %out{i}.M = (mean(repmat(forW,obj.ShapeDim,1).*forM,2)./repmat(sum(forW),obj.ShapeDim,1))';
%            RIP = USNP.updateRIP(DepVar,out.M)';
%            if obj.RIPNorm, RIP = USNP.normalizeRIP(RIP,out.M,out.Var);end
%            out.RIP = RIP;
%            [out.GGInfo,out.Freq] = USNP.getGTInfo(GG,RIP);
%            [out.Balance,out.MinBalance] = USNP.getGGBalance(GG);
%            obj.MOD = out;
%        end
