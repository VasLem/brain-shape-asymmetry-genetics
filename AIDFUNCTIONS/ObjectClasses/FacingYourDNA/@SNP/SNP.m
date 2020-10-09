classdef SNP < superClass
   % General Properties
   properties
      RS;
      GT;
   end
   properties (Dependent = true)
      %AAavgScan;
      %BBavgScan;
      %ABavgScan;
      GTRIP;
      nrAA;
      nrAB;
      nrBB; 
   end
   properties (Hidden = true, Dependent = true)
      Valid;
      Balances;
   end
   properties (Hidden = true, Transient = true)
      ParCont = [];
   end
   properties (Hidden = true)
      TestLabels = {'AAAB BB' 'AA ABBB' 'AABB AB' 'AA AB BB' 'AA BB' 'AA AB' 'AB BB'};
      %GG = [];
   end
   properties (Hidden = true, Dependent = true)
      RegRuns; % Number of Regression runs (SNPBRIM)
      RegSampling; % Bootstrap/None (SNPBRIM)
      RegBalance;  % balance the data yes or no (SNPBRIM)
      RegBalanceTH; % desired balance factor (SNPBRIM)
      RegBalanceMethod; % UpSample/DownSample/ADASYN (SNPBRIM)
      OuterFold; % number of Outer Folds (SNPBRIM)
      InnerFold; % number of Inner Folds (SNPBRIM)
      UseRedShape; % Use reduced space or covariates (SNPBRIM)
      MinFoldSampleSize; % Min samples per Outer Fold and Inner Fold (SNPBRIM)
      MaxIterations; % Maximum number of BRIM iterations (SNPBRIM)
      Htest; % Perform Wilcoxon test om partial regression coefficients (SNPBRIM)
      RIPNormalize; % normalize RIP values based on group averages (SNPBRIM)
      StopCorr; % Stopping correlation between subsequent iterations (SNPBRIM)
      RIPStatPerm; % Number of permutations in testing significance of RIP values (SNPBRIM)
      GroupDef;% Definition of groups, GG = according to Test performed, GT = according to orginal genotypes (SNPBRIM)
      %ShapeSpace;
      %RedShapeSpace;
      %DepVar;
      %RedDepVar;
      %Cov;
      %CovNames;
      AAind;
      ABind;
      BBind;
      %n;
   end
   % Stage 0, BASIC QUALITY CHECKING
   properties (Hidden = true, Dependent = true)
      ST0MinSample;
   end
   properties (Dependent = true)
      ST0Valid;
   end
   % STAGE 1: LOCATION BASE TESTS
   properties
      ST1P = ones(1,7);
      ST1Stat = zeros(1,7);
   end
   properties (Hidden = true, Dependent = true)
      ST1PT;
      ST1MinHits; 
   end
   properties (Dependent = true)
      ST1Code;
      ST1Label;
      ST1TestInd;
      ST1Valid;
   end
   % STAGE 2: SNPBRIM BASE TESTS
   properties
      ST2F = zeros(1,7);
      ST2FP = ones(1,7);
      ST2G = zeros(1,7);
      ST2GP = ones(1,7);
      ST2FFold = zeros(4,7);
      ST2FPFold = ones(4,7);
      ST2GFold = zeros(4,7);
      ST2GPFold = ones(4,7);
   end
   properties (Hidden = true, Dependent = true)
      ST2PFT;
      ST2PGT;
      ST2PFFoldT;
      ST2PGFoldT;
      ST2NrFoldT;
      ST2MinHits;
   end
   properties (Dependent = true)
      ST2FMedFold;
      ST2FMaxFold;
      ST2FMinFold;
      ST2GMedFold;
      ST2GMaxFold;
      ST2GMinFold;
      ST2FCode;
      ST2GCode;
      ST2FFoldCode;
      ST2GFoldCode;
      ST2FTestInd;
      ST2GTestInd;
      ST2FLabel;
      ST2GLabel;
      ST2Valid;
      ST2FValid;
      ST2GValid;
   end
   % STAGE 3: COVARIATES PATTERN TEST
   properties
      ST3TestInd = 0;
      ST3F = 0;
      ST3FP = 1;
      ST3G = 0;
      ST3GP = 1;
      ST3FFold = zeros(4,1);
      ST3FPFold = ones(4,1);
      ST3GFold = zeros(4,1);
      ST3GPFold = ones(4,1);
   end
   properties (Hidden = true, Dependent = true)
      ST3PFT;
      ST3PGT;
      ST3PFFoldT;
      ST3PGFoldT;
      ST3NrFoldT;
      ST3PStrength;
   end
   properties (Dependent = true)
      ST3FMedFold;
      ST3FMaxFold;
      ST3FMinFold;
      ST3GMedFold;
      ST3GMaxFold;
      ST3GMinFold;
      ST3Label;
      ST3Valid;
      ST3FValid;
      ST3GValid;
   end
   % STAGE 4: EXTRACTING PREDICTION MODEL
   properties (Dependent = true)
       ST4TestInd;
       ST4Label;
   end
   properties (Hidden = true)
       ST4MOD = [];
       %ST4Ref = [];
       ST4RedShape = false;
   end
   properties (Hidden = true, Dependent = true)
       %ST4ShapeSpace;
       %ST4RefScan;
   end
   methods % CONSTRUCTOR
        function obj = SNP(varargin)
            obj = obj@superClass(varargin{:});         
        end
   end
   methods % GENERAL: GETTING/SETTING
       function out = get.ParCont(obj)
           out = obj.ParCont;
           if ~superClass.isH(out), out = []; end
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
       function out = get.nrAA(obj)
           if isempty(obj.GT), out = 0; return; end
           out = length(obj.AAind);
       end
       function out = get.nrAB(obj)
           if isempty(obj.GT), out = 0; return; end
           out = length(obj.ABind);
       end
       function out = get.nrBB(obj)
           if isempty(obj.GT), out = 0; return; end
           out = length(obj.BBind);
       end
       function out = get.Balances(obj)
           out = zeros(1,6);
           if isempty(obj.GT), return; end
           out(1) = SNP.getBalance(obj.nrAA+obj.nrAB,obj.nrBB);% Recessive
           out(2) = SNP.getBalance(obj.nrBB+obj.nrAB,obj.nrAA);% Dominant
           out(3) = SNP.getBalance(obj.nrAA+obj.nrBB,obj.nrAB);% Over Dominant
           out(4) = SNP.getBalance(obj.nrAA,obj.nrBB);% AA versus BB
           out(5) = SNP.getBalance(obj.nrAA,obj.nrAB);% AA versus AB
           out(6) = SNP.getBalance(obj.nrAB,obj.nrBB);% AB versus BB
       end
       function out = get.Valid(obj)
          out = obj.ST0Valid*obj.ST1Valid*obj.ST2Valid; 
       end
       function out = get.OuterFold(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.OuterFold;
       end
       function out = get.InnerFold(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.InnerFold;
       end
       function out = get.UseRedShape(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.UseRedShape;
       end
       function out = get.MinFoldSampleSize(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.MinFoldSampleSize;
       end
       function out = get.MaxIterations(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.MaxIterations;
       end
       function out = get.Htest(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.Htest;
       end
       function out = get.RIPNormalize(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.RIPNormalize;
       end
       function out = get.StopCorr(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.StopCorr;
       end
       function out = get.RIPStatPerm(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.RIPStatPerm;
       end
       function out = get.RegRuns(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.RegRuns;
       end
       function out = get.RegSampling(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.RegSampling;
       end
       function out = get.RegBalance(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.RegBalance;
       end
       function out = get.RegBalanceTH(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.RegBalanceTH;
       end
       function out = get.RegBalanceMethod(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.RegBalanceMethod;
       end
       function out = get.GroupDef(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.GroupDef;
       end
       function out = get.GTRIP(obj)
           if isempty(obj.ST4MOD), out = []; return; end
           out = obj.ST4MOD.CompRIP;
       end
   end
   methods % STAGE 0: GETTING/SETTING
       function out = get.ST0Valid(obj)
          if length(obj.AAind)<obj.ST0MinSample, out = 0; return; end
          if length(obj.ABind)<obj.ST0MinSample, out = 0; return; end
          if length(obj.BBind)<obj.ST0MinSample, out = 0; return; end
          out = 1;
       end
       function out = get.ST0MinSample(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST0MinSample;
       end
   end
   methods % STAGE 1: GETTING/SETTING  
       function out = get.ST1Code(obj)
                out = obj.ST1P<=obj.ST1PT;
       end
       function out = get.ST1Label(obj)
                out = retrieveLabel(obj,1,[]);
       end
       function out = get.ST1TestInd(obj)
                out = retrieveTestInd(obj,1,[]);
       end
       function out = get.ST1Valid(obj)
                out = isValid(obj,1);
       end
       function out = get.ST1PT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST1PT;
       end
       function out = get.ST1MinHits(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST1MinHits;
       end
   end
   methods % STAGE 2: GETTING/SETTING
       function out = get.ST2FMedFold(obj)
                out = median(obj.ST2FFold);         
       end
       function out = get.ST2FMaxFold(obj)
                out = max(obj.ST2FFold);         
       end
       function out = get.ST2FMinFold(obj)
                out = min(obj.ST2FFold);         
       end
       function out = get.ST2GMedFold(obj)
                out = median(obj.ST2GFold);         
       end
       function out = get.ST2GMaxFold(obj)
                out = max(obj.ST2GFold);         
       end
       function out = get.ST2GMinFold(obj)
                out = min(obj.ST2GFold);         
       end
       function out = get.ST2FCode(obj)
                out = obj.ST2FP<=obj.ST2PFT;
                out = (out.*obj.ST2FFoldCode>=obj.ST2NrFoldT);
       end
       function out = get.ST2FFoldCode(obj)
                out = sum(obj.ST2FPFold<=obj.ST2PFFoldT);
       end
       function out = get.ST2GCode(obj)
                out = obj.ST2GP<=obj.ST2PGT;
                out = (out.*obj.ST2GFoldCode>=obj.ST2NrFoldT);
       end
       function out = get.ST2GFoldCode(obj)
                out = sum(obj.ST2GPFold<=obj.ST2PGFoldT);
       end
       function out = get.ST2FLabel(obj)
                out = retrieveLabel(obj,2,'F');
       end
       function out = get.ST2FTestInd(obj)
                out = retrieveTestInd(obj,2,'F');
       end
       function out = get.ST2GLabel(obj)
                out = retrieveLabel(obj,2,'G');
       end
       function out = get.ST2GTestInd(obj)
                out = retrieveTestInd(obj,2,'G');
       end
       function out = get.ST2Valid(obj)
                out = obj.ST2FValid*obj.ST2GValid;
       end
       function out = get.ST2FValid(obj)
                out = isValid(obj,2,'F');
       end
       function out = get.ST2GValid(obj)
                out = isValid(obj,2,'G');
       end
       function out = get.ST2PFFoldT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST2PFFoldT;
       end
       function out = get.ST2PFT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST2PFT;
       end
       function out = get.ST2PGFoldT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST2PGFoldT;
       end
       function out = get.ST2PGT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST2PGT;
       end
       function out = get.ST2NrFoldT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST2NrFoldT;
       end
       function out = get.ST2MinHits(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST2MinHits;
       end
   end
   methods % STAGE 3: GETTING/SETTING
       function out = get.ST3Label(obj)
           out = retrieveLabel(obj,3);
       end
       function out = get.ST3FMedFold(obj)
                out = median(obj.ST3FFold);         
       end
       function out = get.ST3FMaxFold(obj)
                out = max(obj.ST3FFold);         
       end
       function out = get.ST3FMinFold(obj)
                out = min(obj.ST3FFold);         
       end
       function out = get.ST3GMedFold(obj)
                out = median(obj.ST3GFold);         
       end
       function out = get.ST3GMaxFold(obj)
                out = max(obj.ST3GFold);         
       end
       function out = get.ST3GMinFold(obj)
                out = min(obj.ST3GFold);         
       end
       function out = get.ST3Valid(obj)
                out = max([obj.ST3FValid obj.ST3GValid]);
       end
       function out = get.ST3FValid(obj)
                out = (obj.ST3FP<=obj.ST3PFT)*(sum(obj.ST3FPFold<=obj.ST3PFFoldT)>=obj.ST3NrFoldT);
       end
       function out = get.ST3GValid(obj)
                out = (obj.ST3GP<=obj.ST3PGT)*(sum(obj.ST3GPFold<=obj.ST3PGFoldT)>=obj.ST3NrFoldT);
       end
       function out = get.ST3PFFoldT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST3PFFoldT;
       end
       function out = get.ST3PFT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST3PFT;
       end
       function out = get.ST3PGFoldT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST3PGFoldT;
       end
       function out = get.ST3PGT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST3PGT;
       end
       function out = get.ST3NrFoldT(obj)
           if isempty(obj.ParCont), out = []; return; end
           out = obj.ParCont.ST3NrFoldT;
       end
       function out = get.ST3PStrength(obj)
                out = SNP.pfast([obj.ST3FP obj.ST3GP]);
       end
   end
   methods % STAGE 4: GETTING/SETTING
       function out = get.ST4TestInd(obj)
                out = obj.ST3TestInd(1);
       end
       function out = get.ST4Label(obj)
                out = obj.ST3Label;
       end
   end
   methods % RUNNING STAGES
       function obj = runStage1(obj,t,DepVar,parcont)
           if nargin==4, obj.ParCont = parcont;end
           for i=1:1:7
               try
                  GG = getGG(obj,i);
                  if ~(i==4) % TWO GROUP TESTS
                     [stat,p] = SNP.DStatTest(DepVar(GG==-1,:),DepVar(GG==1,:),t);  
                  else % THREE GROUP TEST / ADDITIVE
                     [stat,p] = SNP.FStatTest(DepVar(GG==-1,:),DepVar(GG==0,:),DepVar(GG==1,:),t);
                  end
                  obj.ST1Stat(i) = stat;obj.ST1P(i) = p;
               catch
                  obj.ST1Stat(i) = nan;obj.ST1P(i) = nan;
               end
           end
       end
       function obj = runStage2(obj,DepVar,parcont)
           if nargin==3, obj.ParCont = parcont;end
           FFold = nan*zeros(obj.OuterFold,7);FPFold = nan*zeros(obj.OuterFold,7);
           GFold = nan*zeros(obj.OuterFold,7);GPFold = nan*zeros(obj.OuterFold,7);
           Fouter = getST2OuterFolds(obj,DepVar);% you want to use the same folds in each test, to be able to compare them
           for i=1:1:7
               try
                 out = SNPBRIM(obj,i,DepVar,[],Fouter);
                 obj.ST2FP(i) = out.FP;obj.ST2F(i) = out.F;
                 obj.ST2GP(i) = out.GP;obj.ST2G(i) = out.G;
                 FFold(:,i) = out.FoldF(:);FPFold(:,i) = out.FoldFP(:);
                 GFold(:,i) = out.FoldG(:);GPFold(:,i) = out.FoldGP(:);
               catch %#ok<*CTCH>
                 obj.ST2FP(i) = nan;obj.ST2F(i) = nan;
                 obj.ST2GP(i) = nan;obj.ST2G(i) = nan;
                 FFold(:,i) = nan*zeros(obj.OuterFold,1);FPFold(:,i) = nan*zeros(obj.OuterFold,1);
                 GFold(:,i) = nan*zeros(obj.OuterFold,1);GPFold(:,i) = nan*zeros(obj.OuterFold,1);
               end
           end
           obj.ST2FFold = FFold;obj.ST2FPFold = FPFold;
           obj.ST2GFold = GFold;obj.ST2GPFold = GPFold;
       end
       function obj = prepareStage3(obj,stages)
           disp('  ');
           disp(obj.RS)
           disp('  ');
           disp([obj.TestLabels{1} '     ' obj.TestLabels{2} '     ' obj.TestLabels{3} '     ' obj.TestLabels{4} '       ' obj.TestLabels{5} '       ' obj.TestLabels{6} '       ' obj.TestLabels{7}]); 
           disp('--------------------------------------------------------------------------------------');
           for i=1:1:length(stages)
               switch stages(i)
                   case 1
                       disp('SUGGESTION STAGE 1: LOCATION TESTING');
                       tmp = zeros(2,7);tmp(1,:) = obj.ST1P;tmp(2,:) = obj.ST1Code;
                       disp(num2str(tmp));
                       disp('  ');
                       disp(['TEST ID: ' num2str(obj.ST1TestInd(1)) ' ==> ' obj.ST1Label]);
                   case 2
                       disp('SUGGESTION STAGE 2: F STATISTIC');
                       tmp = zeros(5,7);
                       tmp(1,:) = obj.ST2FMedFold;tmp(2,:) = obj.ST2FMaxFold;tmp(3,:) = obj.ST2FP;
                       tmp(4,:) = obj.ST2FFoldCode;tmp(5,:) = obj.ST2FCode;
                       disp(num2str(tmp));
                       disp('  ');
                       disp(['TEST ID: ' num2str(obj.ST2FTestInd(1)) ' ==> ' obj.ST2FLabel]);
                       disp('--------------------------------------------------------------------------------------');
                       disp('  ');
                       disp('SUGGESTION STAGE 2: G STATISTIC');
                       tmp = zeros(5,7);
                       tmp(1,:) = obj.ST2GMedFold;tmp(2,:) = obj.ST2GMaxFold;tmp(3,:) = obj.ST2GP;
                       tmp(4,:) = obj.ST2GFoldCode;tmp(5,:) = obj.ST2GCode;
                       disp(num2str(tmp));
                       disp('  ');
                       disp(['TEST ID: ' num2str(obj.ST2GTestInd(1)) ' ==> ' obj.ST2GLabel]);
                   case 3
                       disp('SUGGESTION STAGE 3: NGBRIM');
               end
               disp('--------------------------------------------------------------------------------------');
               disp('  ');
           end
           disp('nrAA nrAB nrBB');
           tmp = zeros(1,3);tmp(1) = obj.nrAA;tmp(2) = obj.nrAB; tmp(3) = obj.nrBB;
           disp(num2str(tmp));
           disp('--------------------------------------------------------------------------------------');
           disp('  ');
           disp('0      : NO EFFECT');disp('1      : RECESSIVE');disp('2      : DOMINANT');
           disp('3      : OVER DOMINANT');disp('4      : ADDITIVE');disp('5      : AA BB');
           disp('6      : AA AB');disp('7      : AB BB');disp('8      : NON ADDITVE');
           disp('9      : ADDITIVE AA AB');disp('10      : ADDITIVE AB BB');
           str2 = input('Enter DESIRED inheritance pattern [0-10]','s');
           obj.ST3TestInd = [str2double(str2) 1];           
       end
       function obj = runStage3(obj,DepVar,Cov,parcont)
           if nargin==4, obj.ParCont = parcont;end
           try
               out = SNPBRIM(obj,obj.ST3TestInd(1),DepVar,Cov,[]);
               obj.ST3FP = out.FP;obj.ST3F = out.F;
               obj.ST3GP = out.GP;obj.ST3G = out.G;
               obj.ST3FFold = out.FoldF(:);obj.ST3FPFold = out.FoldFP(:);
               obj.ST3GFold = out.FoldG(:);obj.ST3GPFold = out.FoldGP(:);
           catch
               obj.ST3FP = nan;obj.ST3F = nan;
               obj.ST3GP = nan;obj.ST3G = nan;
               obj.ST3FFold = nan*zeros(obj.OuterFold,1);obj.ST3FPFold = nan*zeros(obj.OuterFold,1);
               obj.ST3GFold = nan*zeros(obj.OuterFold,1);obj.ST3GPFold = nan*zeros(obj.OuterFold,1);
           end         
       end
       function obj = runStage4(obj,BaseCont,UseRedShape,parcont)
           try
               if nargin==4, obj.ParCont = parcont; end
               % Determine Shape Space to use
                 if UseRedShape
                    DepVar = BaseCont.RedDepVar;Cov = [];
                    %OrigDepVar = BaseCont.OrigRedDepVar;
                    RefCoeff = BaseCont.ST4Ref.RedCoeff;
                    %ShapeSpace = BaseCont.RedShapeSpace;
                 else
                    DepVar = BaseCont.DepVar;
                    %OrigDepVar = BaseCont.OrigDepVar;
                    RefCoeff = BaseCont.ST4Ref.Coeff;
                    %ShapeSpace = BaseCont.ShapeSpace;
                    switch BaseCont.RedShapeSpaceType
                        case 'X'
                           Cov = [BaseCont.Cov, BaseCont.GB];
                        case 'RIP'
                           Cov = [BaseCont.CovRIP, BaseCont.GBRIP];
                    end
                 end
               % Remember shape space used  
                 obj.ST4RedShape = UseRedShape;
               % Run SNPBRIM 
                 obj.ST4MOD = SNPBRIM(obj,obj.ST4TestInd,DepVar,Cov);
                 obj.ST4MOD = rmfield(obj.ST4MOD,{'nrAdded' 'DepVarAdd' 'IndVarAdd' 'GMAdd' 'BootIter' 'FoldF' 'FoldFP' 'FoldG' 'FoldGP' 'DepVar' 'IndVar','Grp','H'});
               % Setting up reference scan 
                 RIPref = SNP.updateRIP(RefCoeff,obj.ST4MOD.M);
                 RIPref = SNP.normalizeRIP(RIPref,obj.ST4MOD.M,obj.ST4MOD.Var);
                 obj.ST4MOD.RIPRef = RIPref;
               % Complete population
                 %RIP = SNP.updateRIP(DepVar,obj.ST4MOD.M);
                 %RIP = SNP.normalizeRIP(RIP,obj.ST4MOD.M,obj.ST4MOD.Var);
                 %obj.ST4MOD.CompRIP = RIP;  
               % Population used to model things
                 RIP = SNP.updateRIP(DepVar(obj.ST4MOD.UsedInd,:),obj.ST4MOD.M);
                 RIP = SNP.normalizeRIP(RIP,obj.ST4MOD.M,obj.ST4MOD.Var);
                 obj.ST4MOD.RIP = RIP;
                 [obj.ST4MOD.RIPF,~,obj.ST4MOD.RIPG,~,obj.ST4MOD.RIPGT] = SNP.performRIPStat(obj.ST4MOD.GG,RIP,obj.ST4MOD.Var,0);
               % Obtaining regressions to contruct morphs
                 IndVar = [Cov(obj.ST4MOD.UsedInd,:) RIP'];
                 [obj.ST4MOD.MMorphs,obj.ST4MOD.MIMorphs,~] = getRegression(IndVar,DepVar(obj.ST4MOD.UsedInd,:),(1:size(IndVar,2)));
                 %%[obj.ST4MOD.MMorphsOrig,obj.ST4MOD.MIMorphsOrig,obj.ST4MOD.R2MorphsOrig] = getRegression(IndVar,OrigDepVar(obj.ST4MOD.UsedInd,:),(1:size(IndVar,2)));
                 %popShape = repmat(ShapeSpace.AvgVec,1,length(obj.ST4MOD.UsedInd)) + ShapeSpace.EigVec*DepVar(obj.ST4MOD.UsedInd,:)';
                 %[~,~,R2] = getRegression(IndVar,popShape',size(IndVar,2));
                 %R2 = reshape(R2,3,ShapeSpace.Average.nrV);
                 %for i=1:1:3
                 %    R2(i,:) = smoothFunction(ShapeSpace.Average,R2(i,:)',2,'functiondistance');
                 %end
                 %obj.ST4MOD.R2 = R2;
                 %obj.ST4MOD.R2LM = sum(obj.ST4MOD.R2);
               % Obtaining group memberships and additional needed info
                 GTInfo = cell(1,3);% first is AA, second is AB and third is BB
                 GTInfo{1}.GT = -1;GTInfo{2}.GT = 0;GTInfo{3}.GT = 1;
                 switch obj.ST4TestInd
                     case 1 % RECESSIVE
                         GTInfo{1}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{1}.Sign = -1;
                         GTInfo{1}.Opp = 3;
                         GTInfo{2}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{2}.Sign = -1;
                         GTInfo{2}.Opp = 3;
                         GTInfo{3}.Ind = find(obj.ST4MOD.GG==1);
                         GTInfo{3}.Sign = 1;
                         GTInfo{3}.Opp = 1;
                     case 2 % DOMINANT
                         GTInfo{1}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{1}.Sign = -1;
                         GTInfo{1}.Opp = 3;
                         GTInfo{2}.Ind = find(obj.ST4MOD.GG==1);
                         GTInfo{2}.Sign = 1;
                         GTInfo{2}.Opp = 1;
                         GTInfo{3}.Ind = find(obj.ST4MOD.GG==1);
                         GTInfo{3}.Sign = 1;
                         GTInfo{3}.Opp = 1;
                     case 3 % OVER DOMINANT
                         GTInfo{1}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{1}.Sign = -1;
                         GTInfo{1}.Opp = 2;
                         GTInfo{2}.Ind = find(obj.ST4MOD.GG==1);
                         GTInfo{2}.Sign = 1;
                         GTInfo{2}.Opp = 1;
                         GTInfo{3}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{3}.Sign = -1;
                         GTInfo{3}.Opp = 2;
                     case 4 % ADDITIVE
                         GTInfo{1}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{1}.Sign = -1;
                         GTInfo{1}.Opp = [2 3];
                         GTInfo{2}.Ind = find(obj.ST4MOD.GG==0);
                         GTInfo{2}.Sign = 0;
                         GTInfo{2}.Opp = [1 3];
                         GTInfo{3}.Ind = find(obj.ST4MOD.GG==1);
                         GTInfo{3}.Sign = 1;
                         GTInfo{3}.Opp = [1 2];
                     case 5 % AA BB
                         GTInfo{1}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{1}.Sign = -1;
                         GTInfo{1}.Opp = 3;
                         GTInfo{2}.Ind = [];
                         GTInfo{2}.Sign = 0;
                         GTInfo{2}.Opp = [];
                         GTInfo{3}.Ind = find(obj.ST4MOD.GG==1);
                         GTInfo{3}.Sign = 1;
                         GTInfo{3}.Opp = 1;
                     case 6 % AA AB
                         GTInfo{1}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{1}.Sign = -1;
                         GTInfo{1}.Opp = 2;
                         GTInfo{2}.Ind = find(obj.ST4MOD.GG==1);
                         GTInfo{2}.Sign = 1;
                         GTInfo{2}.Opp = 1;
                         GTInfo{3}.Ind = [];
                         GTInfo{3}.Sign = 0;
                         GTInfo{3}.Opp = [];
                     case 7 % AB BB
                         GTInfo{1}.Ind = [];
                         GTInfo{1}.Sign = 0;
                         GTInfo{1}.Opp = [];
                         GTInfo{2}.Ind = find(obj.ST4MOD.GG==-1);
                         GTInfo{2}.Sign = -1;
                         GTInfo{2}.Opp = 3;
                         GTInfo{3}.Ind = find(obj.ST4MOD.GG==1);
                         GTInfo{3}.Sign = 1;
                         GTInfo{3}.Opp = 2;
                     otherwise
                         GTInfo{1}.Ind = [];
                         GTInfo{1}.Sign = 0;
                         GTInfo{1}.Opp = [];
                         GTInfo{2}.Ind = [];
                         GTInfo{2}.Sign = 0;
                         GTInfo{2}.Opp = [];
                         GTInfo{3}.Ind = [];
                         GTInfo{3}.Sign = 0;
                         GTInfo{3}.Opp = [];
                 end
                 GTInfo{1}.Freq = length(GTInfo{1}.Ind)/length(find(~isnan(obj.ST4MOD.GG)));
                 GTInfo{2}.Freq = length(GTInfo{2}.Ind)/length(find(~isnan(obj.ST4MOD.GG)));
                 GTInfo{3}.Freq = length(GTInfo{3}.Ind)/length(find(~isnan(obj.ST4MOD.GG)));
               % Extracting group RIP Statistics  
                 for i=1:1:3
                     GTInfo{i}.Nr = length(GTInfo{i}.Ind);
                     if GTInfo{i}.Nr>0
                        GTInfo{i}.RIP = obj.ST4MOD.RIP(GTInfo{i}.Ind);
                        [GTInfo{i}.AVG,GTInfo{i}.STD,GTInfo{i}.MED,GTInfo{i}.MAD,GTInfo{i}.O] = SNP.extractRobustRIPStatistics(GTInfo{i}.RIP,6);                     
                     end
                 end
                 obj.ST4MOD.GTInfo = GTInfo;
                 obj.ST4MOD.Freq = [GTInfo{1}.Freq GTInfo{2}.Freq GTInfo{3}.Freq];
               % Creating average morphs
                 %obj.ST4MOD.Atyps = nan*zeros(1,3);
                 %for i=1:1:3
                 %    if obj.ST4MOD.GTInfo{i}.Nr==0, continue; end
                 %    [obj.ST4MOD.GTInfo{i}.C,obj.ST4MOD.GTInfo{i}.Atyp,obj.ST4MOD.GTInfo{i}.Morph] = createMorph(obj,BaseCont,obj.ST4MOD.GTInfo{i}.GT,0);
                 %    obj.ST4MOD.Atyps(i) = obj.ST4MOD.GTInfo{i}.Atyp;
                 %end
               % Determination of less typical group
                 %obj.ST4MOD.MAtyp = nanmax(obj.ST4MOD.Atyps);
                 %obj.ST4MOD.mAtyp = nanmin(obj.ST4MOD.Atyps);
                 %obj.ST4MOD.BalAtyp = obj.ST4MOD.mAtyp/obj.ST4MOD.MAtyp;   
           catch
               obj.ST4MOD = [];
           end
       end
   end
   methods % Interfacing Functions
       function out = getGG(obj,TestInd)
           out = obj.GT;% Initialize
           switch TestInd
               case 0 % NO EFFECT
                   out = nan*out;
               case 1 % test 1: AAAB BB (RECESSIVE)
                   out(union(obj.AAind,obj.ABind)) = -1;
                   out(obj.BBind) = 1;
               case 2 % test 2: AA ABBB (DOMINANT)
                   out(obj.AAind) = -1;
                   out(union(obj.BBind,obj.ABind)) = 1;
               case 3 % test 3: AABB AB (OVER DOMINANT)
                   out(union(obj.AAind,obj.BBind)) = -1;
                   out(obj.ABind) = 1;
               case 4 % test 7 AA AB BB (ADDITIVE)
                   return;
               case 5 % test 4: AA BB
                   out(obj.ABind) = nan;
               case 6 % test 5: AA AB
                   out(obj.ABind) = 1;
                   out(obj.BBind) = nan;
               case 7 % test 6: AB BB
                   out(obj.ABind) = -1;
                   out(obj.AAind) = nan;
               case 8 % NON ADDITIVE
                   return;
               case 9 % ADDITIVE AA AB
                   out(obj.AAind) = -1;
                   out(obj.BBind) = 0;
                   out(obj.ABind) = 1;
               case 10 % ADDITIVE AB BB
                   out(obj.ABind) = -1;
                   out(obj.AAind) = 0;
                   out(obj.BBind) = 1;
           end       
       end
       function out = isValid(obj,stage,stat)
           switch stage
               case 1
                   code = obj.ST1Code;
                   MH = obj.ST1MinHits;
               case 2
                   switch stat
                       case 'F' % Using F statistics 
                            code = obj.ST2FCode;
                       case 'G' % Using G statistics 
                            code = obj.ST2GCode;
                   end
                   MH = obj.ST2MinHits;
               case 3
                   code = obj.ST3Code;
                   MH = 1;
           end
           if sum(code)>=MH
              out = 1;
           else
              out = 0;
           end
       end
       function out = retrieveLabel(obj,stage,stat)
           switch stage
               case 1
                   test = obj.ST1TestInd;
               case 2
                   switch stat
                       case 'F'
                           test = obj.ST2FTestInd;
                       case 'G'
                           test = obj.ST2GTestInd;
                   end 
               case 3
                   test = obj.ST3TestInd;
           end
           switch test(1)
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
           %switch test(2)
           %    case 0
           %        out = [out ' / WEAK'];
           %    case 1
           %        out = [out ' / STRONG'];
           %end
       end
       function out = retrieveTestInd(obj,stage,stat)
           out = [0 1];
           switch stage
               case 1
                   code = obj.ST1Code;P = obj.ST1P;
               case 2
                   switch stat
                       case 'F' % Using F statistics
                            code = obj.ST2FCode;P = obj.ST2FP;
                       case 'G' % Using G statistics
                            code = obj.ST2GCode;P = obj.ST2GP;
                   end
           end
           if sum(code)==0, return; end % 0 0 0 0 0 0
           Tests = SNP.interpCodeP1(code(1:4));impos = SNP.interpCodeP2(code(5:7),P(5:7)>0.70);
           Tests = setdiff(Tests,impos);% remove the impossible combinations
           if isempty(Tests), out = [0 0]; return; end % no effect
           if length(Tests)==1, out = [Tests 1]; return; end % single test option
           [~,ind] = min(P(setdiff(Tests,[8 9 10])));out(1) = Tests(ind);out(2) = 0; % mulitple test options
           %out(1) = min(Tests);out(2) = 0;
       end
       function out = getST2OuterFolds(obj,DepVar,parcont)
                if nargin==3, obj.ParCont = parcont; end
                ParInfo = obj.ParCont;
                out = DOB_SCV(ParInfo.OuterFold,DepVar,obj.GT);
       end
       function out = getShapeEffects(obj,BaseCont)
           % TO BE IMPLEMENTED
           tmp = BaseCont.ST4Ref;
           IndVar = [Cov(obj.ST4MOD.UsedInd,:) RIP'];
           [~,R2] = getRegression(IndVar,popShape,size(IndVar,2));
           R2 = reshape(R2,3,obj.ST4Ref.Scan.nrV);
           for i=1:1:3
               R2(i,:) = smoothFunction(obj.ST4Ref.Scan,R2(i,:)',2,'functiondistance');
           end
           obj.ST4MOD.R2 = R2;
           obj.ST4MOD.R2LM = sum(obj.ST4MOD.R2);
           obj.ST4Ref.Scan.Value = obj.ST4MOD.R2LM;
           out = [];
       end
       function obj = reduceSamples(obj,index)
           if nargout==1, obj = clone(obj); end
           obj.GT = obj.GT(index);
           obj.ST4MOD = [];
           try
%             obj.ST4MOD.CompRIP = obj.ST4MOD.CompRIP(index);
%             [ism,loc] = ismember(obj.ST4MOD.UsedInd,index);
%             obj.ST4MOD.RIP = obj.ST4MOD.RIP(find(ism));
%             obj.ST4MOD.UsedInd = loc(find(ism));
%             obj.ST4MOD.GM = obj.ST4MOD.GM(find(ism));
           catch
           end
       end
   end
   methods % NextGeneration BRIM
       function out = SNPBRIM(obj,TestInd,DepVar,Cov,Fouter)
              if nargin<5, Fouter = [];end
           % Make sure to extract the Parameter Container, for when parfor is executed
              ParInfo = obj.ParCont;
           % To bootstrap or not to bootstrap that is the question
              Bootstrap = true;if ParInfo.MaxIterations == 0, Bootstrap = false;end
           % Defining Group Memberhip
              GG = getGG(obj,TestInd); 
              ind = find(~isnan(GG));GG = GG(ind);
              switch ParInfo.GroupDef
                  case 'GG' % Group membership is defined by the test performed
                      GM = GG;
                  case 'GT' % Group membership is defined by the original genotypes
                      GM = obj.GT(ind);
              end
              out.UsedInd = ind;% store the samples used
              out.GG = GG;
              if ~isempty(Fouter), Fouter = Fouter(ind);end
           % Defining Independent and Dependent Variables
              DepVar = DepVar(ind,:);
              if isempty(Cov);
                 IndVar = GG;
              else
                 IndVar = [Cov(ind,:) GG];
              end
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
              [~,tmp] = min(var.nrWEl);var.PosClass = var.El(tmp);% needed for ROC analysis
            % Defining group info
              grp.El = unique(GM);grp.nrEl = length(grp.El);
            % Upsampling of groups if required
              minreq = ParInfo.OuterFold*ParInfo.MinFoldSampleSize(1);
              nrAdded = 0;DepVarAdd = [];IndVarAdd = [];GMAdd = [];
              grp.UpSampled = zeros(1,grp.nrEl);
              if isempty(Fouter)% Cannot upsample with predefined outer folds
                  for i=1:1:grp.nrEl
                     indG = find(GM==grp.El(i));
                     if length(indG)>=minreq; continue; end
                     nrSyn = minreq-length(indG);grp.UpSampled(i) = nrSyn;
                     [DepVarTmp,IndVarTmp] = SNP.foldUpSample(DepVar(indG,:),IndVar(indG,:),nrSyn);
                     nrAdded = nrAdded + nrSyn;
                     DepVarAdd = [DepVarAdd;DepVarTmp];IndVarAdd = [IndVarAdd;IndVarTmp];GMAdd = [GMAdd; grp.El(i)*ones(nrSyn,1)]; %#ok<*AGROW>
                  end
              end
              out.nrAdded = nrAdded;out.DepVar = DepVar;out.IndVar = IndVar;out.GM = GM;
              out.DepVarAdd = DepVarAdd;out.IndVarAdd = IndVarAdd;out.GMAdd = GMAdd;out.Var = var;out.Grp = grp;
              DepVar = [DepVar;DepVarAdd];IndVar = [IndVar;IndVarAdd];GM = [GM;GMAdd];
              [n,dim] = size(DepVar);
            % Outer Partitioning of Data
              if isempty(Fouter), Fouter = DOB_SCV(ParInfo.OuterFold,DepVar,GM);end
              FoldResults = cell(1,ParInfo.OuterFold);
              minreq = ParInfo.InnerFold*ParInfo.MinFoldSampleSize(2);
              parfor fo=1:ParInfo.OuterFold
                 % Extracting training and testing 
                  FoTestInd = find(Fouter==fo);FoTrInd = setdiff((1:n),FoTestInd);nrFoTr = length(FoTrInd);
                  FoldResults{fo}.TestInd = FoTestInd;FoldResults{fo}.TrInd = FoTrInd;
                  FoTrIndVar = IndVar(FoTrInd,:);FoTrDepVar = DepVar(FoTrInd,:);FoTrGM = GM(FoTrInd);
                  FoldResults{fo}.IndVar = IndVar(FoTestInd,end);
                  FoTestDepVar = DepVar(FoTestInd,:);
                  FoldResults{fo}.DepVar = FoTestDepVar;
                  BootProgress = zeros(1,ParInfo.MaxIterations);
                 % Initialize Bootstrapping
                  FoldResults{fo}.UpSampled = zeros(1,grp.nrEl);
                  ContBoot = true;counter = 0;
                  while ContBoot && Bootstrap>0
                    % Keeping track of iterations  
                     counter = counter + 1; 
                    % Upsampling of groups if required
                     nrAdded = 0;DepVarAdd = [];IndVarAdd = [];GMAdd = [];
                     for i=1:1:grp.nrEl
                       indG = find(FoTrGM==grp.El(i));
                       if length(indG)>=minreq; continue; end
                       nrSyn = minreq-length(indG);nrAdded = nrAdded + nrSyn;FoldResults{fo}.UpSamples(i) = nrSyn;
                       [DepVarTmp,IndVarTmp] = SNP.foldUpSample(FoTrDepVar(indG,:),FoTrIndVar(indG,:),nrSyn);
                       GMAdd = [GMAdd;grp.El(i)*ones(nrSyn,1)];DepVarAdd = [DepVarAdd;DepVarTmp];IndVarAdd = [IndVarAdd;IndVarTmp];
                     end
                     FoTrDepVar = [FoTrDepVar;DepVarAdd];FoTrIndVar = [FoTrIndVar;IndVarAdd];FoTrGM = [FoTrGM;GMAdd];
                   % seperate data into inner folds
                     Finner = DOB_SCV(ParInfo.InnerFold,FoTrDepVar,FoTrGM);
                     TMPIndVar = FoTrIndVar(:,end);% Allocate Memory, only last collumn for SNP
                     for fi=1:ParInfo.InnerFold % Fi...
                         FiTestInd = find(Finner==fi);
                         FiTrInd = setdiff(1:nrFoTr+nrAdded,FiTestInd);
                         [M,H] = robustPLSRegression(obj,FoTrIndVar(FiTrInd,:),FoTrDepVar(FiTrInd,:),FoTrGM(FiTrInd),ParInfo);% Compute BootSampled Regression
                         if ParInfo.Htest, M = H.*M;end% Only keep significant regression coefficients
                         rip = SNP.updateRIP(FoTrDepVar(FiTestInd,:),M)';% Get RIP scores
                         if ParInfo.RIPNormalize, rip = SNP.normalizeRIP(rip,M,var);end% rescale RIP scores
                         TMPIndVar(FiTestInd) = rip;% Store updated RIP scores, always in the last collumn
                     end
                   % Removing any addded samples 
                     TMPIndVar = TMPIndVar(1:nrFoTr);
                     FoTrIndVar = FoTrIndVar(1:nrFoTr,:);
                     FoTrDepVar = FoTrDepVar(1:nrFoTr,:);
                     FoTrGM = FoTrGM(1:nrFoTr);
                   % monitoring progress of Booting
                     tmpC = corrcoef(FoTrIndVar(:,end),TMPIndVar);
                     BootProgress(counter) = tmpC(1,2);
                     if BootProgress(counter)>=ParInfo.StopCorr, ContBoot = false; end
                     if counter >= ParInfo.MaxIterations, ContBoot = false; end
                     FoTrIndVar(:,end) = TMPIndVar;% Update Training IndVar for next round
                  end
                  if counter>0, FoldResults{fo}.BootProgress = BootProgress(1:counter);end
                  FoldResults{fo}.BootIter = counter;
                 % Extracting final M from outer training data
                  [M,H] = robustPLSRegression(obj,FoTrIndVar,FoTrDepVar,FoTrGM,ParInfo);
                  if ParInfo.Htest, M = H.*M;end
                  FoldResults{fo}.M = M;FoldResults{fo}.H = H;
                 % Evaluate outer test data
                  rip = updateRIP(FoTestDepVar,M)';% Get RIP Scores 
                  if ParInfo.RIPNormalize, rip = SNP.normalizeRIP(rip,M,var);end% rescale RIP scores
                  FoldResults{fo}.RIP = rip;
                  [FoldResults{fo}.F,FoldResults{fo}.FP,FoldResults{fo}.G,FoldResults{fo}.GP] = SNP.performRIPStat(FoldResults{fo}.IndVar,FoldResults{fo}.RIP,var,ParInfo.RIPStatPerm);
              end
             % Gathering Fold Results
               FoldF = zeros(1,ParInfo.OuterFold);
               FoldFP = zeros(1,ParInfo.OuterFold);
               FoldG = zeros(1,ParInfo.OuterFold);
               FoldGP = zeros(1,ParInfo.OuterFold);
               BootIter = zeros(1,ParInfo.OuterFold);
               M = zeros(ParInfo.OuterFold,dim);
               for fo=1:ParInfo.OuterFold
                   FoldF(fo) = FoldResults{fo}.F; FoldFP(fo) = FoldResults{fo}.FP;
                   FoldG(fo) = FoldResults{fo}.G;FoldGP(fo) = FoldResults{fo}.GP;
                   M(fo,:) = FoldResults{fo}.M;BootIter(fo) = FoldResults{fo}.BootIter;  
               end
               out.FoldF = FoldF;out.FoldFP = FoldFP;out.FoldFP(out.FoldFP==0) = 0.001/ParInfo.RIPStatPerm;
               out.FP = SNP.pfast(out.FoldFP);out.F = nanmedian(FoldF);
               out.FoldG = FoldG;out.FoldGP = FoldGP;out.FoldGP(out.FoldGP==0) = 0.001/ParInfo.RIPStatPerm;
               out.GP = SNP.pfast(out.FoldGP);out.G = nanmedian(FoldG);
               out.BootIter = BootIter;
              % Extracting final M 
                out.M = median(M);
                if ~ParInfo.Htest, return; end
                if ParInfo.OuterFold>=8, 
                   out.H = SNP.wilcoxonM(M);
                else
                   out.H = ones(1,dim); 
                end
                out.M = out.H.*out.M;
       end
       function [M,H,Mfor] = robustPLSRegression(obj,indvar,depvar,gg,parinfo)
             if nargin<5, parinfo = obj.ParCont; end
           % initialize
             sampling = parinfo.RegSampling;runs = parinfo.RegRuns;balance = parinfo.RegBalance;% only read once
             if runs==1, sampling = false; end% in a single run no point in sampling
             if balance % determine group memberships and balances, within test such that not executed when not necessary
                balanceTH = parinfo.RegBalanceTH;balanceMethod = parinfo.RegBalanceMethod;
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
                            [DepVarNew,IndVarNew] = SNP.regUpSample(depvar(balInfo{i}.Ind,:),indvar(balInfo{i}.Ind,:),balInfo{i});
                            SampleDepVar = [SampleDepVar; DepVarNew];
                            SampleIndVar = [SampleIndVar; IndVarNew];
                        end
                    else % Balancing by Decreasing Majority Groups
                        SampleIndVar = [];SampleDepVar = [];
                        for i=1:1:nrgg
                            [DepVarNew,IndVarNew] = SNP.regDownSample(depvar(balInfo{i}.Ind,:),indvar(balInfo{i}.Ind,:),balInfo{i}.nrKeep);
                            SampleDepVar = [SampleDepVar; DepVarNew];
                            SampleIndVar = [SampleIndVar; IndVarNew];
                        end
                    end
                 else % No Balancing
                     SampleIndVar = indvar;SampleDepVar = depvar;
                 end
                 nrS = size(SampleIndVar,1);
                 if sampling % Bootstrap Sampling of given input
                    SampleInd = randsample(nrS,nrS,true);
                 else
                    SampleInd = 1:nrS;
                 end
                 Mfor(s,:) = SNP.getRegression(SampleIndVar(SampleInd,:),SampleDepVar(SampleInd,:));
             end
           % finalizing output
             if ~parinfo.Htest, H=ones(1,nrD);return; end
             if runs==1, H=ones(1,nrD); return; end
             M = median(Mfor); % extract the median of all regression coefficients
             if nargout<2, return; end
             H = SNP.wilcoxonM(Mfor);% wilcoxon test for median
       end
   end
   methods % Model functions
       function [C,Atypicality,scan] = createMorph(obj,BaseCont,genotype,BF,maxAtypicality)
                if nargin<5, maxAtypicality = +inf; end
                if nargin<4, BF = 0; end
                C = []; Atypicality = +inf;scan = [];
                if isempty(obj.ST4MOD), return; end
                switch genotype
                    case -1 % AA
                        GTInfo = obj.ST4MOD.GTInfo{1};
                    case 0  % AB
                        GTInfo = obj.ST4MOD.GTInfo{2};
                    case 1  % BB
                        GTInfo = obj.ST4MOD.GTInfo{3};
                    otherwise
                        return;
                end
                if GTInfo.Nr==0, return; end
                if obj.ST4RedShape
                   Xref = obj.ST4MOD.RIPRef;
                   Yref = BaseCont.ST4Ref.RedCoeff;
                else
                   Yref = BaseCont.ST4Ref.Coeff; 
                   switch BaseCont.RedShapeSpaceType
                        case 'X'
                           Cov = [BaseCont.ST4Ref.AvgCov, BaseCont.ST4Ref.AvgGB];
                        case 'RIP'
                           Cov = [BaseCont.ST4Ref.CovRIP, BaseCont.ST4Ref.GBRIP];
                    end 
                    Xref = [Cov obj.ST4MOD.RIPRef]; 
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
       function [C,Atypicality,scan] = overLaySNP(obj,BaseCont,genotype,Ref,BF,maxAtyp,adjust,TrimBF)
                if nargin<8, TrimBF = 2; end
                if nargin<7, adjust = 'no'; end
                if nargin<6, maxAtyp = +inf; end
                if nargin<5, BF = 0; end
                if nargin<4, Ref.Coeff = BaseCont.ST4Ref.Coeff; Ref.Cov = [BaseCont.ST4Ref.CovRIP, BaseCont.ST4Ref.GBRIP];end
                C = []; Atypicality = 0;scan = [];
                if isempty(obj.ST4MOD), return; end
                switch genotype
                    case -1 % AA
                        GTInfo = obj.ST4MOD.GTInfo{1};
                    case 0  % AB
                        GTInfo = obj.ST4MOD.GTInfo{2};
                    case 1  % BB
                        GTInfo = obj.ST4MOD.GTInfo{3};
                    otherwise
                        return;
                end
                if GTInfo.Nr==0, return; end
                Xref = SNP.updateRIP(Ref.Coeff,obj.ST4MOD.M);
                Xref = SNP.normalizeRIP(Xref,obj.ST4MOD.M,obj.ST4MOD.Var);
                Xref = [Ref.Cov Xref];Yref = Ref.Coeff;
                Atypref = sqrt(sum(Ref.Coeff.^2));
                maxcounter = 100;counter = 0;
                AtypIncrease = +inf;
                while AtypIncrease>=maxAtyp && counter<=maxcounter
                      counter = counter+1;
                      X = Xref;X(end) = GTInfo.AVG+GTInfo.Sign*BF*GTInfo.STD;
                      deltaX = X-Xref;dY = deltaX*obj.ST4MOD.MMorphs;
                      C = Yref+dY;
                      Atypicality = sqrt(sum(C.^2));
                      AtypIncrease = Atypicality-Atypref;
                      BF = BF-0.01;
                end
                if nargout < 3, return; end
                if obj.ST4RedShape
                   ShapeSpace = BaseCont.RedShapeSpace;
                else
                   ShapeSpace = BaseCont.ShapeSpace;
                end
                switch lower(adjust);
                    case 'rescale'
                        C = C.*ShapeSpace.EigStd';
                    case 'trim'
                        B = TrimBF*ShapeSpace.EigStd';
                        index = find(abs(C)>B);
                        C(index) = B(index).*sign(C(index));
                    otherwise
                end
                scan = getScan(ShapeSpace,C);   
       end
       function [IB,pTpF,pT] = matchRIP(obj,rip,genotype,BF)
                 n = length(rip);IB = nan*zeros(1,n);
                 if nargout>1, pTpF = nan*zeros(1,n); end
                 if nargout>2, pT = nan*zeros(1,n); end
                 switch genotype
                     case -1
                         GTInfo = obj.ST4MOD.GTInfo{1};
                     case 0
                         GTInfo = obj.ST4MOD.GTInfo{2};
                     case 1
                         GTInfo = obj.ST4MOD.GTInfo{3};
                     otherwise
                         return;
                 end
                 if GTInfo.Nr==0, return; end
                 mu = GTInfo.AVG+GTInfo.Sign*BF*GTInfo.STD;sigma = GTInfo.STD;
                 pT = normpdf(rip,mu,sigma);
                 pTmax = normpdf(mu,mu,sigma);
                 pF = zeros(length(GTInfo.Opp),n);
                 for i=1:1:length(GTInfo.Opp) 
                     OppInfo = obj.ST4MOD.GTInfo{GTInfo.Opp(i)};
                     pF(i,:) = normpdf(rip,OppInfo.AVG,OppInfo.STD);
                 end
                 if size(pF,1)==1,
                    IB = pT./(pT+pF);
                    pTpF = pT./pF;
                 else
                    %IB = ((pT./(pT+pF(1,:)))+(pT./(pT+pF(2,:))))/2;
                    %pTpF = ((pT./pF(1,:)).*(pT./pF(2,:)));
                    pF = sum(pF);
                    %pF = pF(1,:).*pF(2,:);
                    IB = pT./(pT+pF);
                    pTpF = (2*pT)./pF;
                    %pTpF = (pT.^2)./(pF(1,:).*pF(2,:));
                 end
                 pT = pT./pTmax;
                 %pTIB = pT.*IB;
                 %if size(pF,1)>1, pF = sum(pF);end
                 %IB = pT./(pT+pF);
                 %pTpF = pT./pF;
                 %pT = pT./pTmax; 
       end
       function [IB,pFpT] = inversematchRIP(obj,rip,genotype,BF)
                 n = length(rip);IB = nan*zeros(1,n);pFpT = nan*zeros(1,n);
                 %if nargout>1, pTIB = nan*zeros(1,n); end
                 %if nargout>2, pT = nan*zeros(1,n); end
                 switch genotype
                     case -1
                         GTInfo = obj.ST4MOD.GTInfo{1};
                     case 0
                         GTInfo = obj.ST4MOD.GTInfo{2};
                     case 1
                         GTInfo = obj.ST4MOD.GTInfo{3};
                     otherwise
                         return;
                 end
                 if GTInfo.Nr==0, return; end
                 mu = GTInfo.AVG+GTInfo.Sign*BF*GTInfo.STD;sigma = GTInfo.STD;
                 pT = normpdf(rip,mu,sigma);
                 %pTmax = normpdf(mu,mu,sigma);
                 pF = zeros(length(GTInfo.Opp),n);
                 for i=1:1:length(GTInfo.Opp) 
                     OppInfo = obj.ST4MOD.GTInfo{GTInfo.Opp(i)};
                     pF(i,:) = normpdf(rip,OppInfo.AVG,OppInfo.STD);
                 end
                 npF = size(pF,1);
                 if npF==1, IB = pF./(pT+pF); return; end
                 IB = zeros(1,npF);
                 pFpT = zeros(1,npF);
                 indpF = (1:npF);
                 for i=1:1:npF
                     IB(i) = pF(i,:)./(pT+pF(i,:)+pF(setdiff(indpF,i),:));
                     pFpT(i) = pF(i,:)./pT;
                 end
                 IB = mean(IB);
                 pFpT = mean(pFpT);
                 %if size(pF,1)>1, pF = sum(pF);end
                 %IB = pF./(pT+pF);
                 %if nargout==1, return; end
                 %pT = pF./pFmax;
                 %pTIB = pF.*IB;
       end
       function [IB,pTIB,pT] = matchFaces(obj,faces,genotype,BF)
                rip = SNP.updateRIP(faces,obj.ST4MOD.M);
                rip = SNP.normalizeRIP(rip,obj.ST4MOD.M,obj.ST4MOD.Var);
                [IB,pTIB,pT] = matchRIP(obj,rip,genotype,BF);
                %if nargout == 1, IB = matchRIP(obj,rip,genotype,BF); return; end
                %if nargout == 2, [IB,pTIB] = matchRIP(obj,rip,genotype,BF); return; end
                %if nargout == 3, [IB,pTIB,pT] = matchRIP(obj,rip,genotype,BF); end
       end
       function [class] = classifyRIP(obj,rip)
                n = length(rip); class = nan*ones(1,n);
                if obj.ST4MOD.Var.nrEl==3, return; end % cannot classify
                rip = rip*obj.ST4MOD.Var.PosClass;
                class = -1*ones(1,n);
                class(rip > obj.ST4MOD.RIPGT) = 1;
                class = class*obj.ST4MOD.Var.PosClass;
       end
       function [out1,out2] = classifyErrorRIP(obj,rip,genotype,BF,T)
                IB = matchRIP(obj,rip,genotype,BF);
                out1 = zeros(size(IB));
                out1(IB<T) = 1;
                out2 = IB;out2(IB>=T) = 1;
       end
       function out = illustrateModel(obj,BF)
           disp([obj.RS '  Valid:' num2str(obj.ST3Valid)]);
           disp('   ');
           disp(obj.ST4Label);
           disp('   ');
           disp(['Atyp: ' num2str(obj.ST4MOD.Atyps) '  Balance:' num2str(obj.ST4MOD.BalAtyp)]);
           disp('   ');
           disp(['G: ' num2str(obj.ST4MOD.RIPG)  '  GT: '  num2str(obj.ST4MOD.RIPGT)]);
           disp('   ');
           disp(['FP: ' num2str(obj.ST4MOD.FP)  '  GP: '  num2str(obj.ST4MOD.GP)]);
           disp('   ');
           disp(num2str([obj.nrAA obj.nrAB obj.nrBB]));
           disp('   ');
           rip = min(obj.ST4MOD.RIP):0.01:max(obj.ST4MOD.RIP);
           M = zeros(1,3);Mrip = zeros(1,3);
           gt = [-1 0 1];
           f1 = figure;hold on;
           plot(rip,0.5*ones(size(rip)),'k-');
           %data = zeros(3,4);
           for i=1:1:3
               if obj.ST4MOD.GTInfo{i}.Nr == 0, disp(['NO ' num2str(gt(i))]); continue; end
               [score] = matchRIP(obj,rip,gt(i),BF);
               plot(f1.CurrentAxes,rip,score);
               [M(1,i),ind] = min(score);Mrip(1,i) = rip(ind);
               %[score] = matchRIP(obj,obj.ST4MOD.GTInfo{i}.MED.mu,gt(i),BF);
               %data(i,1) = obj.ST4MOD.GTInfo{i}.MED.mu;
               %data(i,2) = obj.ST4MOD.GTInfo{i}.MAD.mu;
               %data(i,3) = obj.ST4MOD.GTInfo{i}.MED.mu+obj.ST4MOD.GTInfo{i}.Sign*obj.ST4MOD.GTInfo{i}.MAD.mu;
               %data(i,4) = score;
               %scan = obj.ST4MOD.GTInfo{i}.Morph;
               %color = 0.6*ones(1,3);color(i) = 0.8;
               %scan.SingleColor = color;scan.Material = 'Dull';
               %v = viewer(obj.ST4MOD.GTInfo{i}.Morph);
               %v.SceneLightVisible = true;v.SceneLightLinked = true;
           end
           %obj.ST4Ref.Scan.ColorMode = 'Indexed';v = viewer(obj.ST4Ref.Scan);
           %obj.ST4Ref.Scan.Material = 'Dull';
           %v.SceneLightVisible = true;v.SceneLightLinked = true;
           %disp(num2str([M; Mrip]));
           %disp('   ');
           %disp(num2str(data));
           %disp('   ');
           out.M = M;out.Mrip = Mrip;
           out.Atyp = obj.ST4MOD.Atyps;
           out.G = obj.ST4MOD.RIPG;
       end
       function [freq,occ] = GTTypicality(obj,genotype)
           total = obj.nrAA+obj.nrAB+obj.nrBB;
           switch genotype
               case -1
                   occ = (obj.nrAA/total);
                   %atyp = obj.ST4MOD.Atyps(1);
                   freq = obj.ST4MOD.Freq(1);
               case 0
                   occ = (obj.nrAB/total);
                   %atyp = obj.ST4MOD.Atyps(2);
                   freq = obj.ST4MOD.Freq(2);
               case 1
                   occ = (obj.nrBB/total);
                   %atyp = obj.ST4MOD.Atyps(3);
                   freq = obj.ST4MOD.Freq(3);
               otherwise
                   occ = nan;
                   %atyp = nan;
                   freq = nan;
           end
           if freq==0, freq = nan; end
           %matyp = atyp==obj.ST4MOD.MAtyp;
       end
       function [RIP,Distr] = GT2RIP(obj,genotype)
                 RIP = nan;Distr = cell(1,1);
                 if isempty(obj.ST4MOD), return; end
                 switch genotype
                     case -1
                         Info = obj.ST4MOD.GTInfo{1};
                     case 0
                         Info = obj.ST4MOD.GTInfo{2};
                     case 1
                         Info = obj.ST4MOD.GTInfo{3};
                     otherwise
                         return;
                 end
                 if Info.Nr==0, return; end
                 RIP = Info.AVG;
                 Distr{1}.mu = Info.AVG;
                 Distr{1}.sigma = Info.STD;
       end
       function updateGTDistributions(obj,BaseCont,kappa)
                % Extracting group RIP Statistics  
                 for i=1:1:3
                     obj.ST4MOD.GTInfo{i}.Nr = length(obj.ST4MOD.GTInfo{i}.Ind);
                     if obj.ST4MOD.GTInfo{i}.Nr>0
                        obj.ST4MOD.GTInfo{i}.RIP = obj.ST4MOD.RIP(obj.ST4MOD.GTInfo{i}.Ind);
                        [obj.ST4MOD.GTInfo{i}.AVG,obj.ST4MOD.GTInfo{i}.STD,obj.ST4MOD.GTInfo{i}.MED,obj.ST4MOD.GTInfo{i}.MAD,GTInfo{i}.O] = SNP.extractRobustRIPStatistics(obj.ST4MOD.GTInfo{i}.RIP,kappa);                     
                     end
                 end
                % Creating average morphs
                 obj.ST4MOD.Atyps = nan*zeros(1,3);
                 for i=1:1:3
                     if obj.ST4MOD.GTInfo{i}.Nr==0, continue; end
                     [obj.ST4MOD.GTInfo{i}.C,obj.ST4MOD.GTInfo{i}.Atyp,obj.ST4MOD.GTInfo{i}.Morph] = createMorph(obj,BaseCont,obj.ST4MOD.GTInfo{i}.GT,0);
                     obj.ST4MOD.Atyps(i) = obj.ST4MOD.GTInfo{i}.Atyp;
                 end
                % Determination of less typical group
                 obj.ST4MOD.MAtyp = nanmax(obj.ST4MOD.Atyps);
                 obj.ST4MOD.mAtyp = nanmin(obj.ST4MOD.Atyps);
                 obj.ST4MOD.BalAtyp = obj.ST4MOD.mAtyp/obj.ST4MOD.MAtyp;
       end
       function out = alterAdditiveModel(obj)
                if nargout == 1,obj = clone(obj);out = obj;end
                if ~obj.ST4TestInd==4, return; end
                obj.ST4MOD.GTInfo{2}.Nr = 0;
                obj.ST4MOD.GTInfo{1}.Opp = 3;
                obj.ST4MOD.GTInfo{3}.Opp = 1;
                obj.ST4MOD.Freq(2) = nan;
       end
   end
   methods (Static = true)
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
            H = zeros(1,nrD);
            for j=1:1:nrD
               [~,H(j)] = signrank(M(:,j));% wilcoxon test for median  
            end
       end
       function [out] = normalizeRIP(rip,M,var)
             Mrip = dot(var.MDepVar',M'/norm(M'));
             Prip = dot(var.PDepVar',M'/norm(M'));
             out = ((rip-Mrip)/(Prip-Mrip))*var.Range+var.Mel;
       end
       function [out] = updateRIP(in,M)
            out = dot(in',repmat(M'/norm(M'),1,size(in,1)));
       end
       function [F,FP,Gmean,GmeanP,GT] = performRIPStat(g,rip,var,t)
                if var.nrEl<3, doroc = true; else doroc = false; Gmean = 0; GmeanP = 1; GT = nan; end
                index = find(~isnan(g));g = g(index);
                XG = cell(size(g));
                for l=1:var.nrEl % building grouping variable
                    XG((g==var.El(l))) = {num2str(l)};
                end
                rip = rip(index);n = length(index);
                F = SNP.ripAnova(XG,rip);FCount = zeros(1,t);
                if doroc, [Gmean,GT] = SNP.ripROC(g,rip,var.PosClass); GCount = zeros(1,t); end
                if t == 0, FP = 0; GmeanP = 0; return;end
                parfor i=1:t
                    ind = randperm(n);
                    Ffor = SNP.ripAnova(XG,rip(ind));
                    FCount(i) = Ffor>=F;
                    if doroc, Gmeanfor = SNP.ripROC(g,rip(ind),var.PosClass); GCount(i) = Gmeanfor>=Gmean; end
                end
                FP = sum(FCount)/t;
                if doroc, GmeanP = sum(GCount)/t; end
                if ~doroc, Gmean = F;GmeanP = FP; end
       end
       function [F,T] = ripAnova(XG,rip)
                [~,T] = anova1(rip,XG,'off');
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
       function outP = interpCodeP1(code)
                codeS = num2str(code);
                switch sum(code)
                    case 0 % 0 0 0 0
                        outP = [5 6 7];
                    case 1
                        outP = find(code);
                    case 2
                        switch codeS
                            case num2str([1 1 0 0])
                                outP = [1 2 4 5];
                            case num2str([1 0 1 0])
                                outP = [1 3 7];
                            case num2str([1 0 0 1])
                                outP = [1 4 5 7];
                            case num2str([0 1 1 0])
                                outP = [2 3 6];
                            case num2str([0 1 0 1])
                                outP = [2 4 5 6];
                            case num2str([0 0 1 1])
                                outP = [3 4 6 7];
                        end
                    case 3
                        switch codeS
                            case num2str([1 1 1 0])
                                outP = [8 1 2 3];
                            case num2str([1 1 0 1])
                                outP = [1 2 4];
                            case num2str([1 0 1 1])
                                outP = [1 3 4 5 6];
                            case num2str([0 1 1 1])
                                outP = [2 3 4 5 7];
                        end
                    case 4 % 1 1 1 1
                        outP = [1 2 3 4 8];
                end
       end
       function [outIP,outP] = interpCodeP2(code1,code0)
                % outP are the possible tests, outIP are the impossible tests
                code1S = num2str(code1);code0S = num2str(code0);
                outP = (1:10);
                switch sum(code1) % All possible differences
                       case 0 % 0 0 0 % no suggestion
                           outIP1 = [];
                       case 1
                           switch code1S
                               case num2str([1 0 0]) % Contrast AA with BB, additive outP = [1 2 4 5 8]; % possible tests
                                    outIP1 = [3 9 10];
                               case num2str([0 1 0]) % Contrast AA with AB
                                    outIP1 = [1 10];
                               case num2str([0 0 1]) % Contrast BB with AB
                                    outIP1 = [2 9];
                           end
                       case 2
                           switch code1S
                               case num2str([0 1 1]) % Over Dom or Add
                                   outIP1 = [1 2];
                               case num2str([1 0 1]) % Rec or add
                                   outIP1 = [2 3 9];
                               case num2str([1 1 0]) % Dom or add
                                   outIP1 = [1 3 10];
                           end
                       case 3 % 1 1 1 Non add or add % strong suggestion
                             outIP1 = [1 2 3 9 10];
                end
                switch sum(code0) % All possible non-differences
                    case 0 % 0 0 0 % no suggestion
                           outIP0 = [];
                       case 1
                           switch code0S
                               case num2str([1 0 0]) % DO NOT contrast AA with BB
                                    outIP0 = [1 2 4 5 8];
                               case num2str([0 1 0]) % DO NOT Contrast AA with AB
                                    outIP0 = [2 3 4 6 8 9]; 
                               case num2str([0 0 1]) % DO NOT Contrast BB with AB
                                    outIP0 = [1 3 4 7 8 10];
                           end
                       case 2
                           switch code0S
                               case num2str([0 1 1]) % DO NOT contrast AA with AB and BB with AB
                                   outIP0 = [1 2 3 4 6 7 8 9 10];
                               case num2str([1 0 1]) % DO NOT contrast AA with BB and BB with AB
                                   outIP0 = [1 2 3 4 5 7 8 10];
                               case num2str([1 1 0]) % DO NOT contrast AA with AB and BB with AA
                                   outIP0 = [1 2 3 4 5 6 8 9];
                           end
                       case 3 % DO NOT CONTRAST ANYTHING
                             outIP0 = [1 2 3 4 5 6 7 8 9 10];
                end
                outIP = union(outIP1,outIP0);
                outP = setdiff(outP,outIP);
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
       function [AVG,STD,MED,MAD,O] = extractRIPStatistics(rip,runs,nrs)
                % Initialize
                  nrO = length(rip);
                  tmpAVG = nan*zeros(1,runs);tmpSTD = nan*zeros(1,runs);tmpMED = nan*zeros(1,runs);tmpMAD = nan*zeros(1,runs);
                % bootstrap sampling  
                  for i=1:1:runs 
                      s = randsample(nrO,nrs,true);
                      tmpAVG(i) = mean(rip(s));
                      tmpSTD(i) = std(rip(s));
                      tmpMED(i) = median(rip(s));
                      tmpMAD(i) = mad(rip(s));
                  end
                  AVG.mu = mean(tmpAVG);AVG.sigma = std(tmpAVG);
                  STD.mu = mean(tmpSTD);STD.sigma = std(tmpSTD);
                  MED.mu = mean(tmpMED);MED.sigma = std(tmpMED);
                  MAD.mu = mean(tmpMAD);MAD.sigma = std(tmpMAD);
                  if nargout < 5, return; end
                  O.avg = mean(rip);O.std = std(rip);O.med = median(rip);O.mad = mad(rip);
       end
   end    
end


% function out = get.ST4ShapeSpace(obj)
%            if isempty(obj.ParCont), out = []; return; end
%            if obj.ST4RedShape, out = obj.ParCont.RedShapeSpace; return; end 
%            out = obj.ParCont.ShapeSpace;
% end
% function out = get.ST4RefScan(obj)
%                 out = obj.ST4ShapeSpace.Average;
% end

%        function out = get.n(obj)
%            out = length(obj.GT);
%        end

%        function out = get.AAavgScan(obj)
%            if isempty(obj.AAind), out = []; return; end
%            if isempty(obj.RedShapeSpace), out = []; return; end
%            out = getScan(obj.RedShapeSpace,mean(obj.RedDepVar(obj.AAind,:)));
%        end
%        function out = get.BBavgScan(obj)
%            if isempty(obj.BBind), out = []; return; end
%            if isempty(obj.RedShapeSpace), out = []; return; end
%            out = getScan(obj.RedShapeSpace,mean(obj.RedDepVar(obj.BBind,:)));
%        end
%        function out = get.ABavgScan(obj)
%            if isempty(obj.ABind), out = []; return; end
%            if isempty(obj.RedShapeSpace), out = []; return; end
%            out = getScan(obj.RedShapeSpace,mean(obj.RedDepVar(obj.ABind,:)));
%        end

% function out = get.Cov(obj)
%            if isempty(obj.ParCont), out = []; return; end
%            out = obj.ParCont.Cov;
%        end
%        function out = get.CovNames(obj)
%            if isempty(obj.ParCont), out = []; return; end
%            out = obj.ParCont.CovNames;
%        end



%        function out = get.ShapeSpace(obj)
%            if isempty(obj.ParCont), out = []; return; end
%            out = obj.ParCont.ShapeSpace;
%        end
%        function out = get.RedShapeSpace(obj)
%            if isempty(obj.ParCont), out = []; return; end
%            out = obj.ParCont.RedShapeSpace;
%        end
%        function out = get.DepVar(obj)
%            if isempty(obj.ParCont), out = []; return; end
%            out = obj.ParCont.DepVar;
%        end
%        function out = get.RedDepVar(obj)
%            if isempty(obj.ParCont), out = []; return; end
%            out = obj.ParCont.RedDepVar;
%        end






%            str = input('Do you want to DISCARD (0), KEEP first (1), second (2) or third(3) or ALTER (4) sugestion? [0,1,2,3,4]: ','s');
%            disp('  ');
%            switch str
%                case '0'
%                  obj.ST3TestInd = [0 0];
%                case '1'
%                  obj.ST3TestInd = obj.ST1TestInd;
%                case '2'
%                  obj.ST3TestInd = obj.ST2FTestInd;
%                case '3'
%                  obj.ST3TestInd = obj.ST2GTestInd;
%                case '4'
%                  disp('0      : NO EFFECT');disp('1      : RECESSIVE');disp('2      : DOMINANT');
%                  disp('3      : OVER DOMINANT');disp('4      : ADDITIVE');disp('5      : AA BB');
%                  disp('6      : AA AB');disp('7      : AB BB');disp('8      : NON ADDITVE');
%                  disp('9      : ADDITIVE AA AB');disp('10      : ADDITIVE AB BB');
%                  str2 = input('Enter DESIRED inheritance pattern [0-10]','s');
%                  obj.ST3TestInd = [str2double(str2) 0];  
%            end           




%        function tmpcodeInterpretation(obj,check)
%            % stage 2 is a control of the suggested testindex
%            obj.ST3TestInd = obj.ST1TestInd;
%            if check          
%                disp(obj.RS)
%                disp('  ');
%                disp('SUGGESTION STAGE 1: ');
%                disp([obj.Labels{1} '     ' obj.Labels{2} '     ' obj.Labels{3} '     ' obj.Labels{4} '       ' obj.Labels{5} '       ' obj.Labels{6}]); 
%                tmp = zeros(2,6);tmp(1,:) = obj.ST1P;tmp(2,:) = obj.ST1Code;
%                disp(num2str(tmp));
%                disp('  ');
%                disp(['TEST ID: ' num2str(obj.ST1TestInd(1)) ' ==> ' obj.ST1Label]);
%                disp('  ');
%                disp('SUGGESTION STAGE 2: ');
%                disp([obj.Labels{1} '     ' obj.Labels{2} '     ' obj.Labels{3} '     ' obj.Labels{4} '       ' obj.Labels{5} '       ' obj.Labels{6}]); 
%                tmp = zeros(2,6);tmp(1,:) = obj.ST2P;tmp(2,:) = obj.ST2Code;
%                disp(num2str(tmp));
%                disp('  ');
%                disp(['TEST ID: ' num2str(obj.ST2TestInd(1)) ' ==> ' obj.ST2Label]);
%                disp('  ');
%                disp('nrAA nrAB nrBB');
%                tmp = zeros(1,3);tmp(1) = obj.nrAA;tmp(2) = obj.nrAB; tmp(3) = obj.nrBB;
%                disp(num2str(tmp));
%                disp('  ');
%                str = input('Do you want to ALTER (0) or KEEP first (1) or second (2) sugestion? [0,1,2]: ','s');
%                disp('  ');
%                switch str
%                    case '0'
%                        disp('0      : NO EFFECT');disp('1      : RECESSIVE');disp('2      : DOMINANT');
%                        disp('3      : OVER DOMINANT');disp('4      : ADDITIVE');disp('5      : NON ADDITIVE');
%                        disp('6      : AA BB AB');disp('7      : AB AA BB');
%                        str2 = input('Enter DESIRED inheritance pattern [0-7]','s');
%                        obj.ST3TestInd = [str2double(str2) 0];
%                    case '1'
%                        obj.ST3TestInd = obj.ST1TestInd;
%                    case '2'
%                        obj.ST3TestInd = obj.ST2TestInd;
%                end
%            end
%        end




%            switch sum(codeP1)
%                case 0 % 0 0 0 0 NO SUGGESTION
%                    switch sum(codeP2)
%                        case 0 % 0 0 0 % no suggestion
%                            out = [0 1];%'No Effect';
%                        case 1
%                             switch codeP2S
%                                 case num2str([1 0 0]) % Contrast AA with BB
%                                     out = [5 1];
%                                 case num2str([0 1 0]) % Contrast AA with AB
%                                     out = [6 1];
%                                 case num2str([0 0 1]) % Contrast BB with AB
%                                     out = [7 1];
%                             end
%                        case 2
%                             switch codeP2S
%                                 case num2str([0 1 1]) % Over Dom or Add
%                                      out = [3 0];% Weak Over dominant
%                                 case num2str([1 0 1]) % Rec or add
%                                      out = [1 0];% Weak rec
%                                 case num2str([1 1 0]) % Dom or add
%                                      out = [2 0];%'Weak Dom';
%                             end
%                        case 3 % 1 1 1 % non add or add % strong suggestion
%                            if P(4)<=P(3)
%                               out = [4 0];% Weak Add
%                            else
%                               out = [8 0];% Weak Non Add
%                            end
%                    end
%                case 1
%                    switch codeP1S
%                        case num2str([1 0 0 0]) % Rec
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [1 0];%'Weak Rec';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [1 0];%'Weak Rec';
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [0 0];% Contradiction
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [1 0];%'Weak Rec';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 1];%'Rec';
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add
%                                    out = [1 0];%'Contradiction';
%                            end
%                        case num2str([0 1 0 0]) % Dom
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [2 0];%'Weak Dom';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [2 0];%'Weak Dom';
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [2 0];%'Weak Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 1];%'Dom';
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [0 0];%'Contradiction';
%                            end
%                        case num2str([0 0 1 0]) % Over Dom
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [3 0];%'Weak Over Dom';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [0 0];% Contradiction
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 1];%'Over Dom';
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [0 0];%'Contradiction';
%                            end
%                        case num2str([0 0 1 0]) % Additief
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [4 0];%'Weak Add';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [4 0];% 'Weak Add'
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [4 0];% 'Weak Add';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [4 0];% 'Weak Add';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [4 0];% 'Weak Add';
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [4 0];% 'Weak Add'
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [4 0];% 'Weak Add'
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 1];% 'Additive';
%                            end
%                    end
%                case 2
%                    switch codeP1S
%                        case num2str([1 1 0 0]) % recessive/dominant/additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(1)<P(2)
%                                       out = [1 0]; % Weak rec
%                                    elseif P(1)>P(2)
%                                       out = [2 0]; % Weak dom
%                                    else
%                                       out = [5 0];% Contrast AA with BB
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [5 1];% Contrast AA with BB
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [2 0];% 'Weak Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [1 0];% 'Weak Rec';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 1];% Rec
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 1];% Dom
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 0];% 'Weak Add';
%                            end
%                        case num2str([1 0 1 0]) % recessive/overdominant/non additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(1)<P(3)
%                                       out = [1 0]; % Weak rec
%                                    elseif P(1)>P(3)
%                                       out = [3 0]; % Over dominant Weak dom
%                                    else
%                                       out = [7 0];% Contrast AB with BB
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [1 0];
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [3 0];
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             if P(1)<P(3)
%                                               out = [1 0]; % Weak rec
%                                             elseif P(1)>P(3)
%                                               out = [3 0]; % Over dominant Weak dom
%                                             else
%                                               out = [7 0];% Contrast AB with BB
%                                             end
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 0];% Rec
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [8 0];% 'Weak non Add';
%                            end
%                        case num2str([1 0 0 1]) % recessive/additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(1)<P(4)
%                                       out = [1 0]; % Weak rec
%                                    else
%                                       out = [4 0];% Add
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             if P(1)<P(4)
%                                                out = [1 0]; % Weak rec
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [4 0];
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             if P(1)<P(4)
%                                                out = [1 0]; % Weak rec
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 0];% Rec
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [4 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 0];% 'Weak Add';
%                            end
%                        case num2str([0 1 1 0]) % Dominant/ Over Dominant 
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(2)<P(3)
%                                       out = [2 0]; % Weak rec
%                                    else
%                                       out = [3 0];% Over Dom
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [2 0];
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             if P(2)<P(3)
%                                                out = [2 0]; % Weak rec
%                                             else
%                                                out = [3 0];% Over Dom
%                                             end
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [3 0];
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 0];% Over Dom
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 0];% Dom
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [8 0];% 'Weak Non Add';
%                            end
%                        case num2str([0 1 0 1]) % Dominant/additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(2)<P(4)
%                                       out = [2 0]; % Weak Dom
%                                    else
%                                       out = [4 0];% Add
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             if P(2)<P(4)
%                                                out = [2 0]; % Weak Dom
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             if P(2)<P(4)
%                                                out = [2 0]; % Weak Dom
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [4 0];
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [4 0];% Add
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [4 0];% Add
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 0];% Dom
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 0];% 'Weak Add';
%                            end
%                        case num2str([0 0 1 1]) % Over dominance/Additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(3)<P(4)
%                                       out = [3 0]; % Over Dom
%                                    else
%                                       out = [4 0];% Add
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [4 0];
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             if P(3)<P(4)
%                                                out = [3 0]; % Over Dom
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             if P(3)<P(4)
%                                                out = [3 0]; % Over Dom
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 0];% Over dom
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [4 0];% Add
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [4 0];% Add
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 0];% 'Weak Add'; could also be non-additive
%                            end
%                    end
%                case 3
%                    switch codeP1S
%                        case num2str([1 1 1 0])
%                            switch codeP2S
%                                case num2str([1 1 1])
%                                    out = [8 1];% Strong Non Additive
%                                otherwise
%                                    out = [8 0];% Weak Non Additive
%                            end
%                        case num2str([1 1 0 1])
%                            switch codeP2S
%                                case num2str([1 1 1])
%                                    out = [4 1];% Strong Additive
%                                otherwise
%                                    out = [4 0];% Weak Additive
%                            end
%                        case num2str([1 0 1 1])
%                            ind = min(P(1:4))
%                            switch ind
%                                case 1
%                                    switch codeP2S
%                                        case num2str([1 0 1])
%                                            out = [1 1];% Strong Rec
%                                        otherwise
%                                            out = [1 0];% Weak Rec
%                                    end
%                                case 3
%                                    switch codeP2S
%                                        case num2str([0 1 1])
%                                            out = [1 1];% Strong Rec
%                                        otherwise
%                                            out = [1 0];% Weak Rec
%                                    end
%                                case 4
%                            end
%                            
%                        case num2str([0 1 1 1])
%                    end
%                case 4 % 1 1 1 1
%                    switch codeP2S
%                           case num2str([1 1 1])
%                                out = [8 1];% Strong Non Additive
%                           otherwise
%                                out = [8 0];% Weak Non Additive
%                    end   
%            end




%        function out = retrieveTestInd(obj,stage,stat)
%            switch stage
%                case 1
%                    code = obj.ST1Code;P = obj.ST1P;
%                case 2
%                    switch stat
%                        case 'F' % Using F statistics only
%                             code = obj.ST2FCombCode;P = obj.STFPFish;
%                        case 'G' % Using G statistics only
%                             code = obj.ST2GCombCode;P = obj.STGPFish;
%                        case 'FG' % Using F and G statistcs combined
%                             code = obj.ST2Code;P = obj.ST2P;
%                    end
%            end
%            if sum(code)==0, out = [0 1]; return; end % 0 0 0 0 0 0
%            codeP1 = code(1:4);codeP1S = num2str(codeP1);% interpretation part
%            codeP2 = code(5:7);codeP2S = num2str(codeP2);% suggestive or confirmation part
%            out = [0 0];
%            switch sum(codeP1)
%                case 0 % 0 0 0 0 NO SUGGESTION
%                    switch sum(codeP2)
%                        case 0 % 0 0 0 % no suggestion
%                            out = [0 1];%'No Effect';
%                        case 1
%                             switch codeP2S
%                                 case num2str([1 0 0]) % Contrast AA with BB
%                                     out = [5 1];
%                                 case num2str([0 1 0]) % Contrast AA with AB
%                                     out = [6 1];
%                                 case num2str([0 0 1]) % Contrast BB with AB
%                                     out = [7 1];
%                             end
%                        case 2
%                             switch codeP2S
%                                 case num2str([0 1 1]) % Over Dom or Add
%                                      out = [3 0];% Weak Over dominant
%                                 case num2str([1 0 1]) % Rec or add
%                                      out = [1 0];% Weak rec
%                                 case num2str([1 1 0]) % Dom or add
%                                      out = [2 0];%'Weak Dom';
%                             end
%                        case 3 % 1 1 1 % non add or add % strong suggestion
%                            if P(4)<=P(3)
%                               out = [4 0];% Weak Add
%                            else
%                               out = [8 0];% Weak Non Add
%                            end
%                    end
%                case 1
%                    switch codeP1S
%                        case num2str([1 0 0 0]) % Rec
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [1 0];%'Weak Rec';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [1 0];%'Weak Rec';
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [0 0];% Contradiction
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [1 0];%'Weak Rec';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 1];%'Rec';
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add
%                                    out = [1 0];%'Contradiction';
%                            end
%                        case num2str([0 1 0 0]) % Dom
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [2 0];%'Weak Dom';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [2 0];%'Weak Dom';
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [2 0];%'Weak Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 1];%'Dom';
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [0 0];%'Contradiction';
%                            end
%                        case num2str([0 0 1 0]) % Over Dom
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [3 0];%'Weak Over Dom';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [0 0];% Contradiction
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 1];%'Over Dom';
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [0 0];%'Contradiction';
%                            end
%                        case num2str([0 0 1 0]) % Additief
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [4 0];%'Weak Add';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [4 0];% 'Weak Add'
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [4 0];% 'Weak Add';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [4 0];% 'Weak Add';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [4 0];% 'Weak Add';
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [4 0];% 'Weak Add'
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [4 0];% 'Weak Add'
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 1];% 'Additive';
%                            end
%                    end
%                case 2
%                    switch codeP1S
%                        case num2str([1 1 0 0]) % recessive/dominant/additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(1)<P(2)
%                                       out = [1 0]; % Weak rec
%                                    elseif P(1)>P(2)
%                                       out = [2 0]; % Weak dom
%                                    else
%                                       out = [5 0];% Contrast AA with BB
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [5 1];% Contrast AA with BB
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [2 0];% 'Weak Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [1 0];% 'Weak Rec';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 1];% Rec
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 1];% Dom
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 0];% 'Weak Add';
%                            end
%                        case num2str([1 0 1 0]) % recessive/overdominant/non additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(1)<P(3)
%                                       out = [1 0]; % Weak rec
%                                    elseif P(1)>P(3)
%                                       out = [3 0]; % Over dominant Weak dom
%                                    else
%                                       out = [7 0];% Contrast AB with BB
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [1 0];
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [3 0];
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             if P(1)<P(3)
%                                               out = [1 0]; % Weak rec
%                                             elseif P(1)>P(3)
%                                               out = [3 0]; % Over dominant Weak dom
%                                             else
%                                               out = [7 0];% Contrast AB with BB
%                                             end
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 0];% Rec
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [8 0];% 'Weak non Add';
%                            end
%                        case num2str([1 0 0 1]) % recessive/additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(1)<P(4)
%                                       out = [1 0]; % Weak rec
%                                    else
%                                       out = [4 0];% Add
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             if P(1)<P(4)
%                                                out = [1 0]; % Weak rec
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [4 0];
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             if P(1)<P(4)
%                                                out = [1 0]; % Weak rec
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 0];% Rec
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [4 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 0];% 'Weak Add';
%                            end
%                        case num2str([0 1 1 0]) % Dominant/ Over Dominant 
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(2)<P(3)
%                                       out = [2 0]; % Weak rec
%                                    else
%                                       out = [3 0];% Over Dom
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [2 0];
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             if P(2)<P(3)
%                                                out = [2 0]; % Weak rec
%                                             else
%                                                out = [3 0];% Over Dom
%                                             end
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [3 0];
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 0];% Over Dom
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 0];% Dom
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [8 0];% 'Weak Non Add';
%                            end
%                        case num2str([0 1 0 1]) % Dominant/additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(2)<P(4)
%                                       out = [2 0]; % Weak Dom
%                                    else
%                                       out = [4 0];% Add
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             if P(2)<P(4)
%                                                out = [2 0]; % Weak Dom
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             if P(2)<P(4)
%                                                out = [2 0]; % Weak Dom
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [4 0];
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [4 0];% Add
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [4 0];% Add
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 0];% Dom
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 0];% 'Weak Add';
%                            end
%                        case num2str([0 0 1 1]) % Over dominance/Additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    if P(3)<P(4)
%                                       out = [3 0]; % Over Dom
%                                    else
%                                       out = [4 0];% Add
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [4 0];
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             if P(3)<P(4)
%                                                out = [3 0]; % Over Dom
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             if P(3)<P(4)
%                                                out = [3 0]; % Over Dom
%                                             else
%                                                out = [4 0];% Add
%                                             end
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 0];% Over dom
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [4 0];% Add
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [4 0];% Add
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 0];% 'Weak Add'; could also be non-additive
%                            end
%                    end
%                case 3
%                    switch codeP1S
%                        case num2str([1 1 1 0])
%                            switch codeP2S
%                                case num2str([1 1 1])
%                                    out = [8 1];% Strong Non Additive
%                                otherwise
%                                    out = [8 0];% Weak Non Additive
%                            end
%                        case num2str([1 1 0 1])
%                            switch codeP2S
%                                case num2str([1 1 1])
%                                    out = [4 1];% Strong Additive
%                                otherwise
%                                    out = [4 0];% Weak Additive
%                            end
%                        case num2str([1 0 1 1])
%                            ind = min(P(1:4))
%                            switch ind
%                                case 1
%                                    switch codeP2S
%                                        case num2str([1 0 1])
%                                            out = [1 1];% Strong Rec
%                                        otherwise
%                                            out = [1 0];% Weak Rec
%                                    end
%                                case 3
%                                    switch codeP2S
%                                        case num2str([0 1 1])
%                                            out = [1 1];% Strong Rec
%                                        otherwise
%                                            out = [1 0];% Weak Rec
%                                    end
%                                case 4
%                            end
%                            
%                        case num2str([0 1 1 1])
%                    end
%                case 4 % 1 1 1 1
%                    switch codeP2S
%                           case num2str([1 1 1])
%                                out = [8 1];% Strong Non Additive
%                           otherwise
%                                out = [8 0];% Weak Non Additive
%                    end   
%            end
%        end





%        function out = retrieveTestIndST1(obj,stage)
%            switch stage
%                case 1
%                    code = obj.ST1Code;P = obj.ST1P;
%                case 2
%                    code = obj.ST2Code;P = obj.ST2P;
%            end
%            if sum(code)==0, out = [0 1]; return; end % 0 0 0 0 0 0
%            codeP1 = code(1:3);codeP1S = num2str(codeP1);% interpretation part
%            codeP2 = code(4:6);codeP2S = num2str(codeP2);% suggestive or confirmation part
%            out = [0 0];
%            switch sum(codeP1)
%                case 0 % 0 0 0 weak Additive
%                    switch sum(codeP2)
%                        case 0 % 0 0 0 % no suggestion
%                            out = [0 0];%'No Effect';
%                        case 1
%                             switch codeP2S
%                                 case num2str([1 0 0]) % Contrast AA with BB, additive
%                                     out = [4 0];%'Weak Add';
%                                 case num2str([0 1 0]) % Contrast AA with AB
%                                     out = [6 0];%'AA BB AB';
%                                 case num2str([0 0 1]) % Contrast BB with AB
%                                     out = [7 0];%'AB AA BB';
%                             end
%                        case 2
%                             switch codeP2S
%                                 case num2str([0 1 1]) % Over Dom or Add
%                                     out = [4 0];%'Weak Add';
%                                 case num2str([1 0 1]) % Rec or add
%                                     out = [4 0];%'Weak Add';
%                                 case num2str([1 1 0]) % Dom or add
%                                     out = [4 0];%'Weak Add';
%                             end
%                        case 3 % 1 1 1 % non add or add % strong suggestion
%                            out = [4 1];% Add
%                    end
%                case 1
%                    switch codeP1S
%                        case num2str([1 0 0]) % Rec
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [1 0];%'Weak Rec';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [1 0];%'Weak Rec';
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [0 0];% Contradiction
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [1 0];%'Weak Rec';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [1 1];%'Rec';
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 1];%'Add';
%                            end
%                        case num2str([0 1 0]) % Dom
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [2 0];%'Weak Dom';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [2 0];%'Weak Dom';
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [2 0];%'Weak Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [2 1];%'Dom';
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 1];%'Add';
%                            end
%                        case num2str([0 0 1]) % Over Dom
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [3 0];%'Weak Over Dom';
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [0 0];% Contradiction
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add
%                                             out = [3 1];%'Over Dom';
%                                         case num2str([1 0 1]) % Rec or add
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add
%                                             out = [0 0];% Contradiction
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [5 0];%'Non Add';
%                            end
%                    end
%                case 2
%                    switch codeP1S
%                        case num2str([0 1 1])% Non add or Dom or Over Dom % contrast AA with AB
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    %out = 'AA BB AB';% it does not matter where we put BB, let BRIM figure it out
%                                    % further testing is possible
%                                    if P(2)<=P(3)
%                                       out = [2 0];%'Weak Dom';% preference by occurence in Nature
%                                    else
%                                       if P(4)>=P(6)
%                                          out = [3 0];% 'Weak Over Dom';
%                                       else
%                                          out = [2 0];%'Weak Dom';% preference by occurence in Nature
%                                       end
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Contrast AA with BB, additive
%                                             out = [2 0];%'Weak Dom';% no restriction on BB with AB
%                                         case num2str([0 1 0]) % Contrast AA with AB, no new intell
%                                             % further testing is possible
%                                             if P(2)<=P(3)
%                                                 out = [2 0];%'Weak Dom';% preference by occurence in Nature
%                                             else
%                                                 if P(4)>=P(6)
%                                                    out = [3 0];%'Weak Over Dom';
%                                                 else
%                                                    out = [2 0];%'Weak Dom';% preference by occurence in Nature
%                                                 end
%                                             end
%                                         case num2str([0 0 1]) % Contrast BB with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add Contrast AA with AB, contrast BB with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                         case num2str([1 0 1]) % Rec or add, contrast AA with BB, contrast BB with AB
%                                             out = [0 0];% Contradiction
%                                         case num2str([1 1 0]) % Dom or add. contrast AA with BB, contrast AA with AB
%                                             out = [2 0];%'Weak Dom';
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [5 0];%'Non Add';
%                            end
%                        case num2str([1 0 1]) % Non add or Rec or Over Dom, contrast BB with AB
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion, no new intel
%                                    %out = 'AB AA BB';% it does not matter where we put AA
%                                    % further testing is possible
%                                    if P(1)<=P(3)
%                                       out = [1 0];%'Weak Rec';% preference by occurence in Nature
%                                    else
%                                       if P(4)>=P(5)
%                                          out = [3 0];%'Weak Over Dom';
%                                       else
%                                          out = [1 0];%'Weak Rec';% preference by occurence in Nature
%                                       end
%                                    end
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Additive, contrast AA with BB
%                                             out = [1 0];%'Weak Rec';
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [3 0];%'Weak Over Dom';
%                                         case num2str([0 0 1]) % Contrast BB with AB, no new intel
%                                             % further testing is possible
%                                             if P(1)<=P(3)
%                                                 out = [1 0];%'Weak Rec';% preference by occurence in Nature
%                                             else
%                                                 if P(4)>=P(5)
%                                                    out = [3 0];%'Weak Over Dom';
%                                                 else
%                                                    out = [1 0];%'Weak Rec';% preference by occurence in Nature
%                                                 end
%                                             end
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add, contrast AA with AB, contrast BB with AB 
%                                             out = [3 0];%'Weak Over Dom';
%                                         case num2str([1 0 1]) % Rec or add, contrast AA with BB, contrast BB with AB
%                                             out = [1 0];%'Weak Rec';
%                                         case num2str([1 1 0]) % Dom or add, contrats AA with BB, contrast AA with AB
%                                             out = [5 0];%'Weak Non Add';
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [5 0];%'Non Add';
%                            end
%                        case num2str([1 1 0]) % additive
%                            switch sum(codeP2)
%                                case 0 % 0 0 0 % no suggestion
%                                    out = [4 0];%'Weak Add';% it does not matter where we put AB
%                                case 1
%                                     switch codeP2S
%                                         case num2str([1 0 0]) % Additive, contrast AA with BB
%                                             out = [4 0];%'Weak Add';
%                                         case num2str([0 1 0]) % Contrast AA with AB
%                                             out = [4 0];%'Weak Add';
%                                         case num2str([0 0 1]) % Contrast BB with AB, no new intel
%                                             out = [4 0];%'Weak Add';
%                                     end
%                                case 2
%                                     switch codeP2S
%                                         case num2str([0 1 1]) % Over Dom or Add, contrast AA with AB, contrast BB with AB 
%                                             out = [4 0];%'Weak Add';
%                                         case num2str([1 0 1]) % Rec or add, contrast AA with BB, contrast BB with AB
%                                             out = [4 0];%'Weak Add';
%                                         case num2str([1 1 0]) % Dom or add, contrats AA with BB, contrast AA with AB
%                                             out = [4 0];%'Weak Add';
%                                     end
%                                case 3 % 1 1 1 Non add or add % strong suggestion
%                                    out = [4 1];%'Add';
%                            end
%                            
%                    end
%                case 3% 1 1 1 Non Add
%                    switch sum(codeP2)
%                        case 0 % 0 0 0 % no suggestion
%                             out = [5 0];%'Weak Non Add';% it does not matter where we put AB
%                        case 1
%                             switch codeP2S
%                                    case num2str([1 0 0]) % Additive, contrast AA with BB
%                                         out = [5 0];%'Weak Non Add';
%                                    case num2str([0 1 0]) % Contrast AA with AB
%                                         out = [5 0];%'Weak Non Add';
%                                    case num2str([0 0 1]) % Contrast BB with AB, no new intel
%                                         out = [5 0];%'Weak Non Add';
%                             end
%                        case 2
%                             switch codeP2S
%                                    case num2str([0 1 1]) % Over Dom or Add, contrast AA with AB, contrast BB with AB 
%                                         out = [3 0];%'Weak Over Dom';
%                                    case num2str([1 0 1]) % Rec or add, contrast AA with BB, contrast BB with AB
%                                         out = [1 0];%'Weak Rec';
%                                    case num2str([1 1 0]) % Dom or add, contrats AA with BB, contrast AA with AB
%                                         out = [2 0];%'Weak Dom';
%                             end
%                        case 3 % 1 1 1 Non add or add % strong suggestion
%                             out = [5 1];%'Non Add';
%                    end
%            end
%        end 
