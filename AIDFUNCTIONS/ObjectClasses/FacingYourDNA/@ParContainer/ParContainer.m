classdef ParContainer < superClass
   % A container class containing all the parameter values
   properties
      %ShapeSpace = [];
      %RedShapeSpace = [];
      %Cov = [];
      %CovNames = {};
   end
   properties (Dependent = true, Hidden = true)
      %DepVar;
      %RedDepVar;
      %AvgCov;
      %n;
   end
   % Stage 0 QUALITY CHECKING
   properties
      ST0MinSample = 10;
   end
   % Stage 1 LOCATION BASE TESTS 
   properties
      ST1PT = 0.05;
      ST1MinHits = 1;
   end
   % Stage 2 SNPBRIM BASE TESTS
   properties
      ST2PFFoldT = 0.12;
      ST2PGFoldT = 0.12;
      ST2PFT = 0.05;
      ST2PGT = 0.05;
      ST2NrFoldT = 1;
      ST2MinHits = 1;
   end
   % Stage 3 SNPBRIM COVARIATES TEST
   properties
      ST3PFFoldT = 0.12;
      ST3PGFoldT = 0.12;
      ST3PFT = 0.12;
      ST3PGT = 0.12;
      ST3NrFoldT = 1;
   end
   properties % SNPBRIM
      RegRuns = 50; % Number of Regression runs (SNPBRIM)
      RegSampling = true;% Bootstrap/None (SNPBRIM)
      RegBalance = true; % balance the data yes or no (SNPBRIM)
      RegBalanceTH = 0.3; % desired balance factor (SNPBRIM)
      RegBalanceMethod = 'UpSample'; % UpSample/DownSample/ADASYN (SNPBRIM)
      OuterFold = 4; % number of Outer Folds (SNPBRIM)
      InnerFold = 10; % number of Inner Folds (SNPBRIM)
      UseRedShape = true; % Use reduced space or covariates (SNPBRIM)
      MinFoldSampleSize = [0 0]; % Min samples per Outer Fold and Inner Fold (SNPBRIM)
      MaxIterations = 5; % Maximum number of BRIM iterations (SNPBRIM)
      Htest = true; % Perform Wilcoxon test om partial regression coefficients (SNPBRIM)
      RIPNormalize = true; % normalize RIP values based on group averages (SNPBRIM)
      StopCorr = 0.98; % Stopping correlation between subsequent iterations (SNPBRIM)
      RIPStatPerm = 1000; % Number of permutations in testing significance of RIP values (SNPBRIM)
      GroupDef = 'GG';% Definition of groups, GG = according to Test performed, GT = according to orginal genotypes (SNPBRIM)
   end
   methods % Constructor
        function obj = ParContainer(varargin)
            obj = obj@superClass(varargin{:});         
        end
   end
end

%        function out = get.AvgCov(obj)
%            if isempty(obj.Cov), out = []; return; end
%            out = nanmean(obj.Cov);
%        end

%        function out = get.ShapeSpace(obj)
%            out = obj.ShapeSpace;
%            if ~superClass.isH(out), out = []; end
%        end
%        function out = get.RedShapeSpace(obj)
%            out = obj.RedShapeSpace;
%            if ~superClass.isH(out), out = []; end
%        end
%        function out = get.n(obj)
%           if isempty(obj.ShapeSpace), out = 0; return; end
%           out = obj.ShapeSpace.n;
%        end 
%        function out = get.DepVar(obj)
%             if isempty(obj.ShapeSpace), out = [];return; end
%             out = obj.ShapeSpace.Tcoeff./repmat(obj.ShapeSpace.EigStd',obj.n,1);
%        end
%        function out = get.RedDepVar(obj)
%             if isempty(obj.RedShapeSpace), out = [];return; end
%             out = obj.RedShapeSpace.Tcoeff./repmat(obj.RedShapeSpace.EigStd',obj.n,1);
%        end

%    methods % Interface functions
%        function createReducedSpace(obj)
%            if isempty(obj.ShapeSpace), return; end
%            if isempty(obj.Cov), return; end
%            A = obj.Cov;
%            for i=1:1:size(A,2)% setting nanvalues to average
%                A(isnan(A(:,i)),i) = nanmean(A(:,i));
%            end
%            rm = PLSRShapeModel;
%            rm.Model = clone(obj.ShapeSpace);
%            rm.X = A;
%            update(rm);
%            % changing faces
%            RedDepVar = obj.ShapeSpace.Tcoeff;
%            avgA = mean(A);
%            for i=1:1:obj.n
%                RedDepVar(i,:) = changeX(rm,avgA,A(i,:),obj.ShapeSpace.Tcoeff(i,:));
%            end
%            tmpSpace = clone(obj.ShapeSpace);
%            tmpSpace.Tcoeff = RedDepVar;
%            shape = reconstructTraining(tmpSpace);
%            obj.RedShapeSpace = shapePCA;
%            obj.RedShapeSpace.RefScan = clone(obj.ShapeSpace.Average);
%            getAverage(obj.RedShapeSpace,shape);
%            getModel(obj.RedShapeSpace,shape);
%            reduceNrPc(obj.RedShapeSpace,obj.ShapeSpace.nrEV-size(A,2));
%        end
%    end
