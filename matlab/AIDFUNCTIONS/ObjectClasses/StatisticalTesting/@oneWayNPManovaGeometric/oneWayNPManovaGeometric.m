classdef oneWayNPManovaGeometric < superClass % V2 is much faster than previous implementation (oneWayNPManova)
    properties
        G = [];% cell with the default two groups, observations (rows) X variables (columns)
        SST = []% Total Sums of Squares
        SSW = [];% Within Sums of Squares
        SSB = [];% Between Sums of Squares
        F = [];% F-statistic
        pF = [];% P-value of F statistic
        t = 100;% number of permutations
        DistType = 'euclidean';% Type of distances between observations
    end
    properties (Dependent = true)
        nrWG;% Array of number of observations per group
        nrG;% Number of groups
        N;% Total Number of observations
        V;% Dimension of variables
    end
    properties (Hidden = true, Dependent = true);
        GT;% All groups accumulated
        GAvg;% Group averages
        GTAvg;% Total average
    end
    properties (Hidden = true); % computationally friendly hidden properties
        GTH = [];% All groups accumulated
        nrWGH = [];% Array of number of observations per group
        nrGH = [];% Number of groups
        NH = [];% Total Number of observations
        VH = [];% Dimension of Variables
        GAvgH = [];% within group averages
        GTAvgH = [];% Total average
    end
    methods % Constructor
        function obj = oneWayNPManovaGeometric(varargin)
          obj = obj@superClass(varargin{:});
        end
    end
    methods % Special Getting and Setting
        function out = get.nrG(obj)
           if ~isempty(obj.nrGH), out = obj.nrGH; return; end
           out = length(obj.G);
           if out>0, obj.nrGH = out; end
        end
        function out = get.nrWG(obj)
            if ~isempty(obj.nrWGH), out = obj.nrWGH; return; end
            if obj.nrG==0, out = []; return; end
            out = zeros(1,obj.nrG);
            for i=1:1:obj.nrG
                out(i) = size(obj.G{i},1);
            end
            obj.nrWGH = out;
        end
        function out = get.N(obj)
           if ~isempty(obj.NH), out = obj.NH; return; end 
           out = sum(obj.nrWG);
           if out>0, obj.NH = out; end
        end
        function out = get.V(obj)
            if ~isempty(obj.VH), out = obj.VH; return; end
            if obj.nrG==0, out = []; return; end
            if isempty(obj.G{1}), out = []; return; end
            out = size(obj.G{1},2);
            obj.VH = out;
        end
        function out = get.GT(obj)
           if ~isempty(obj.GTH), out = obj.GTH; return; end 
           out = [];
           if obj.nrG == 0, return; end
           for i=1:1:obj.nrG
               out = [out; obj.G{i}]; %#ok<*AGROW>
           end
           obj.GTH = out;
        end
        function out = get.GAvg(obj)
           if ~isempty(obj.GAvgH), out = obj.GAvgH; return;end
           if obj.nrG==0; out = []; return; end
           out = cell(obj.nrG,obj.V);
           for i=1:1:obj.nrG
               out{i} = mean(obj.G{i});
           end
           obj.GAvgH = out;
        end
        function out = get.GTAvg(obj)
           if ~isempty(obj.GTAvgH), out = obj.GTAvgH; return;end
           if obj.nrG==0; out = []; return; end
           out = mean(obj.GT);
           obj.GTAvgH = out;
        end
    end
    methods % Interface functions
        function clear(obj)
             obj.GTH = [];% All groups accumulated
             obj.nrWGH = [];% Array of number of observations per group
             obj.nrGH = [];% Number of groups
             obj.NH = [];% Total Number of observations
             obj.VH = [];% Dimension of variables
             obj.GAvgH = [];% group averages
             obj.GTAvgH = [];% Total average
        end
        function initialize(obj)
             clear(obj);
        end
        function out = getSST(obj)
                 out = oneWayNPManovaGeometric.compSST(obj.GT,obj.GTAvg,obj.DistType,obj.N-1);
                 if nargout==1, return; end
                 obj.SST = out;
        end
        function out = getSSW(obj)
                 out = oneWayNPManovaGeometric.compSSW(obj.G,obj.GAvg,obj.nrG,obj.DistType,obj.N-obj.nrG);
                 if nargout==1, return; end
                 obj.SSW = out;
        end
        function out = getF(obj)
                 if isempty(obj.SSB), getSSB(obj); end
                 if isempty(obj.SSW), getSSW(obj); end
                 out = obj.SSB/obj.SSW;
                 if nargout==1, return; end
                 obj.F = out;
        end
        function out = getSSB(obj)
                 out = oneWayNPManovaGeometric.compSSB(obj.GAvg,obj.GTAvg,obj.nrG,obj.nrWG,obj.DistType,obj.nrG-1);
                 if nargout==1, return; end
                 obj.SSB = out;
        end
        function getPValue(obj)
                 obj.pF = 'TO BE IMPLEMENTED';
%                  FCount = false(1,obj.t);
%                  disp('Permuting');
%                  tic;
%                  parfor i=1:obj.t
%                      ind = randperm(obj.N);
%                      Objperm = oneWayNPManova;
%                      Objperm.D = obj.D(ind,ind);
%                      Objperm.n = obj.n;
%                      getF(Objperm);
%                      FCount(i) = Objperm.F>=obj.F; 
%                  end
%                  obj.pF = sum(FCount)/obj.t;
%                  toc;          
        end
        function run(obj)
                 % Initialize, make sure hidden variables are cleared and
                 % compute Squared Distance Matrix
                 initialize(obj);
                 % Compute total Sum Of Squared Differences
                 getSST(obj);
                 % Compute Within Group Sum Of Squared Differences
                 getSSW(obj);
                 % Compute F statistic
                 getF(obj);
                 % Compute P Value Under Permutation
                 getPValue(obj);        
        end
    end
    methods (Static = true)
        function out = compSST(GT,GTAvg,DistType,degrees)
            out = sum(pdist2(GT,GTAvg,DistType).^2)/degrees;
        end
        function out = compSSW(G,GAvg,nrG,DistType,degrees)
            out = zeros(1,nrG);
            for i=1:1:nrG
                 out(i) = sum(pdist2(G{i},GAvg{i},DistType).^2);
            end
            out = sum(out)/degrees;
        end
        function out = compSSB(GAvg,GTAvg,nrG,nrWG,DistType,degrees)
            out = zeros(1,nrG);
            for i=1:1:nrG
                   out(i) = nrWG(i)*(pdist2(GAvg{i},GTAvg,DistType).^2);
            end
            out = sum(out)/degrees;
        end
    end
end