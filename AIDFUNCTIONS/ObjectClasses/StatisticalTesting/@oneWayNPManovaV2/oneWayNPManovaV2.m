classdef oneWayNPManovaV2 < superClass % V2 is much faster than previous implementation (oneWayNPManova)
    properties
        G = [];% cell with the default two groups, observations (rows) X variables (columns)
        F = [];% F-statistic
        pF = [];% P-value of F statistic
        t = 1000;% number of permutations
        DistType = 'euclidean';% Type of distances between observations
    end
    properties (Dependent = true)
        SST;% Total Sums of Squares
        SSW;% Within Sums of Squares
        SSB;% Between Sums of Squares
        nrWG;% Array of number of observations per group
        nrG;% Number of groups
        N;% Total Number of observations
        V;% Dimension of variables
    end
    properties (Hidden = true, Dependent = true);
        GT;% All groups accumulated
        arr;% Aid array to extract within group D
    end
    properties (Hidden = true); % computationally friendly hidden properties
        SSTH = [];% Non normalized SST 
        SSWH = [];% Non normalized SSW
        SSBH = [];% Non normalized SSB
        arrH = [];% Aid array to extract within group D
        GTH = [];% All groups accumulated
        nrWGH = [];% Array of number of observations per group
        nrGH = [];% Number of groups
        NH = [];% Total Number of observations
        VH = [];% Dimension of Variables
        DT = [];% Total squared Distance Matrix
    end
    methods % Constructor
        function obj = oneWayNPManovaV2(varargin)
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
        function out = get.arr(obj)
           if ~isempty(obj.arrH), out = obj.arrH; return; end
           if isempty(obj.nrWG), out = []; return; end
           out = [0 obj.nrWG 0];
           obj.arrH = out;
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
        function out = get.SST(obj)
            if isempty(obj.SSTH), out = []; return; end
            out = obj.SSTH/(obj.N-1);
        end
        function out = get.SSW(obj)
            if isempty(obj.SSWH), out = []; return; end
            out = obj.SSWH/(obj.N-obj.nrG);
        end
        function out = get.SSB(obj)
            if isempty(obj.SSBH), out = []; return; end
            out = obj.SSBH/(obj.nrG-1);
        end
    end
    methods % Interface functions
        function clear(obj)
             obj.arrH = [];% Aid array to extract within group D
             obj.GTH = [];% All groups accumulated
             obj.nrWGH = [];% Array of number of observations per group
             obj.nrGH = [];% Number of groups
             obj.NH = [];% Total Number of observations
             obj.VH = [];% Dimension of variables
             obj.DT = [];% Total squared Distance Matrix
             obj.SSTH = [];% Non normalized SST 
             obj.SSWH = [];% Non normalized SSW
             obj.SSBH = [];% Non normalized SSB
        end
        function initialize(obj)
             clear(obj);
             getDT(obj);
        end
        function getDT(obj)
            obj.DT = squareform(pdist(obj.GT,obj.DistType).^2);
        end
        function out = getSST(obj,DT)
                 if nargin<2, DT = obj.DT; end
                 UD = triu(DT);
                 out = sum(UD(UD>0))/obj.N;
                 if nargout==1, return; end
                 obj.SSTH = out;
        end
        function out = getSSW(obj,DT)
                 if nargin<2, DT = obj.DT; end
                 out = zeros(1,obj.nrG);
                 for k=1:obj.nrG
                     s1 = sum(obj.arr(1:k)); %#ok<*PFBNS>
                     s2 = sum(obj.arr(1:k+1));
                     Dred = triu(DT(s1+1:s2,s1+1:s2));
                     out(k) = sum(Dred(Dred>0))/obj.nrWG(k);  
                 end
                 out = sum(out);
                 if nargout==1, return; end
                 obj.SSWH = out;
        end
        function getF(obj)
                 if isempty(obj.SSB), getSSB(obj); end
                 obj.F = obj.SSB/obj.SSW;
        end
        function out = getSSB(obj)
                 out = (obj.SSTH-obj.SSWH);
                 if nargout==1, return; end
                 obj.SSBH = out;
        end
        function getPValue(obj)
                 if obj.t==0, return; end
                 FCount = false(1,obj.t);
                 parfor i=1:obj.t
                     ind = randperm(obj.N);
                     DTfor = obj.DT(ind,ind);
                     SSWFor = getSSW(obj,DTfor);
                     Ffor = ((obj.SSTH-SSWFor)/(obj.nrG-1))/(SSWFor/(obj.N-obj.nrG));
                     FCount(i) = Ffor>=obj.F;
                 end
                 obj.pF = sum(FCount)/obj.t;         
        end
        function run(obj)
                 % Initialize, make sure hidden variables are cleared and
                 % compute Squared Distance Matrix
                 initialize(obj);
                 % Compute total Sum Of Squared Differences
                 getSST(obj);
                 % Compute Within Group Sum Of Squared Differences
                 getSSW(obj);
                 % Compute Between Group Sum of Squared Differences
                 getSSB(obj);
                 % Compute F statistic
                 getF(obj);
                 % Compute P Value Under Permutation
                 getPValue(obj);        
        end
    end
end