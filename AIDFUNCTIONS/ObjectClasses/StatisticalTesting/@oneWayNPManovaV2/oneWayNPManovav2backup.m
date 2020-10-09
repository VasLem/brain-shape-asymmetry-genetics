classdef oneWayNPManovaV2 < superClass % V2 is much faster than previous implementation (oneWayNPManova)
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
    end
    properties (Hidden = true, Dependent = true);
        GT;% All groups accumulated
        arr;% Aid array to extract within group D
    end
    properties (Hidden = true); % computationally friendly hidden properties
        arrH = [];% Aid array to extract within group D
        GTH = [];% All groups accumulated
        nrWGH = [];% Array of number of observations per group
        nrGH = [];% Number of groups
        NH = [];% Total Number of observations
        DT = [];% Total squared Distance Matrix
        DW = [];% cell with within group squared distance matrices
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
            if obj.nrG==0; out = []; return; end
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
    end
    methods % Interface functions
        function clear(obj)
             obj.nrH = [];% Aid array to extract within group D
             obj.GTH = [];% All groups accumulated
             obj.nrWGH = [];% Array of number of observations per group
             obj.nrGH = [];% Number of groups
             obj.NH = [];% Total Number of observations
             obj.DT = [];% Total squared Distance Matrix
             obj.DW = [];% cell with within group squared distance matrices
        end
        function getDT(obj)
            obj.DT = squareform(pdist(obj.GT,obj.DistType).^2);
        end
        function out = getSST(obj)
                 UD = triu(obj.DT);
                 out = sum(UD(UD>0))/obj.N;
                 if nargout==1, return; end
                 obj.SST = out;
        end
        function out = getSSW(obj)
                 g = obj.g;
                 nr = obj.nr;
                 out = zeros(1,obj.g);
                 for k=1:g
                     s1 = sum(nr(1:k)); %#ok<*PFBNS>
                     s2 = sum(nr(1:k+1));
                     Dred = obj.D(s1+1:s2,s1+1:s2);
                     ng = obj.n(k);
                     for i=1:1:(ng-1)
                        for j=(i+1):1:ng
                            out(k) = out(k)+Dred(i,j)^2;
                        end
                     end  
                 end
                 obj.SSW = sum(out)/obj.n(1);
        end
        function getSSB(obj)
                 obj.SSB = obj.SST-obj.SSW;
        end
        function getF(obj)
                 if isempty(obj.SST), getSST(obj); end
                 if isempty(obj.SSW), getSSW(obj); end
                 obj.F = ((obj.SST-obj.SSW)/(obj.g-1))/...
                         (obj.SSW/(obj.N-obj.g));
        end
        function getPValue(obj)
                 FCount = false(1,obj.t);
                 disp('Permuting');
                 tic;
                 parfor i=1:obj.t
                     ind = randperm(obj.N);
                     Objperm = oneWayNPManova;
                     Objperm.D = obj.D(ind,ind);
                     Objperm.n = obj.n;
                     getF(Objperm);
                     FCount(i) = Objperm.F>=obj.F; 
                 end
                 obj.pF = sum(FCount)/obj.t;
                 toc;
                 
        end
        function perform(obj)
                 % Initialize, make sure hidden variables are cleared
                 clear(obj);
                 % Compute squared distance matrix
                 getDT(obj);
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
end