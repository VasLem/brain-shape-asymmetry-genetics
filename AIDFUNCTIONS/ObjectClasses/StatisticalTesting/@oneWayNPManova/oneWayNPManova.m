classdef oneWayNPManova < superClass
    properties
        D = [];% Distance Matrix
        n = [];% Array of number of observations per group
        SST = []% Total Sums of Squares
        SSW = [];% Within Sums of Squares
        SSB = [];% Between Sums of Squares
        F = [];% F-statistic
        pF = [];% P-value of F statistic
        t = 100;% number of permutations
    end
    properties (Dependent = true)
        g;% Number of groups
        N;% Total Number of observations
        D2;% squared distance matrix
    end
    properties (Hidden = true, Dependent = true);
        nr;
    end
    methods % Constructor
        function obj = oneWayNPManova(varargin)
          obj = obj@superClass(varargin{:});
        end
    end
    methods % Special Getting and Setting
        function out = get.g(obj)
           out = length(obj.n); 
        end
        function out = get.N(obj)
           out = sum(obj.n); 
        end
        function out = get.nr(obj)
           out = [0 obj.n 0]; 
        end
        function out = get.D2(obj)
            out = obj.D.^2;
        end
    end
    methods % Interface functions
        function getSST(obj)
                 UD = triu(obj.D.^2);
                 obj.SST = sum(UD(UD>0))/obj.N;
        end
        function getSSW(obj)
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
                     %out(k) = out(k)/obj.n(k);
                 end
                 %obj.SSW = sum(out);
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
                 getSST(obj);
                 getSSW(obj);
                 getF(obj);
                 getPValue(obj);        
        end
    end    
end