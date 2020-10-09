classdef twoWayNPManova < superClass
    properties
        D = [];% Distance Matrix
        lA = [];% number of A levels
        lB = [];% number of B levels
        SST = []% Total Sums of Squares
        SSA = [];% Within Sums of Squares of Factor A
        SSB = [];% Within Sums of Squares of Factor B
        SSR = [];% residual sums of squares
        SSAB = [];% Interaction Sums of Squares
        t = 1000;% number of permutations
    end
    properties (Dependent = true)
        g;% Number of groups
        N;% Total Number of observations
        nA;
        dfA;
        dfB;
        dfAB;
        dfR;
        dfT;
        nB;
        MSA;
        MSB;
        MSAB;
        MSR;
        FA;
        FB;
        FAB;
        %nR; 
    end
    properties (Hidden = true, Dependent = true);
        nr;
    end
    methods % Constructor
        function obj = twoWayNPManova(varargin)
          obj = obj@superClass(varargin{:});
        end
    end
    methods % Special Getting and Setting
        function out = get.g(obj)
           out = length(obj.n); 
        end
        function out = get.N(obj)
           out = size(obj.D,1); 
        end
        function out = get.nr(obj)
           out = [0 obj.n 0]; 
        end
        function out = get.nA(obj)
            out = obj.N/obj.lA;
        end
        function out = get.nB(obj)
            out = (obj.N/obj.lA)/obj.lB;
        end
        function out = get.dfA(obj)
                 out = obj.lA-1; 
        end
        function out = get.dfB(obj)
                 out = obj.lB-1;
        end
        function out = get.dfAB(obj)
                 out = obj.dfA*obj.dfB;
        end
        function out = get.dfT(obj)
                 out = obj.N-1; 
        end
        function out = get.dfR(obj)
                 out = obj.N-obj.lA*obj.lB; 
        end
        function out = get.MSA(obj)
                 out = obj.SSA/obj.dfA; 
        end
        function out = get.MSB(obj)
                 out = obj.SSB/obj.dfB; 
        end
        function out = get.MSAB(obj)
                 out = obj.SSAB/obj.dfAB; 
        end
        function out = get.MSR(obj)
                 out = obj.SSR/obj.dfR; 
        end
        function out = get.FA(obj)
                 out = obj.MSA/obj.MSR;
        end
        function out = get.FB(obj)
                 out = obj.MSB/obj.MSR;
        end
        function out = get.FAB(obj)
                 out = obj.MSAB/obj.MSR;
        end
            
            
    end
    methods % Interface functions
        function out = inA(obj,i)
                 out = ceil(i/obj.nA);
        end
        function out = inB(obj,i)
                 i = i-((inA(obj,i)-1)*obj.nA);
                 out = ceil(i/obj.nB);
        end
        function getSS(obj)
                 T = 0;
                 A = 0;
                 B = 0;
                 R = 0;
                 N = obj.N;
                 for i=1:1:(N-1)
                    for j=(i+1):1:N
                        d = obj.D(i,j)^2;
                        T = T+d;
                        if inA(obj,i)==inA(obj,j)
                           A = A+d;
                        end
                        if inB(obj,i)==inB(obj,j)
                           B = B+d;
                           if inA(obj,i)==inA(obj,j)
                              R = R+d;
                           end
                        end
                    end
                 end
                 obj.SST = T/N;
                 obj.SSA = A/(obj.nB*obj.lB);
                 obj.SSB = B/(obj.nB*obj.lA);
                 obj.SSR = R/obj.nB;
                 obj.SSAB = obj.SST-obj.SSA-obj.SSB-obj.SSR;
        end
        function getSST(obj)
                 out = 0;
                 N = obj.N;
                 for i=1:1:(N-1)
                    for j=(i+1):1:N
                        out = out+obj.D(i,j)^2;
                    end
                 end
                 obj.SST = out/N;
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
                 getSST(obj);
                 getSSW(obj);
                 getF(obj);
                 getPValue(obj);        
        end
    end    
end