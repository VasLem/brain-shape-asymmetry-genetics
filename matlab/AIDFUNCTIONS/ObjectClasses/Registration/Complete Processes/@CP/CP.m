classdef CP < mySuperClass
    % This is the abstract interface class Complete Processes
    properties
        Evaluation = [];% Complete Process negative log evaluations over all floating surface points
        Derivative = [];% Complete Process gradients in floating surface points
    end
    properties (Abstract = true)
        InlierP; % link to Inlier Process
        OutlierP; % link to Outlier Process
        LatentV; % link to Latent variable
        Target;% Target handle is stored in CompleteProcess.InlierP.Smeasure.Target;
        Smeasure;% Smeasure handle is stored in CompleteProcess.InlierP
        Kappa;% Statistical significance of Inlierp vs Outlierp
        Fields;% Fields on which completeP operates
        B;% B values = LatentV.Value
        Correspondence;% Dependent on Correspondence stored in CompleteProcess.InlierProces.Smeasure.Correspondence
        Temp;
        FastEval; % Fast evaluation of Inlierprocess
    end
    properties (Dependent = true)
        nrFields;
        %nrV;
    end
    methods %Constructor
        function obj = CP(varargin)
          obj = obj@mySuperClass(varargin{:});
          if nargin > 0
             Input = find(strcmp(varargin, 'InlierP'));if ~isempty(Input), obj.InlierP = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Kappa'));if ~isempty(Input), obj.Kappa = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'OutlierP'));if ~isempty(Input), obj.OutlierP = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'LatentV'));if ~isempty(Input), obj.LatentV = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Evaluation'));if ~isempty(Input), obj.Evaluation = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Derivative'));if ~isempty(Input), obj.Derivative = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Smeasure'));if ~isempty(Input), obj.Smeasure = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Target'));if ~isempty(Input), obj.Target = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Correspondence'));if ~isempty(Input), obj.Correspondence = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Temp'));if ~isempty(Input), obj.Correspondence = varargin{Input+1}; end
          end
        end
    end
   methods % Special Getting & Setting
       function nr = get.nrFields(obj)
            nr = length(obj.Fields);
       end
   end
   methods % InterFace functions 
        function copy(obj,cobj)
                 copy@mySuperClass(obj,cobj);
                 cobj.Evaluation = obj.Evaluation;
                 cobj.Derivative = obj.Derivative;          
        end
        function out = clear(obj)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 obj.Evaluation = [];
                 obj.Derivative = [];
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@mySuperClass(obj);
            struc.Evaluation = obj.Evaluation;
            struc.Derivative = obj.Derivative;
        end
        function obj = struc2obj(obj,struc)
                 obj.Evaluation = struc.Evaluation;
                 obj.Derivative = struc.Derivative; 
        end
   end
   methods (Abstract = true)
       out = eval(obj,Tmodel);
       out = derivate(obj,Tmodel);
       out = update(obj,Tmodel);
       out = initialize(obj,Tmodel);
       out = updateCorrespondence(obj,Tmodel);       
   end
end % classdef

