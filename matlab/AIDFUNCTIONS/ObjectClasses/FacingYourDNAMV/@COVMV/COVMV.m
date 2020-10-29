classdef COVMV < CONTMV
    properties
       VarName;
    end
    properties (Dependent = true)
       VarType;
    end
    properties (Hidden = true, Dependent = true)
       ID; 
    end
    methods % CONSTRUCTOR
        function obj = COVMV(varargin)
            obj = obj@CONTMV(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.VarType(obj) %#ok<*MANU>
           out = 'COVARIATE';
        end
        function out = get.ID(obj)
           if isempty(obj.VarName), out = 2; return; end
           switch obj.VarName
               case 'AGE'
                   out = 2;
               case 'W'
                   out = 3;
               case 'H'
                   out = 4;
           end
        end
    end
    methods % INTERFACING
    end
    methods (Static = true)    
    end
end