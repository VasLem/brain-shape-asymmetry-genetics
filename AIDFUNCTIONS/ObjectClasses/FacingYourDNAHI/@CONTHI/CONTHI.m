classdef CONTHI < PPMHI
    properties
       VarName;
    end
    properties (Dependent = true)
       CONTV; 
       VarType;
    end
    properties (Hidden = true)
        NanVal = nan;
    end
    properties(Hidden = true, Dependent = true)
        ID;
    end
    methods % CONSTRUCTOR
        function obj = CONTHI(varargin)
            obj = obj@PPMHI(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function set.CONTV(obj,in)
                 obj.XX = in;
        end
        function out = get.CONTV(obj)
                 out = obj.XX;
        end
        function out = get.VarName(obj) %#ok<*MANU>
           if isempty(obj.VarName), out = 'BIN'; return; end
           out = obj.VarName; 
        end
        function out = get.VarType(obj)
           out = 'CONTINUOUS VARIABLE';
        end
        function out = get.ID(obj)
            out = 999;
        end
    end
    methods % CLASSIFIER
    end
    methods % BIOMETRICS
    end
    methods % INTERFACING
        function out = XX2XG(obj,in)
                 out = in;
        end
        function out = XG2GG(obj,in)
                 out = in;
        end
    end
    methods (Static = true)    
    end
end