classdef GBMOD < CONTMOD
    properties
       GBID = 0;
    end
    properties (Dependent = true)
       VarName;
       VarType;
    end
    properties (Hidden = true, Dependent = true)
       ID; 
    end
    methods % CONSTRUCTOR
        function obj = GBMOD(varargin)
            obj = obj@CONTMOD(varargin{:});         
        end
    end
    methods % GENERAL SETTING/GETTING
        function out = get.VarType(obj) %#ok<*MANU>
           out = 'GENETICBACKGROUND';
        end
        function out = get.VarName(obj)
            out = ['GB' num2str(obj.GBID)];
        end
        function out = get.ID(obj)
           out = 100; 
        end
    end
    methods % INTERFACING
    end
    methods (Static = true)    
    end
end