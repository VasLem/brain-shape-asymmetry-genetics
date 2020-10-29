classdef PLSRC2X < PLSRSuper
    properties
        Parent = [];
    end
    properties (Dependent = true)
        X;% Predictors
        XNames;% PredictorNames
        Y;
    end
    methods % Constructor
        function obj = PLSRC2X(varargin)
            obj = obj@PLSRSuper(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.X(obj)
            out = obj.Parent.C;
        end
        function out = get.XNames(obj)
            out = obj.Parent.CNames;
        end
        function out = get.Y(obj)
           out = obj.Parent.X; 
        end
    end
    methods % Concrete Basic Interface & Conversion functions
       % nothing to implement     
    end
end