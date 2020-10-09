classdef PLSRReduced < PLSRSuper
    properties
        Parent = [];
    end
    properties (Dependent = true)
        Y;% Model to regress against
        X;% Predictors
        XNames;% PredictorNames
    end
    methods % Constructor
        function obj = PLSRReduced(varargin)
            obj = obj@PLSRSuper(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.Y(obj)
           out = obj.Parent.YCResiduals; 
        end
        function out = get.X(obj)
           out = obj.Parent.XCResiduals; 
        end
        function out = get.XNames(obj)
            out = obj.Parent.XNames;
        end
    end
    methods % Concrete Basic Interface & Conversion functions
       % nothing to implement
    end
end