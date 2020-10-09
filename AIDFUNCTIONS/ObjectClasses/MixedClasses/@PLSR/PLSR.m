classdef PLSR < PLSRSuper
    properties
        Y = [];% Model to regress against
        X = [];% Predictors
        XNames = [];% PredictorNames
    end
    methods % Constructor
        function obj = PLSR(varargin)
            obj = obj@PLSRSuper(varargin{:});         
        end
    end
end