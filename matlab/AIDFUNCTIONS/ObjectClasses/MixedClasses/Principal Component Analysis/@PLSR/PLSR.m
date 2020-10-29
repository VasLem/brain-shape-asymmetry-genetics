classdef PLSR < superClass
    properties
        Resp = [];% Responses
        Pred = [];% Predictors
        PredNames = [];% PredictorNames
        nrcmp = [];% nrcmp
        M = [];
        PcVar = [];
    end
    properties (Dependent = true)
        RespAvg;
        RespStd;
        PredAvg;
        PredStd;
    end
    methods % Constructor
        function obj = PLSR(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.nrEV(obj)
            out = length(obj.EigVal);
        end
        function obj = set.EigStd(obj,in)
            obj.EigVal = in.^2;
        end
    end
    methods % Concrete Basic Interface & Conversion functions
       function out = Coeff2Vec(obj,in)
            % Convert given PCA coeff into a Vector
            if size(in,1)==1, in = in'; end
            if obj.Centering
                out = obj.AvgVec + obj.EigVec*in;
            else
                out = obj.EigVec*in;
            end
       end
    end
end