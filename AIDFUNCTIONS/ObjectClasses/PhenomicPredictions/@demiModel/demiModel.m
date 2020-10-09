classdef demiModel < superClass
    properties
        Tag;% Tag on which the demiModel operates, needed to link to input data in an excell file
        TrX;% Training X values, these are tyically, original predictor variables
        XType;% The type of X values, Categorical or Continous
        TrY;% Training Y values, these are the Facial parameters (best to use PCA coefficients)
        Acc;% Accuracy of the demiModel using standard tester
        Level;% Genomic or Genotype, what is the level of DNA that the demi model is active (used to discriminate demimodels used in the creation of the base face or not)
    end
    properties (Dependent = true)
       nrTr;% nr Training faces
    end
    methods % Constructor
        function obj = demiModel(varargin)
            obj = obj@superClass(varargin{:});         
        end
    end
    methods % Special Setting & Getting
         function out = get.nrTr(obj)
            out = length(obj.TrX);
        end
    end
    methods (Abstract = true) % Abstract Basic Interface & Conversion functions
        % Need to be made explicit depending on the type of data, eg shape, texture,...
        out = getScore(obj,in,val);
    end
    methods % General Interface Functions        
        function getAcc(obj)
            errors = zeros(1,obj.nrTr);
            for i=1:obj.nrTr
               errors(i) = getScore(obj,obj.TrY(i,:)',obj.TrX(i));
            end
            errors = errors(~isnan(obj.TrX));
            obj.Acc = sqrt(mean(errors.^2));
        end
    end
end