classdef detInputLV < detLV
    % This is the abstract interface class for Transformation models
    properties
        TH = 0;
        PriorValues = [];
    end
    properties (Dependent = true)
        Function;% the info of the deterministic function to update LV values
    end
    methods %Constructor
        function obj = detInputLV(varargin)
          obj = obj@detLV(varargin{:});
        end
    end        
    methods % Special Setting and Getting
        function F = get.Function(obj)
            F = obj.PriorValues;
        end
        function obj = set.Function(obj,input) %#ok<INUSD>
            % Dummy      
        end       
    end
end % classdef

