classdef distanceSM < SM
    % This is the abstract interface class for similarity measures working
    % on any type of distances to the target surface
    properties
        Fields = {'Vertices' 'Distance'};% Cannot be changed
        Evaluation = [];
        Derivative = [];
        nrSM = 1;% Cannot be changed
    end
    methods %Constructor
        function obj = distanceSM(varargin)
            obj = obj@SM(varargin{:});
        end
    end
    methods %Special Setting & Getting
        function obj = set.Fields(obj,fields) %#ok<INUSD>
            %dummy
        end
        function obj = set.nrSM(obj,nr) %#ok<INUSD>
            %dummy
        end
    end
end % classdef