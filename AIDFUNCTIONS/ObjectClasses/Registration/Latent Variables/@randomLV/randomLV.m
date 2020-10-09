classdef randomLV < LV
    % This is the abstract interface class Random LV
    methods %Constructor
        function obj = randomLV(varargin)
          obj = obj@LV(varargin{:});
        end
    end
end % classdef

