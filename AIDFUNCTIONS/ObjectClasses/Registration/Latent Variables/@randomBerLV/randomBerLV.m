classdef randomBerLV < randomLV
    % This is the abstract interface class random Bernoulli Distributed
    % Lantent Variables
    methods %Constructor
        function obj = randomBerLV(varargin)
          obj = obj@randomLV(varargin{:});
        end
    end        
end % classdef

