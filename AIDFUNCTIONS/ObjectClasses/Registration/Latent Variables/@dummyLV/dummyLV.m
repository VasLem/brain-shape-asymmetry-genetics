classdef dummyLV < LV
    % This is an interface class for Dummy Latent Variables, never
    % change!!!
    methods %Constructor
        function obj = dummyLV(varargin)
          obj = obj@LV(varargin{:});
        end
    end       
    methods % InterFace functions
        function out = update(obj,Tmodel) %#ok<INUSD>
            % the dummy LV does not update the LV values
            obj.Value = ones(1,Tmodel.nrV);
            if nargout == 1, out = obj.Value;end
        end
        
    end
end % classdef

