classdef detLV < LV
    % This is the abstract interface class for Deterministic latent
    % Variables
    properties (Abstract = true) 
        TH;% Threshold value on function evaluations to decidce inlier vs outlier
    end
    properties (Abstract = true, Dependent = true)
        Function;% the info of the deterministic function to update LV values, typically dependent on Target Surface (stored in CompleteP->InlierP->SM)
    end
    methods %Constructor
        function obj = detLV(varargin)
          obj = obj@LV(varargin{:});
          if nargin > 0
             Input = find(strcmp(varargin, 'TH'));if ~isempty(Input), obj.TH = varargin{Input+1}; end
          end
        end
    end        
    methods % InterFace functions 
        function copy(obj,cobj)
                 copy@LV(obj,cobj);
                 cobj.TH = obj.TH;
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@LV(obj);
            struc.TH = obj.TH;
        end
        function obj = struc2obj(obj,struc)
                 obj = struc2obj@LV(obj,struc);
                 obj.TH = struc.TH;
        end
    end
end % classdef

