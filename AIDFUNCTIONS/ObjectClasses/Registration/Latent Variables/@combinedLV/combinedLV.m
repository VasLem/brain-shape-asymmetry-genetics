classdef combinedLV < LV
    % This is the abstract interface class for Transformation models
    properties
        LatentV;% the info of the deterministic function to update LV values
    end
    properties (Dependent = true)
        nrLV;
    end
    methods %Constructor
        function obj = combinedLV(varargin)
          obj = obj@LV(varargin{:});
          if nargin > 0
             Input = find(strcmp(varargin, 'LatentV'));if ~isempty(Input), obj.LatentV = varargin{Input+1}; end
          end
        end
    end
    methods% special setting & getting
        function nr = get.nrLV(obj)
            nr = length(obj.LatentV);
        end
        function out = get.LatentV(obj)
                 out = obj.LatentV;
                 if isempty(out), return; end
                 index = (1:length(out));
                 badindex = [];
                 for i=1:1:length(out)
                     if ~mySuperClass.isH(out{i}), badindex = [badindex i]; end %#ok<AGROW>
                 end
                 restindex = setdiff(index,badindex);
                 out = out(restindex);
        end
        function obj = set.LatentV(obj,in)
            if isempty(obj.LatentV), obj.LatentV = in; return; end
            for i=1:1:max(obj.nrLV,length(in))
                obj.LatentV{i} = in{i};
            end
        end
    end
    methods % InterFace functions
        function delete(obj)
            if isempty(obj.LatentV), return; end
            latentv = obj.LatentV;
            for i=1:1:obj.nrLV
                delete(latentv{i});
            end
        end
        function copy(obj,cobj)
                 copy@LV(obj,cobj);
                 if isempty(obj.LatentV), return; end
                 for i=1:1:obj.nrLV
                     cobj.LatentV{i} = clone(obj.LatentV{i});
                 end
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@LV(obj);
            if obj.nrLV == 0, struc.LatentV = {}; return; end
            for i=1:1:obj.nrLV
                struc.LatentV{i} = obj2struc(obj.LatentV{i});
            end
        end
        function obj = struc2obj(obj,struc)
                 obj = struc2obj@LV(obj,struc);
                 if isempty(struc.LatentV), obj.LatentV = {}; return; end
                 for i=1:1:length(struc.LatentV)
                     obj.LatentV{i} = struc2obj(eval(struc.LatentV{i}.Type),struc.LatentV{i});
                 end
        end
        function out = initialize(obj,Tmodel)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 initialize@LV(obj,Tmodel);
                 if isempty(obj.LatentV), return; end
                 for i=1:1:obj.nrLV
                     initialize(obj.LatentV{i},Tmodel);
                 end       
        end
        function out = update(obj,Tmodel)
                 if isempty(obj.LatentV), return; end
                 tmpout = ones(1,Tmodel.nrV);
                 for i=1:1:obj.nrLV
                     update(obj.LatentV{i},Tmodel);
                     tmpout = tmpout.*obj.LatentV{i}.Value;
                 end
                 if nargout == 1, out = tmpout; return; end
                 obj.Value = tmpout;
        end
        function out = clear(obj)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 clear@LV(obj);
                 if isempty(obj.LatentV), return; end
                 for i=1:1:obj.nrLV
                     clear(obj.LatentV{i});
                 end 
        end
        function postCompletePSetting(obj)
            if isempty(obj.LatentV), return; end
            for i=1:1:obj.nrLV;
                obj.LatentV{i}.CompleteP = obj.CompleteP;
            end
            %dummy
        end
    end
end % classdef

