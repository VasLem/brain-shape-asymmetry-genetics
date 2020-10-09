classdef LV < mySuperClass
    % This is the abstract interface class for Random LV, typically
    % combined with a InlierP and OutlierP in a CompleteP, to improve
    % robustness against outliers
    properties 
        CompleteP = []; %Link to Parent being a complete Process
        Value = [];% LV Value
    end
    properties (Dependent = true)
        sumV;
        nrV;
    end
    methods %Constructor
        function obj = LV(varargin)
          if nargin > 0
             Input = find(strcmp(varargin, 'CompleteP'));if ~isempty(Input), obj.CompleteP = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Value'));if ~isempty(Input), obj.Value = varargin{Input+1}; end
          end
        end
    end  
    methods % Special Setting and Getting
        function s = get.sumV(obj)
            s = sum(obj.Value);
        end
        function nr = get.nrV(obj)
            nr = length(obj.Value);
        end
        function cp = get.CompleteP(obj)
            cp = obj.CompleteP;
            if~mySuperClass.isH(cp), cp = []; end
        end
        function obj = set.CompleteP(obj,completep)
            obj.CompleteP = completep;
            postCompletePSetting(obj);% relevant for combined LV
        end
        function postCompletePSetting(obj) %#ok<INUSD>
            %dummy
        end
    end       
    methods % InterFace functions 
        function copy(obj,cobj)
                 copy@mySuperClass(obj,cobj);
                 cobj.CompleteP = obj.CompleteP;% parents are not cloned
                 cobj.Value = obj.Value;
        end
        function out = clear(obj)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 obj.Value = [];
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@mySuperClass(obj);
            struc.Value = obj.Value;
        end
        function obj = struc2obj(obj,struc)
                 obj.Value = struc.Value;
        end
        function out = initialize(obj,Tmodel)
                 if nargout == 1
                    obj = clone(obj);
                    out = obj;
                 end
                 obj.Value = ones(1,Tmodel.nrV);            
        end
    end
    methods (Abstract = true) %Abstract Interfaces Functions
            out = update(obj,Tmodel);
    end
end % classdef

