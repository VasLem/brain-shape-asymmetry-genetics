classdef ClPtsDistanceSM < distanceSM
    properties
        Update = true;% Are Closest Points updated or not
        Index = [];% Index in the Target surface of the points being the closest points of the Floating surface 
        Distance = [];% Distance between closest Points ||p1-p2||
        Difference = [];% Difference between closest Points (p1-p2)
        %Direction = [];% Direction between closest Points
        RelatedSM = 'ClPtsAngleSM';% Related SM, used to reduce computation when SM is combined with related SM
        RelatedFields = {'Index' 'Distance' 'Difference'};% Related Fields wit related SM
    end
    properties (Dependent = true)
        TargetInfo;
        Direction;% Dependent on Difference;
    end
    methods %Constructor
        function obj = ClPtsDistanceSM(varargin)
            obj = obj@distanceSM(varargin{:});
            if nargin > 0
               Input = find(strcmp(varargin, 'Update'));if ~isempty(Input), obj.Update = varargin{Input+1}; end
               Input = find(strcmp(varargin, 'Index'));if ~isempty(Input), obj.Index = varargin{Input+1}; end
               Input = find(strcmp(varargin, 'Distance'));if ~isempty(Input), obj.Distance = varargin{Input+1}; end
               Input = find(strcmp(varargin, 'Difference'));if ~isempty(Input), obj.Difference = varargin{Input+1}; end
            end
        end
    end
    methods% special setting & getting
        function info = get.TargetInfo(obj)
            switch class(obj.Target)
                case {'meshObj' 'LMObj'}
                    info = obj.Target.Vertices;
                otherwise
                    info = obj.Target;% Target must be list of points
            end          
        end
        function out = get.Direction(obj)
            out = obj.Difference./repmat(obj.Distance,3,1);
        end
    end
    methods % InterFace functions      
        function struc = obj2struc(obj)
            % converting relevant information
            struc = obj2struc@distanceSM(obj);
            struc.Update = obj.Update;
            struc.Index = obj.Index;
            struc.Distance = obj.Distance;
            struc.Difference = obj.Difference;            
        end
        function obj = struc2obj(obj,struc)
            obj = struc2obj@distanceSM(obj,struc);
            obj.Update = struc.Update;
            obj.Index = struc.Index;
            obj.Distance = struc.Distance;
            obj.Difference = struc.Difference;
        end
        function copy(obj,cobj) %#ok<INUSD>
                 copy@distanceSM(obj,cobj);
                 cobj.Update = obj.Update;
                 cobj.Index = obj.Index;
                 cobj.Distance = obj.Distance;
                 cobj.Difference = obj.Difference;
        end
        function out = clear(obj)
            if nargout == 1
               obj = clone(obj);
               out = obj;
            end
            clear@distanceSM(obj);
            if obj.Update, obj.Index = []; end
            obj.Distance = [];
            obj.Difference = [];
        end
    end
end % classdef

