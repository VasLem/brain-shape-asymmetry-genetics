classdef floatingS < hgsetget
    % This is the abstract interface class for similarity measures
    properties
        Type = 'floatingS';% Type of object
        Surface = [];% Floating surface information (meshObj,LMObj or a pointlist);
        Tmodel = [];% Transformation model
        B = [];% Floating surface confidence values
    end
    properties (Dependent = true)
        Vertices;
        nrV;
        SurfaceClass;% can be a meshObj,LMObj or a pointlist
        sumB;
    end
    methods %Constructor
        function obj = floatingS(varargin)
            if nargin > 0
             Input = find(strcmp(varargin, 'Surface'));if ~isempty(Input), obj.Surface = varargin{Input+1}; end
             Input = find(strcmp(varargin, 'Tmodel'));if ~isempty(Input), obj.Tmodel = varargin{Input+1}; end
            end
        end
    end
    methods % Special Setting and Getting
        function out = get.SurfaceClass(obj)
                 out = class(obj.Target);
        end
        function points = get.Vertices(obj)
                 switch obj.SurfaceClass
                     case {'meshObj' 'LMObj'}
                         points = obj.Surface.Vertices;
                     case 'double'
                         points = obj.Surface;
                 end
        end
        function obj = set.Vertices(obj,points)
                 if ~size(points,1)==3, points = points';end
                 switch obj.SurfaceClass
                     case {'meshObj' 'LMObj'}
                         obj.Surface.Vertices = points;
                     case 'double'
                         obj.Surface = points;
                 end
        end
        function nr = get.nrV(obj)
            nr = size(obj.Vertices,2);
        end
        function s = get.sumB(obj)
            s = sum(obj.B);
        end
    end
    methods % InterFace functions
        function delete(obj) 
           if ~isempty(obj.Tmodel),delete(obj.Tmodel);end
        end
        function struc = saveobj(obj)
            struc = obj2struc(obj);
        end
        function struc = obj2struc(obj)
            % converting relevant information
            struc.Type = obj.Type;
            struc.SurfaceClass = obj.SurfaceClass;
            switch obj.SurfaceClass
                case {'meshObj' 'LMObj'}
                    struc.Surface = obj2struc(obj.Surface);
                otherwise
                    struc.Surface = obj.Surface;
            end
            if isempty(obj.Tmodel), struc.Tmodel = {}; return; end
            struc.Tmodel = obj2struc(obj.Tmodel);
        end
        function obj = struc2obj(obj,struc)
            obj.Type = struc.Type;
            switch struc.SurfaceClass
                case {'meshObj' 'LMObj'}
                    obj.Surface = struc2obj(eval(struc.TargetClass),struc.Surface);
                otherwise
                    obj.Surface = struc.Surface;
            end
            if isempty(struc.Tmodel),obj.Tmodel = {}; return; end
            obj.Tmodel = struc2obj(eval(struc.Tmodel.Type),struc.Tmodel);
        end
        function cobj = clone(obj)
                 cobj = eval(class(obj));
                 switch obj.SurfaceClass
                     case {'meshObj' 'LMObj'}
                         cobj.Surface = clone(obj.Surface);
                     otherwise
                         cobj.Surface = obj.Surface;
                 end
                 if isempty(obj.Tmodel), return; end
                 for t=1:1:obj.nrTM
                     cobj.Tmodel{t} = clone(obj.Tmodel{t});
                 end
        end
        function varargout = initialize(obj)
            if nargout == 1
               obj = clone(obj);
               varargout{1} = obj;
            end
            obj.B = ones(1,obj.nrV);
            initialize(obj.Tmodel);
        end
    end
    methods (Static = true) % Static interface functions
        function obj = loadobj(struc)
            obj = struc2obj(eval(struc.Type),struc);
        end        
    end
end % classdef

