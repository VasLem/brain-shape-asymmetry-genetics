classdef shapeBasePCA < shapePCA
    properties
        ShapeBase = [];
    end
    methods % Constructor
        function obj = shapeBasePCA(varargin)
            obj = obj@shapePCA(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.ShapeBase(obj)
           out = obj.ShapeBase;
           if ~superClass.isH(out), out = [];return; end
        end
        function obj = set.ShapeBase(obj,in)
           if ~isempty(obj.ShapeBase)&&~(in==obj.ShapeBase), delete(obj.ShapeBase); end
           obj.ShapeBase = in; 
        end
    end
    methods % Specific Interface Functions
        function out = Struc2Vec(obj,in) %#ok<INUSL>
           switch class(in)
                  case {'meshObj' 'LMObj'}
                       out = in.Vertices(:);
                  case 'double'
                       if size(in,2)==3, in = in'; end
                       if ~size(in,1)==3, error('Input must contain 3D vertices'); end
                       out = in(:);
                  otherwise
                       error('wrong input for input struc2vec shape model');
           end
        end
        function out = Vec2Struc(obj,in) %#ok<INUSL>
            % converts a vector 2 structure
            out = reshape(in,3,length(in)/3);
        end
        function out = IndexStruc2Vec(obj,in) %#ok<INUSL>
            if isempty(obj.Average), error('average or refscan needs to be set first'); end
            tmp = zeros(3,obj.Average.nrV);
            tmp(:,in) = 1;tmp = tmp(:);
            out = find(tmp==1);
        end
        function out = IndexVec2Struc(obj,in)
            tmp = zeros(size(obj.AvgVec));
            tmp(in) = 1;
            tmp = reshape(tmp,3,length(tmp)/3);
            out = find(tmp(1,:));
        end
        function out = WeightStruc2Vec(obj,in) %#ok<INUSL>
            if size(in,1)==3, out = in(:); return; end
            in = repmat(in,3,1);
            out = in(:);
        end
        function out = WeightVec2Struc(obj,in) %#ok<INUSL>
            out = reshape(in,3,length(in)/3);
        end
        function out = getData(obj,in) %#ok<INUSL>
            switch class(in)
                case 'double'
                    out = in;
                case 'struct'
                    if ~isfield(in,'Shape'), error('Shape Field is missing in data structure'); end
                    out = in.Shape;
            end
        end
        function [f,scan] = initializeShowPC(obj,coeff)
                 f = viewer3DObj;
                 scan = clone(obj.Average);
                 scan.Vertices = reconstruct(obj,coeff);
                 scan.Axes = f.RenderAxes;scan.Visible = true;scan.Selected = true;
                 scan.SingleColor = [0.8 0.8 0.8];scan.ColorMode = 'Single';
                 scan.Material = 'Dull';
                 set(f.Toolbar.link_toggle,'State','on');
                 set(f.Toolbar.light_toggle,'State','on');
        end
        function updateShowPC(obj,scan,coeff)
                 scan.Vertices = reconstruct(obj,coeff);           
        end
        function objout = embedInScan(obj,in)
                 checkAverage(obj);
                 objout = meshObj('Vertices',in,'Faces',obj.Average.Faces);
        end
    end
end