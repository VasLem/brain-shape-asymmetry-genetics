classdef texturePCA < scanFieldPCA
    properties(Dependent = true)
        AvgTexture; % Average textureColor values
        AvgVec;
        RefScan;
    end
    methods % Constructor
        function obj = texturePCA(varargin)
            obj = obj@scanFieldPCA(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function obj = set.RefScan(obj,in)
            checkAverage(obj);
            if isempty(obj.Average.TextureColor), obj.Average.TextureColor = in.TextureColor;end
            obj.Average.Vertices = in.Vertices;
            switch class(in)
                case 'LMObj'
                    return;
                case 'meshObj'
                    obj.Average.Faces = in.Faces;
            end
        end
        function out = get.RefScan(obj)
            out = obj.Average;
        end
        function out = get.AvgVec(obj)
            if isempty(obj.Average), out = []; return; end
            out = obj.Average.TextureColor(:)*255;
        end
        function obj = set.AvgVec(obj,in)
            in = Vec2Struc(obj,in);
            checkAverage(obj);
            if (size(in,2)~=size(obj.Average.TextureColor,2)),obj.Average = meshObj;end
            obj.Average.TextureColor = in;
        end
        function out = get.AvgTexture(obj)
            if isempty(obj.Average), out = []; return;end
            out = obj.Average.TextureColor;
        end
        function obj = set.AvgTexture(obj,in)
            checkAverage(obj);
            obj.Average.TextureColor = in;
        end
    end
    methods % Specific Interface Functions
        function out = Struc2Vec(obj,in) %#ok<INUSL>
           switch class(in)
                  case {'meshObj' 'LMObj'}
                       out = in.TextureColor(:);
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
            if nargout < 2, return; end
            objout = embedInScan(obj,out);
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
                    if ~isfield(in,'TextureColor'), error('TextureColor Field is missing in data structure'); end
                    out = in.TextureColor;
            end
        end
        function [f,scan] = initializeShowPC(obj,coeff)
                 f = viewer3DObj;
                 scan = clone(obj.Average);
                 scan.TextureColor = reconstruct(obj,coeff);
                 scan.Axes = f.RenderAxes;scan.Visible = true;scan.Selected = true;
                 scan.ColorMode = 'Texture';
                 scan.Material = 'Dull';
                 set(f.Toolbar.light_toggle,'State','off');
                 set(f.Toolbar.link_toggle,'State','off');
        end
        function updateShowPC(obj,scan,coeff)
                 scan.TextureColor = reconstruct(obj,coeff);           
        end
        function objout = embedInScan(obj,in)
                 checkAverage(obj);
                 objout = meshObj('Vertices',obj.Average.Vertices,'Faces',obj.Average.Faces,'TextureColor',in);
        end
    end
end