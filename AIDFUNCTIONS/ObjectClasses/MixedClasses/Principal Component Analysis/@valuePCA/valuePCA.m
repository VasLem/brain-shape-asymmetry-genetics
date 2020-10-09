classdef valuePCA < scanFieldPCA
    properties(Dependent = true)
        AvgValue; % Average textureColor values
        AvgVec;
        RefScan;
    end
    properties 
        Range = [0 1];
    end
    methods % Constructor
        function obj = valuePCA(varargin)
            obj = obj@scanFieldPCA(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function obj = set.RefScan(obj,in)
            checkAverage(obj);
            if isempty(obj.Average.Value), obj.Average.Value = in.Value;end
            obj.Average.Vertices = in.Vertices;
            switch class(in)
                case 'LMObj'
                    return;
                case 'meshObj'
                    obj.Average.Faces = in.Faces;
            end
            obj.Average.ColorMode = 'Indexed';
        end
        function out = get.RefScan(obj)
            out = obj.Average;
        end
        function out = get.AvgVec(obj)
            if isempty(obj.Average), out = []; return; end
            out = obj.Average.Value(:);
        end
        function obj = set.AvgVec(obj,in)
            in = Vec2Struc(obj,in);
            checkAverage(obj);
            obj.Average.Value = in;
        end
        function out = get.AvgValue(obj)
            if isempty(obj.Average), out = []; return;end
            out = obj.Average.Value;
        end
        function obj = set.AvgValue(obj,in)
            checkAverage(obj);
            obj.Average.Value = in;
        end
    end
    methods % Specific Interface Functions
        function out = Struc2Vec(obj,in) %#ok<INUSL>
           switch class(in)
                  case {'meshObj' 'LMObj'}
                       out = in.Value(:);
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
            out = in;
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
                    if ~isfield(in,'Value'), error('Value Field is missing in data structure'); end
                    out = in.Value;
            end
        end
        function [f,scan] = initializeShowPC(obj,coeff)
                 f = viewer3DObj;
                 scan = clone(obj.Average);
                 scan.Value = reconstruct(obj,coeff);
                 scan.Axes = f.RenderAxes;scan.Visible = true;scan.Selected = true;
                 scan.ColorMode = 'Indexed';
                 scan.Material = 'Dull';
                 set(f.Toolbar.light_toggle,'State','off');
                 set(f.Toolbar.link_toggle,'State','off');
                 set(f.RenderAxes,'clim',obj.Range);
        end
        function updateShowPC(obj,scan,coeff)
                 scan.Value = reconstruct(obj,coeff);           
        end
        function objout = embedInScan(obj,in)
                 checkAverage(obj);
                 objout = meshObj('Vertices',obj.Average.Vertices,'Faces',obj.Average.Faces,'Value',in);
        end
    end
end