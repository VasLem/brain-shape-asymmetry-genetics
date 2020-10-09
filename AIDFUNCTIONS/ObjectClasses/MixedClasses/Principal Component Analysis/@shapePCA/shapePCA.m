classdef shapePCA < scanFieldPCA
    properties(Dependent = true)
        AvgVertices;% The average vertices 
        AvgVec;% The average vertices in vector notation
        %RefScan;% The reference scan defining the Faces
    end
    properties 
        RefScan = [];
    end
    methods % Constructor
        function obj = shapePCA(varargin)
            obj = obj@scanFieldPCA(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function obj = set.RefScan(obj,in)
            checkAverage(obj);
            if isempty(obj.Average.Vertices), obj.Average.Vertices = in.Vertices; end
            switch class(in)
                case 'LMObj'
                    return;
                case 'meshObj'
                    obj.Average.Faces = in.Faces;
            end
            obj.RefScan = in;
        end
        function out = get.RefScan(obj)
            %out = obj.Average;
            out = obj.RefScan;
        end
        function out = get.AvgVec(obj)
            if isempty(obj.Average), out = []; return; end
            out = obj.Average.Vertices(:);
%             if obj.Centering
%                 out = obj.Average.Vertices(:);
%             else
%                 out = obj.RefScan.Vertices(:);
%             end
        end
        function obj = set.AvgVec(obj,in)
            vert = Vec2Struc(obj,in);
            checkAverage(obj);
            if (size(vert,2)~=obj.Average.nrV); obj.Average = meshObj; end
            obj.Average.Vertices = vert;
        end
        function out = get.AvgVertices(obj)
            if isempty(obj.Average), out = []; return;end
            out = obj.Average.Vertices;
%             if obj.Centering
%                out = obj.Average.Vertices;
%             else
%                out = obj.RefScan.Vertices;
%             end
        end
        function obj = set.AvgVertices(obj,in)
            checkAverage(obj);
            obj.Average.Vertices = in;
        end
    end
    methods % Specific Interface Functions
        function out = Struc2Vec(obj,in) %#ok<INUSL>
           switch class(in)
                  case {'meshObj' 'LMObj'}
                       out = in.Vertices(:);
                  case {'double' 'single'}
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
                case {'double' 'single'}
                    out = in;
                case 'struct'
                    if ~isfield(in,'Shape'), error('Shape Field is missing in data structure'); end
                    out = in.Shape;
            end
        end
        function [f,scan] = initializeShowPC(obj,coeff)
                 f = viewer3DObj;
                 scan = clone(obj.RefScan);
                 if ~isempty(scan.PoseLM), delete(scan.PoseLM); end
                 if obj.Centering
                     scan.Vertices = reconstruct(obj,coeff);
                 else
                    scan.Vertices = reconstruct(obj,coeff) + obj.RefScan.Vertices;
                 end
                 scan.Axes = f.RenderAxes;scan.Visible = true;scan.Selected = true;
                 scan.SingleColor = [0.8 0.8 0.8];scan.ColorMode = 'Single';
                 scan.ViewMode = 'Solid';
                 scan.Material = 'Dull';
                 %disp('initializing');
                 f.SceneLightVisible = true;f.SceneLightLinked = true;
                 set(f.Toolbar.light_toggle,'State','on');
                 set(f.Toolbar.link_toggle,'State','on');
        end
        function updateShowPC(obj,scan,coeff)
                 if obj.Centering
                    scan.Vertices = reconstruct(obj,coeff);           
                 else
                    scan.Vertices = reconstruct(obj,coeff) + obj.RefScan.Vertices; 
                 end
        end
        function objout = embedInScan(obj,in)
                 checkAverage(obj);
                 %objout = meshObj('Vertices',in,'Faces',obj.Average.Faces);
                 if ~obj.Centering, in = in + obj.RefScan.Vertices; end
                 objout = meshObj('Vertices',in,'Faces',obj.RefScan.Faces);
        end
    end
end