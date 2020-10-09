classdef textureMapPCA < scanFieldPCA
    properties
        Mask = [];
    end
    properties(Dependent = true)
        AvgTexture;
        AvgVec;
        RefScan;
        nrP;% Number of Pixels in the Mask
        UV; % UV coordinated in TextureMap
        Pindex;% Index of Pixels belonging to the Mask if present
        Psize; % Matrix Size of texture Map
        Psub;% Matrix coordinates of Pixels belonging to the white area
    end
    methods % Constructor
        function obj = textureMapPCA(varargin)
            obj = obj@scanFieldPCA(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.Mask(obj)
            out = obj.Mask;
            if ~mySuperClass.isH(out), out = []; return; end
        end
        function obj = set.Mask(obj,in)
            if isempty(in), obj.Mask= []; return; end
            switch class(in)
                case 'image2D'
                     if ~isempty(obj.Mask)&&~(in==obj.Mask), delete(obj.Mask); end
                     obj.Mask = clone(in);
                case {'uint8' 'double'}
                     if isempty(obj.Mask), obj.Mask = image2D; end
                     obj.Mask.Image = in;
                case 'meshObj'
                     if isempty(obj.Mask), obj.Mask = image2D; end
                     obj.Mask.Image = in.TextureMap.Image;
                otherwise
                    error('Wrong mask input');
            end
        end
        function obj = set.RefScan(obj,in)
            checkAverage(obj);
            if isempty(obj.Average.TextureMap), obj.Average.TextureMap = clone(in.TextureMap);end
            obj.Average.Vertices = in.Vertices;
            switch class(in)
                case 'LMobj'
                    return;
                case 'meshObj'
                    obj.Average.Faces = in.Faces;
                    obj.Average.UV = in.UV;
            end
        end
        function out = get.RefScan(obj)
            out = obj.Average;
        end
        function out = get.nrP(obj)
            out = [];
            if ~isempty(obj.Mask)
               out = numel(obj.Mask);
               return;
            end
            if ~isempty(obj.Average)
               out = numel(obj.Average.TextureMap);
               return;
            end
        end  
        function out = get.Pindex(obj)
            out = [];
            if ~isempty(obj.Mask)
               out = find(obj.Mask.Image);
               return; 
            end
            if ~isempty(obj.Average)
               out = (1:obj.nrP);
               return;
            end           
        end
        function out = get.UV(obj)
           if isempty(obj.Average), out = []; return; end
           out = obj.Average.UV; 
        end
        function out = get.Psub(obj)
            if isempty(obj.Pindex)||isempty(obj.Psize), out = []; return; end
            [I,J,K] = ind2sub(obj.Psize,obj.Pindex);
            out = [I;J;K];
        end
        function out = get.Psize(obj)
            out = [];
            if ~isempty(obj.Mask)
               out = size(obj.Mask.Image);
               return; 
            end
            if ~isempty(obj.Average)
               out = size(obj.Average.TextureMap.Image);
            end
        end
        function out = get.AvgVec(obj)
            if isempty(obj.Average), out = []; return; end
            out = obj.Average.TextureMap.Image(:);
            out = out(obj.Pindex)*255;
        end
        function obj = set.AvgVec(obj,in)
            im = Vec2Struc(obj,in);
            checkAverage(obj);
            obj.Average.TextureMap.Image = im;
        end
        function out = get.AvgTexture(obj)
            if isempty(obj.Average),out = []; return; end
            out = obj.Average.TextureMap;
        end
        function obj = set.AvgTexture(obj,in)
            checkAverage(obj);
            obj.Average.TextureMap = in;
        end
    end
    methods % Specific Interface Functions
        function out = Struc2Vec(obj,in) %#ok<INUSL>
           switch class(in)
                  case {'meshObj' 'LMObj'}
                       out = in.TextureMap.Image;
                  case 'double'
                       out = in*255;
                  case 'uint8'
                       out = double(in);
                  case 'image2D'
                       out = in.Image*255;
                  otherwise
                       error('wrong input for input struc2vec texture map model');
           end
           out = out(obj.Pindex);
        end
        function out = Vec2Struc(obj,in) %#ok<INUSL>
                 if ~isempty(obj.Mask)
                    out = obj.Mask.Image;
                 elseif ~isempty(obj.Average);
                    out = obj.Average.TextureMap.Image;
                 else
                    error('Not able to convert vector, lack of reference data (Mask or Average or Reference Scan');
                 end
                 out(obj.Pindex) = in/255;
                 if nargout < 2, return; end
                 objout = embedInScan(obj,out);
        end
        function out = IndexStruc2Vec(obj,in) %#ok<INUSL>
            % To be improved!
            out = intersect(obj.Pindex,in);
        end
        function out = IndexVec2Struc(obj,in) %#ok<INUSL>
            % To be improved
            out = in;
        end
        function out = WeightStruc2Vec(obj,in) %#ok<INUSL>
            % To be improved
            switch class(in)
                case 'image2D'
                    out = in.Image;
                case 'double'
                    out = in;
            end
            out = out(obj.Pindex);
        end
        function out = WeightVec2Struc(obj,in) %#ok<INUSL>
            % To be improved
            out = Vec2Struc(obj,in);
        end
        function out = getData(obj,in) %#ok<INUSL>
            switch class(in)
                case 'double'
                    out = in;
                case 'struct'
                    if ~isfield(in,'TextureMap'), error('TextureMap Field is missing in data structure'); end
                    out = in.TextureMap;
            end
            out = out(obj.Pindex,:);
        end
        function [f,scan] = initializeShowPC(obj,coeff)
%                  f = viewer3DObj;
%                  scan.Axes = f.RenderAxes;scan.Visible = true;scan.Selected = true;
%                  scan.ColorMode = 'Texture';
%                  subdivideTriangles(scan,'runs',2);
%                  scan.Material = 'Dull';
%                  set(f.Toolbar.light_toggle,'State','off');
%                  set(f.Toolbar.link_toggle,'State','off');
                   f = figure;
                   scan = clone(obj.Average);
                   scan.TextureMap.Image = reconstruct(obj,coeff);
                   scan.UserData = f;
                   clf(f);
                   imshow(scan.TextureMap.Image);
        end
        function updateShowPC(obj,scan,coeff)
                 scan.TextureMap.Image = reconstruct(obj,coeff);
                 %scan.ColorMode = 'Texture';
                 clf(scan.UserData);
                 imshow(scan.TextureMap.Image);
        end
        function objout = embedInScan(obj,in)
                 if isempty(obj.Average), objout = image2D('Image',out); return; end
                 objout = clone(obj.Average);
                 objout.TextureMap.Image = in;
        end
    end
end
