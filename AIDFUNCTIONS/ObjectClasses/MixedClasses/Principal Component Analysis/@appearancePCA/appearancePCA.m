classdef appearancePCA < scanPCA
    properties
        Average = [];
        Shape = [];
        Texture = [];
    end
    properties(Dependent = true)
        AvgVertices;% The Average vertices
        AvgTexture;% the Average Texture
        AvgVec;% The average vertices in vector notation
        RefScan;% The reference scan defining the Faces
        nrSC;% Number of Shape Coefficients
        nrTC;% Number of Texture Coefficients
        nrC;% Total number of Coefficients 
        WS;% Shape versus Texture Normalisation
    end
    methods % Constructor
        function obj = appearancePCA(varargin)
            obj = obj@scanPCA(varargin{:});         
        end
    end
    methods % Special Setting & Getting
        function out = get.Average(obj)
            out = obj.Average;
            if ~superClass.isH(out), out = []; return; end
            if ~isempty(obj.Shape)
               checkAverage(obj.Shape);
               out.Vertices = obj.Shape.Average.Vertices;
               out.Faces = obj.Shape.Average.Faces;
            end
            if ~isempty(obj.Texture)
               checkAverage(obj.Texture)
               switch class(obj.Texture)
                   case 'texturePCA'
                       out.TextureColor = obj.Texture.Average.TextureColor;
                   case 'textureMapPCA'
                       out.TextureMap = obj.Texture.Average.TextureMap.Image;
                       out.UV = obj.Texture.Average.UV;
               end
            end
        end
        function obj = set.Average(obj,in)
           if ~isempty(obj.Average)&&~(in==obj.Average), delete(obj.Average); end
           obj.Average = in; 
        end
        function out = get.Shape(obj)
            out = obj.Shape;
            if ~mySuperClass.isH(out), out = []; end
        end
        function obj = set.Shape(obj,in)
            if ~isempty(obj.Shape)&&~(in==obj.Shape), delete(obj.Shape); end
            obj.Shape = in;
            checkAverage(obj);
        end
        function out = get.AvgVertices(obj)
            if isempty(obj.Shape), out = []; return; end
            out = obj.Shape.AvgVertices;
        end
        function obj = set.AvgVertices(obj,in)
            if isempty(obj.Shape), return; end
            obj.Shape.AvgVertices = in;
        end
        function out = get.Texture(obj)
            out = obj.Texture;
            if ~mySuperClass.isH(out), out = []; end
        end
        function obj = set.Texture(obj,in)
            if ~isempty(obj.Texture)&&~(in==obj.Texture), delete(obj.Texture); end
            obj.Texture = in;
            checkAverage(obj);
            obj.Average.ColorMode = 'Texture';
        end
        function out = get.AvgTexture(obj)
            if isempty(obj.Texture), out = []; return; end
            out = obj.Texture.AvgTexture;
        end
        function obj = set.AvgTexture(obj,in)
            if isempty(obj.Texture), return; end
            obj.Texture.AvgTexture = in;
        end        
        function obj = set.RefScan(obj,in)
            checkShape(obj);obj.Shape.RefScan = in;
            checkTexture(obj);obj.Texture.RefScan = in;
        end
        function out = get.RefScan(obj)
            out = obj.Average;
        end
        function out = get.AvgVec(obj)
            out = zeros(obj.nrC,1);
        end
        function obj = set.AvgVec(obj,in)
            %dummy does not change
        end
        function out = get.nrSC(obj)
            if isempty(obj.Shape), out  =[]; return;end
            out = obj.Shape.nrEV;
        end
        function out = get.nrTC(obj)
            if isempty(obj.Texture), out  =[]; return;end
            out = obj.Texture.nrEV;
        end
        function out = get.nrC(obj)
            out = obj.nrSC+obj.nrTC;
        end
        function out = get.WS(obj)
            if isempty(obj.Shape)||isempty(obj.Texture), out = 1; return; end
            % the ratio of complete variances
            %out = 1;
            out = (sum(obj.Texture.EigVal)/sum(obj.Shape.EigVal));
        end
    end
    methods % Specific Interface Functions
        function checkShape(obj)
            if isempty(obj.Shape), obj.Shape = shapePCA; end
        end
        function checkTexture(obj,type)   
            if ~isempty(obj.Texture), return; end
            if nargin < 2, type = 'texturePCA'; end% default texture as rgb per vertex
            obj.Texture = eval(type);
        end
        function out = Coeff2Vec(obj,in)
            % Turn given PCA coeff into a Vector
            if size(in,1)==1, in = in'; end
            out = obj.EigVec*in;
            if obj.nrSC > 0, out(1:obj.nrSC) = out(1:obj.nrSC)/obj.WS;end
        end
        function out = Struc2Vec(obj,in) %#ok<INUSL>
           out = [];
           if ~isempty(obj.Shape), out = [out Struc2Coeff(obj.Shape,in)];end
           if ~isempty(obj.Texture), out = [out Struc2Coeff(obj.Texture,in)];end
        end
        function out = Vec2Struc(obj,in) %#ok<INUSL>
                 if obj.nrSC>0
                    shapein = in(1:obj.nrSC);
                    out.Vertices = Coeff2Struc(obj.Shape,shapein);
                 end
                 if obj.nrTC>0
                    texturein = in(end-obj.nrTC+1:end);
                    switch class(obj.Texture)
                        case 'texturePCA'
                            out.TextureColor = Coeff2Struc(obj.Texture,texturein);
                        case 'textureMapPCA'
                            out.TextureMap = Coeff2Struc(obj.Texture,texturein);
                    end
                 end            
        end
        function out = IndexStruc2Vec(obj,in) %#ok<INUSL>
            % To be improved, but probably fine like this
            out = in;
        end
        function out = IndexVec2Struc(obj,in) %#ok<INUSL>
            % To be improved, but probably fine like this
            out = in;
        end
        function out = WeightStruc2Vec(obj,in) %#ok<INUSL>
            % To be improved, but probably fine like this
            out = in;
        end
        function out = WeightVec2Struc(obj,in) %#ok<INUSL>
            % To be improved, but probably fine like this
            out = in;
        end
        function out = getData(obj,in) %#ok<INUSL>
           out = in;
        end
        function [f,scan] = initializeShowPC(obj,coeff)
                 f = viewer3DObj;
                 scan = clone(obj.Average);
                 tmp = reconstruct(obj,coeff);
                 scan.Vertices = tmp.Vertices;
                 switch class(obj.Texture)
                     case 'texturePCA'
                         scan.TextureColor = tmp.TextureColor;
                     case 'textureMapPCA'
                         scan.TextureMap = tmp.TextureMap;
                 end
                 scan.Axes = f.RenderAxes;scan.Visible = true;scan.Selected = true;
                 scan.SingleColor = [0.8 0.8 0.8];scan.ColorMode = 'Texture';
                 scan.Material = 'Dull';
                 set(f.Toolbar.light_toggle,'State','on');
                 set(f.Toolbar.link_toggle,'State','on');
        end
        function updateShowPC(obj,scan,coeff)
                 tmp = reconstruct(obj,coeff);
                 scan.Vertices = tmp.Vertices;
                 switch class(obj.Texture)
                     case 'texturePCA'
                         scan.TextureColor = tmp.TextureColor;
                     case 'textureMapPCA'
                         scan.TextureMap = tmp.TextureMap;
                 end         
        end
        function objout = embedInScan(obj,in)
                 checkAverage(obj);
                 objout = clone(obj.Average);
                 struc2obj(objout,in);
        end
        function out = getShapeCoeff(obj,in)
            out = getVec(obj,in);
            out = out(1:obj.nrSC);
        end
        function out = getTextureCoeff(obj,in)
            out = getVec(obj,in);
            out = out(obj.nrSC+1:end);
        end
    end
end