classdef ShapeValuePCA < scanPCA
    properties
        Average = [];
        Shape = [];
        Value = [];
        Range = [0 1];
        Normalize = true;
    end
    properties(Dependent = true)
        AvgVertices;% The Average vertices
        AvgValue;% the Average Texture
        AvgVec;% The average vertices in vector notation
        RefScan;% The reference scan defining the Faces
        nrSC;% Number of Shape Coefficients
        nrVC;% Number of Texture Coefficients
        nrC;% Total number of Coefficients 
        WS;% Shape versus Texture Normalisation
    end
    methods % Constructor
        function obj = ShapeValuePCA(varargin)
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
            if ~isempty(obj.Value)
               checkAverage(obj.Value)
               out.Value = obj.Value.Average.Value;
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
        function out = get.Value(obj)
            out = obj.Value;
            if ~mySuperClass.isH(out), out = []; end
        end
        function obj = set.Value(obj,in)
            if ~isempty(obj.Value)&&~(in==obj.Value), delete(obj.Value); end
            obj.Value = in;
            checkAverage(obj);
            obj.Average.ColorMode = 'Indexed';
        end
        function out = get.AvgValue(obj)
            if isempty(obj.Value), out = []; return; end
            out = obj.Value.AvgValue;
        end
        function obj = set.AvgValue(obj,in)
            if isempty(obj.Value), return; end
            obj.Value.AvgValue = in;
        end        
        function obj = set.RefScan(obj,in)
            checkShape(obj);obj.Shape.RefScan = in;
            checkValue(obj);obj.Value.RefScan = in;
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
        function out = get.nrVC(obj)
            if isempty(obj.Value), out  =[]; return;end
            out = obj.Value.nrEV;
        end
        function out = get.nrC(obj)
            out = obj.nrSC+obj.nrVC;
        end
        function out = get.WS(obj)
            if isempty(obj.Shape)||isempty(obj.Value), out = 1; return; end
            % the ratio of complete variances
            if obj.Normalize
                out = sum(obj.Value.EigVal)/sum(obj.Shape.EigVal);
            else
                out = 1;
            end
        end
    end
    methods % Specific Interface Functions
        function checkShape(obj)
            if isempty(obj.Shape), obj.Shape = shapePCA; end
        end
        function checkValue(obj,type)   
            if ~isempty(obj.Value), return; end
            if nargin < 2, type = 'valuePCA'; end% default texture as rgb per vertex
            obj.Value = eval(type);
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
           if ~isempty(obj.Value), out = [out Struc2Coeff(obj.Value,in)];end
        end
        function out = Vec2Struc(obj,in) %#ok<INUSL>
                 if obj.nrSC>0
                    shapein = in(1:obj.nrSC);
                    out.Vertices = Coeff2Struc(obj.Shape,shapein);
                 end
                 if obj.nrVC>0
                    valuein = in(end-obj.nrVC+1:end);
                    out.Value = Coeff2Struc(obj.Value,valuein);
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
                 scan.Value = tmp.Value;
                 scan.Axes = f.RenderAxes;scan.Visible = true;scan.Selected = true;
                 scan.SingleColor = [0.8 0.8 0.8];scan.ColorMode = 'Texture';
                 scan.Material = 'Dull';
                 set(f.Toolbar.light_toggle,'State','on');
                 set(f.Toolbar.link_toggle,'State','on');
                 set(f.RenderAxes,'clim',obj.Range);
        end
        function updateShowPC(obj,scan,coeff)
                 tmp = reconstruct(obj,coeff);
                 scan.Vertices = tmp.Vertices;
                 scan.Value = tmp.Value;       
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
        function out = getValueCoeff(obj,in)
            out = getVec(obj,in);
            out = out(obj.nrSC+1:end);
        end
    end
end