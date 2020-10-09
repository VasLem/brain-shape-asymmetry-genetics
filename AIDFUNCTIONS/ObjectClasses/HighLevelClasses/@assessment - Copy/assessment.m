classdef assessment < superClass
    properties
        Norm = [];
        Scan = [];
        Average = [];
        T = [];
        Dysmorphogram = [];
        DysmorphogramShell = [];
        PercDysmorph = [];
        Significance = 2;
        ThresholdMap = [];
        PercThresh = [];
        Threshold = 3;
        DistanceMap = [];
        Displacements = [];
        RMSE = [];
        Signed = false;
        VectorField = [];
        DistanceRange = [0 5];
        CompensateScale = false;
        OutliersOnly = false;
        OutliersBinary = false;
        OutliersThreshold = 0.8;
        Tag = 'Assessement';
        UpdateSuperimposition = true;
        Update = true;
        NoiseLevel;
    end
    methods %Constructor
        function obj = assessment(varargin)
            obj = obj@superClass(varargin{:});
        end %Constructor
    end  
    methods % Special Setting and Getting
        function out = get.Norm(obj)
            out = obj.Norm;
            if ~superClass.isH(out), out = []; end
        end
        function out = get.Scan(obj)
            out = obj.Scan;
            if ~superClass.isH(out), out = []; end
        end
        function obj = set.Scan(obj,in)
            if ~(obj.Scan==in), obj.UpdateSuperimposition = true; end
            obj.Scan = in;
        end
        function obj = set.Norm(obj,in)
            if ~(obj.Norm==in), obj.UpdateSuperimposition = true; end
            obj.Norm = clone(in);
        end
        function obj = set.Significance(obj,in)
            if ~(obj.Significance==in),obj.UpdateSuperimposition = true; end
            obj.Significance = in;
        end
        function obj = set.CompensateScale(obj,in)
            if ~(obj.CompensateScale==in),obj.UpdateSuperimposition = true; end
            obj.CompensateScale = in;
        end
        function obj = set.Signed(obj,in)
            obj.Signed = in;
            %updateDistanceMap(obj);
        end
        function obj = set.UpdateSuperimposition(obj,in)
            obj.UpdateSuperimposition = in;
            if in == true, obj.Update = true; end
        end
        function obj = set.Threshold(obj,in)
            if ~(numel(obj.Threshold)==numel(in))||~(obj.Threshold==in),obj.Update = true; end
            obj.Threshold= in;
        end
        function obj = set.OutliersOnly(obj,in)
            if ~(obj.OutliersOnly==in),obj.Update = true; end
            obj.OutliersOnly= in;
        end
        function obj = set.OutliersBinary(obj,in)
            if ~(obj.OutliersBinary==in),obj.Update = true; end
            obj.OutliersBinary= in;
        end
        function obj = set.OutliersThreshold(obj,in)
            if ~(obj.OutliersThreshold==in),obj.Update = true; end
            obj.OutliersThreshold= in;
        end
    end
    methods % Interface Functions
        function out = update(obj,norm,scan)
            if nargout == 1, out = obj; end
            if nargin>1,obj.Norm = norm;end
            if nargin>2,obj.Scan = scan;end
            if isempty(obj.Norm), msgbox('Norm is empty, returning from action'); return;end
            if isempty(obj.Scan), msgbox('Scan is empty, returning from action'); return; end
            if ~obj.Update, return; end
            updateSuperimposition(obj);
            updateDysmorphogram(obj);
            updateVectorField(obj);
            updateDistanceMap(obj);
            updateThresholdMap(obj);
            updateAverage(obj);
            obj.Update = false;
        end
        function updateSuperimposition(obj)     
            if ~obj.UpdateSuperimposition, return; end
            if isempty(obj.Norm), msgbox('Norm is empty, returning from action'); return;end
            if isempty(obj.Scan), msgbox('Scan is empty, returning from action'); return; end
            % Pose alignement
            old = clone(obj.Norm);
            if obj.CompensateScale
               alignPoseAndScale(obj.Norm,obj.Scan,obj.Significance);
               obj.T = scaledRigidTM;              
            else
               alignPose(obj.Norm,obj.Scan,obj.Significance);
               obj.T = rigidTM;
            end
            obj.NoiseLevel = obj.Norm.NoiseLevel;
            match(obj.T,obj.Norm,old);
            delete(old);
            obj.UpdateSuperimposition = false;
        end
        function updateDysmorphogram(obj)
            if obj.UpdateSuperimposition, updateSuperimposition(obj); end
            obj.Dysmorphogram = 1-obj.Norm.Value;
            if obj.OutliersBinary
               values = zeros(1,length(obj.Dysmorphogram)); 
               values(obj.Dysmorphogram>=obj.OutliersThreshold) = 1;
               obj.Dysmorphogram = values;
               clear values;
            end
            total = length(obj.Dysmorphogram);
            outliers = sum(obj.Dysmorphogram);
            obj.PercDysmorph = (outliers/total)*100;
            clear outliers total;
        end
        function updateDistanceMap(obj)
            % Distance Map
            if obj.UpdateSuperimposition, updateSuperimposition(obj); end
            obj.Displacements = obj.Norm.Vertices-obj.Scan.Vertices;
            obj.DistanceMap = sqrt(sum(obj.Displacements.^2));           
            if obj.OutliersOnly
               obj.DistanceMap = obj.DistanceMap.*(obj.Dysmorphogram);
            end
            if obj.Signed
               if isempty(obj.VectorField), updateVectorField(obj); end
               angle = vectorAngle(obj.Scan.Gradient,obj.VectorField.Direction);
               signs = ones(1,obj.Scan.nrV);
               signs(find(angle>90)) = -1; %#ok<FNDSB>
               obj.DistanceMap = signs.*obj.DistanceMap; 
            end
            obj.RMSE = sqrt(mean(obj.DistanceMap.^2));       
        end
        function updateThresholdMap(obj)
            % ThresholdMap
            if obj.UpdateSuperimposition, updateSuperimposition(obj); end
            obj.ThresholdMap = zeros(size(obj.Dysmorphogram));
            obj.ThresholdMap(find(obj.DistanceMap>=obj.Threshold)) = 1;
            total = length(obj.ThresholdMap);
            outliers = sum(obj.ThresholdMap);
            obj.PercThresh = (outliers/total)*100;
        end
        function updateVectorField(obj)
            % VectorField
            if obj.UpdateSuperimposition, updateSuperimposition(obj); end
            obj.VectorField = vectorField3D;
            obj.VectorField.StartPoints = obj.Scan.Vertices;
            obj.VectorField.EndPoints = obj.Norm.Vertices;
%             if obj.Signed
%                angle = vectorAngle(obj.Scan.Gradient,obj.VectorField.Direction);
%                signs = ones(1,obj.Scan.nrV);
%                signs(find(angle>90)) = -1;
%                obj.VectorField.Signs = signs;
%             end
        end
        function updateAverage(obj)
            obj.Average = clone(obj.Scan);
            obj.Average.Vertices = (obj.Scan.Vertices+obj.Norm.Vertices)/2;
        end
        function showDysmorphogram(obj,v)
            if isempty(obj.Dysmorphogram), return; end
            if isempty(obj.Scan), return; end
            shell = clone(obj.Scan);
            shell.Value = obj.Dysmorphogram;
            if nargin < 2, v = viewer3DObj;end
            viewer(shell,'Viewer',v);
            v.Tag = 'Outlier Map';
            shell.ColorMode = 'Indexed';
            colormap(v.RenderAxes, jet);
            set(v.RenderAxes,'clim',[0 1]);
            if strcmp(shell.Type,'LMObj'), return; end
            delete(shell.PoseLM);
            if ~isempty(shell.Border), shell.Border.Visible = false; end
        end
        function showThresholdMap(obj,v)
            if isempty(obj.ThresholdMap), return; end
            if isempty(obj.Scan), return; end
            shell = clone(obj.Scan);
            shell.Value = obj.ThresholdMap;
            if nargin < 2, v = viewer3DObj;end
            viewer(shell,'Viewer',v);
            v.Tag = 'Threshold Map';
            shell.ColorMode = 'Indexed';
            colormap(v.RenderAxes, winter);
            set(v.RenderAxes,'clim',[0 1]);
            if strcmp(shell.Type,'LMObj'), return; end
            delete(shell.PoseLM);
            if ~isempty(shell.Border), shell.Border.Visible = false; end  
        end
        function showDistanceMap(obj,v)
            if isempty(obj.DistanceMap), return; end
            if isempty(obj.Scan), return; end
            shell = clone(obj.Scan);
            shell.Value = obj.DistanceMap;
            if nargin < 2, v = viewer3DObj;end
            viewer(shell,'Viewer',v);
            v.Tag = 'Distance Map';
            shell.ColorMode = 'Indexed';
            if ~obj.Signed
                set(v.RenderAxes,'clim',obj.DistanceRange);
            else
                set(v.RenderAxes,'clim',[-1*obj.DistanceRange(2) obj.DistanceRange(2)]);
            end
            if strcmp(shell.Type,'LMObj'), return; end
            delete(shell.PoseLM);
            if ~isempty(shell.Border), shell.Border.Visible = false; end
        end
        function showDisplacements(obj,v)
            if isempty(obj.Displacements), return; end
            if isempty(obj.Scan), return; end
            for i=1:1:3
                shell = clone(obj.Scan);
                shell.Value = obj.Displacements(i,:);
                if nargin < 2, v = viewer3DObj;end
                viewer(shell,'Viewer',v);
                switch i
                    case 1
                        v.Tag = 'X Disp.';
                    case 2
                        v.Tag = 'Y Disp.';
                    case 3
                        v.Tag = 'Z Disp.';
                end
                shell.ColorMode = 'Indexed';
                delete(shell.PoseLM);
                if ~isempty(shell.Border), shell.Border.Visible = false; end
            end
        end
        function showVectorField(obj,v)
            if isempty(obj.VectorField), return; end
            if isempty(obj.Scan), return; end
            if nargin < 2, v = viewer3DObj;end
            VF = clone(obj.VectorField);
            viewer(VF,'Viewer',v);
            v.Tag = 'Vector Field';
            shell = clone(obj.Scan);
            shell.Axes = v.RenderAxes;
            shell.SingleColor = [0.5 0.5 0.5];
            if strcmp(shell.Type,'LMObj')
                shell.Visible = true;
                shell.Value = obj.DistanceMap;
                shell.ColorMode = 'Single';
                shell2 = clone(obj.Norm);
                shell2.Axes = v.RenderAxes;
                shell2.Value = obj.DistanceMap;
                shell2.ColorMode = 'Indexed';
                shell2.Visible = true;
            else    
                delete(shell.PoseLM);
                reduceTriangles(shell,0.2);
                shell.ViewMode = 'Solid';
                shell.Visible = true;
                if ~isempty(shell.Border), shell.Border.Visible = false; end
                v.Renderer = 'opengl';
                shell.Alpha = 0.6;
                v.SceneLightVisible = true;
                v.SceneLightLinked = true;
                shell.Material = 'Dull';
            end
            set(v.RenderAxes,'clim',obj.DistanceRange);
        end
        function showScan(obj,v)
            if nargin < 2, v = viewer3DObj;end
            if isempty(obj.Scan), return; end
            viewer(obj.Scan,'Viewer',v);
            v.Tag = 'Scan';
            if strcmp(obj.Scan.Type,'LMObj'), return; end
            if ~isempty(obj.Scan.PoseLM), obj.Scan.PoseLM.Visible = false;end
            if ~isempty(obj.Scan.Border), obj.Scan.Border.Visible = false; end
        end
        function showNorm(obj,v)
            if nargin < 2, v = viewer3DObj;end
            if isempty(obj.Norm), return; end
            viewer(obj.Norm,'Viewer',v);
            v.Tag = 'Norm';
            if strcmp(obj.Norm.Type,'LMObj'), return; end
            if ~isempty(obj.Norm.PoseLM), obj.Norm.PoseLM.Visible = false; end
            if ~isempty(obj.Norm.Border), obj.Norm.Border.Visible = false; end
        end  
        function superimpose(obj)
            cNorm = clone(obj.Norm);
            cScan = clone(obj.Scan);
            v = viewer(cNorm);
            cNorm.SingleColor = [0.2039 0.3020 0.4941];
            cNorm.ViewMode = 'Points';
            cNorm.ColorMode = 'Single';
            v.Tag = 'Superimposition';
            cScan.Axes = v.RenderAxes;
            cScan.SingleColor = [0 0 1];
            cScan.Visible = true;
        end
        function show(obj)
            showScan(obj);
            showNorm(obj);
            showDysmorphogram(obj);
            showDistanceMap(obj);
            showThresholdMap(obj);
            showVectorField(obj);
            %superimpose(PreVsNe);
        end
    end
end % classdef